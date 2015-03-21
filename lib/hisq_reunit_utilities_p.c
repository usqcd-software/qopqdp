/* ******************* hisq_reunit_utilities_p.c **********************
   Auxiliary routines needed for reunitarization 
   A. Bazavov, 02/2010
*/

#include <stdlib.h>
#include <math.h>
#include <qop_internal.h>

//AB TODO:
/* 1) rewrite all operations through QLA calls,
      do not access matrix elements directly, etc.
   2) count how many SVDs and force filters are applied
      (flags FFFILTER, LLSVD)
   3) make verbosity flag for outputting matrix state for SVD 
      and/or force filter: FFFILTEROUT, LLSVDOUT
   4) whether to use SVD and/or force filter, and also their 
      parameters should be options of QOP, not compile-time flags;
      overall we need the following parameters:
      unitarization group
      unitarization method
      two values for criteria for SVD
      one value for force filter
      flag if to use SVD
      flag if to use force filter
*/


#if 0
static void 
PrintMone(QLA_ColorMatrix m){
  printf("Field: %e %e %e\n", m.e[0][0].real,
	 m.e[0][1].real, m.e[0][2].real);
  printf("       %e %e %e\n", m.e[1][0].real,
	 m.e[1][1].real, m.e[1][2].real);
  printf("       %e %e %e\n\n", m.e[2][0].real,
	 m.e[2][1].real, m.e[2][2].real);
}

#endif

#if QOP_Precision == 'F'
#define QOP_U3_UNIT_ANALYTIC_EPS 1.0e-6
//#define SU3_ROOT_INV_NORM_EPS 1.0e-6
#else
#define QOP_U3_UNIT_ANALYTIC_EPS 1.0e-14
//#define SU3_ROOT_INV_NORM_EPS 1.0e-9
#endif
#define QOP_PI 3.14159265358979323846

// Forward declaration for SVD functions
static int QOP_D3_svd2x2bidiag(QOP_info_t *info, QLA_Real *a00, QLA_Real *a01, 
			       QLA_Real *a11, QLA_Real U2[2][2], QLA_Real V2[2][2]);
static int QOP_D3_svd3x3(QOP_info_t *info, QLA_ColorMatrix *A, QLA_Real *sigma, 
			 QLA_ColorMatrix *U, QLA_ColorMatrix *V);

#define Uelem(a,b) QLA_elem_M(*U,a,b)

#if 1
// Determinant of 3x3 complex matrix
QLA_Complex QOPPC(su3_mat_det)( QLA_ColorMatrix *U) {
  QLA_Complex a, b, m0, m1, m2, cdet;

  // brute-force calculation
  QLA_c_eq_c_times_c( a, Uelem(1,1), Uelem(2,2) );
  QLA_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,1) );
  QLA_c_eq_c_minus_c( m0, a, b );
  
  QLA_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,2) );
  QLA_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,0) );
  QLA_c_eq_c_minus_c( m1, a, b );
  
  QLA_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,1) );
  QLA_c_eq_c_times_c( b, Uelem(1,1), Uelem(2,0) );
  QLA_c_eq_c_minus_c( m2, a, b );

  QLA_c_eq_c_times_c( a, Uelem(0,0), m0 );
  QLA_c_eq_c_times_c( b, Uelem(0,1), m1 );
  QLA_c_eq_c_minus_c( cdet, a, b );
  QLA_c_eq_c_times_c( a, Uelem(0,2), m2 );
  QLA_c_eq_c_plus_c( cdet, cdet, a );
  
  
  return cdet;
}
#else
// Determinant of 3x3 complex matrix
QLA_Complex
QOPPC(su3_mat_det)( QLA_ColorMatrix *U)
{
  QLA_Complex d;
  QLA_C_eq_det_M(&d, U);
  return d;
}
#endif

// Analytic reunitarization
int
QOPPC(u3_un_analytic)(QOP_info_t *info, QLA_ColorMatrix *V, QLA_ColorMatrix *W)
{
  QLA_Real c0, c1, c2, S, S3, R, /*R2, CQ3,*/ RoS, theta, theta3, pi23, denom;
  QLA_Real g0, g1, g2, g0sq, g1sq, g2sq, f0, f1, f2, us, vs, ws, det_check=0;
  QLA_Real sigma[3];
  QLA_Complex det;
  QLA_ColorMatrix Q, Q2, Q3, Uleft, Vright, S1, S2;
  int i, j, perform_svd;
  int svd_calls = 0;
  size_t nflops = 0;

  // No SVD initially
  perform_svd=0;

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//  printf("Enter QOPPC(u3_un_analytic) in hisq_reunit_utilities_p.c\n");
#endif

  if(QOP_hisq_links.reunit_allow_svd) {
    // get determinant for future comparison
    det = QOPPC(su3_mat_det)( V );
    // det_check = |det|^2
    det_check = det.real*det.real + det.imag*det.imag;
  }

  nflops += 67;

  /* Hermitian matrix: Q=V^+ V */
  QLA_M_eq_Ma_times_M( &Q, V, V );
  nflops += 198;

  /* Q^2 */
  QLA_M_eq_M_times_M( &Q2, &Q, &Q );
  nflops += 198;

  /* Q^3 */
  QLA_M_eq_M_times_M( &Q3, &Q2, &Q );
  nflops += 198;

  /* (real) traces */
  QLA_R_eq_re_trace_M(&c0, &Q);
  QLA_R_eq_re_trace_M(&c1, &Q2);
  c1 /= 2;
  QLA_R_eq_re_trace_M(&c2, &Q3);
  c2 /= 3;

  S = c1/3 - c0 * (c0/18);
  nflops += 11;
  if( fabs(S)<QOP_U3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
    nflops += 3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    //R2 = R*R;
    //CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;
    nflops += 15;
    if( !( fabs(RoS)<1.0 ) ) {
      if( R>0 ) {
        theta = 0.0;
      }
      else {
        theta = QOP_PI;
      }
    } 
    else {
      theta = acos( RoS );
      if(isnan(theta)){
        printf("Hit NaN in QOPPC(u3_un_analytic)\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for( i=0; i<3; i++ ) {
          for( j=0; j<3; j++ ) {
            printf( "%24.18g %24.18g\n", QLA_real(QLA_elem_M(*V,i,j)),
		    QLA_imag(QLA_elem_M(*V,i,j)) );
          }
        }
        QDP_abort(0);
      }
    }

    /* eigenvalues of Q */
    theta3 = theta/3;
    pi23 = (2*QOP_PI) / 3;
    g0 = c0/3 + 2 * S * cos( theta3 );
    g1 = c0/3 + 2 * S * cos( theta3 + pi23 );
    g2 = c0/3 + 2 * S * cos( theta3 + 2*pi23 );

    nflops += 21;
  }

  if(QOP_hisq_links.reunit_allow_svd){
    
    if(!QOP_hisq_links.reunit_svd_only){
      /* conditions to call SVD */
      if(det_check!=0) {
	if( fabs(det_check-g0*g1*g2)/fabs(det_check)
	    > QOP_hisq_links.reunit_svd_rel_error) {
	  perform_svd=1;
	  nflops += 4;
	}
      }
      if( det_check<QOP_hisq_links.reunit_svd_abs_error ) {
	perform_svd=1;
      }
    } else {
      /* exclusively use SVD for finding eigenvalues and reunitarization,
	 this is slow since Q, Q^2 and Q^3 are calculated anyway;
	 this option is for testing: under normal circumstances SVD is
	 rarely used, which makes it harder to test, therefore one can
	 SVD with this switch */
      perform_svd=1;
    }
  } /* reunit_allow_svd */

  if( 0!=perform_svd ) {
    /* call SVD */
    QOPPC(svd3x3)( info, V, sigma, &Uleft, &Vright);
    
    svd_calls++;
    
    if(!QOP_hisq_links.reunit_svd_only && QOP_hisq_links.svd_values_info){
      printf("*** QOPQDP reunitarization while calculating fat links ***\n");
      printf("*** Resort to using SVD ***\n");
      printf("*** printing from node %d ***\n",QDP_this_node);
      printf("Eigenvalues from cubic equation:\n");
      printf( "  g0 = %28.18f\n", g0 );
      printf( "  g1 = %28.18f\n", g1 );
      printf( "  g2 = %28.18f\n", g2 );
      printf("Eigenvalues from singular value decomposition:\n");
      printf( "  g0 = %28.18f\n", sigma[0]*sigma[0] );
      printf( "  g1 = %28.18f\n", sigma[1]*sigma[1] );
      printf( "  g2 = %28.18f\n", sigma[2]*sigma[2] );
    }
    
    // W = Uleft * Vright^+
    QLA_M_eq_M_times_Ma( W, &Uleft, &Vright );
    nflops += 198;
    
  }  else { 
    // SVD is not performed, just use eigenvalues from cubic equation
    
    /* roots of eigenvalues */
    g0sq = sqrt( g0 );
    g1sq = sqrt( g1 );
    g2sq = sqrt( g2 );
    
    /* symmetric combinations */
    us = g1sq + g2sq;
    ws = g1sq * g2sq;
    vs = g0sq * us + ws;
    us += g0sq;
    ws *= g0sq;
    
    if( ws < QOP_U3_UNIT_ANALYTIC_EPS ) {
      printf( "QOPQCD WARNING: QOPPC(u3_un_analytic): ws is too small!\n" );
      printf( "  g0 = %28.18f\n", g0 );
      printf( "  g1 = %28.18f\n", g1 );
      printf( "  g2 = %28.18f\n", g2 );
    }
    
    denom = ws * ( us*vs - ws );
    
    /* constants in inverse root expression */
    f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
    f1 = ( 2*us*vs - ws - us*us*us ) / denom;
    f2 = us / denom;
    
    nflops += 30;

    /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
    QLA_M_eq_r_times_M( &S1, &f2, &Q2 );
    QLA_M_eq_r_times_M_plus_M( &S2, &f1, &Q, &S1 );
    QLA_c_peq_r(QLA_elem_M(S2,0,0), f0);
    QLA_c_peq_r(QLA_elem_M(S2,1,1), f0);
    QLA_c_peq_r(QLA_elem_M(S2,2,2), f0);

    nflops += 18 + 36 + 3;
    
    /* W = V*S2 */
    QLA_M_eq_M_times_M( W, V, &S2 );
    
  } // end of SVD related if

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//  printf("Exit  QOPPC(u3_un_analytic) in hisq_reunit_utilities_p.c\n");
#endif

  info->final_flop += nflops;

  return svd_calls;
}



// Analytic derivative of the unitarized matrix with respect to the original:
//   dW/dV and d(W^+)/dV, where W=V(V^+V)^-1/2 */
void QOPPC(u3_un_der_analytic)( QOP_info_t *info, QLA_ColorMatrix *V, 
				QLA_ColorTensor4 *dwdv, QLA_ColorTensor4 *dwdagdv,
				int *svd_calls, int *ff_counter) {
  int i, j, m, n, perform_svd;
  QLA_Complex det, der, ctmp, ctmp2;
  QLA_Real det_check=0, c0, c1, c2, S, g0, g1, g2, R, /*R2, CQ3,*/ S3, RoS, theta, theta3, pi23;
  QLA_Real g0sq, g1sq, g2sq, us, vs, ws, f0, f1, f2, denom;
  QLA_Real u2, u3, u4, u5, u6, u7, u8, v2, v3, v4, v5, v6, w2, w3, w4, w5;
  QLA_Real b00, b01, b02, b11, b12, b22, denom3;
  QLA_ColorMatrix Q12, Q, Q2, Q3, Uleft, Vright, W, S1, S2;
  QLA_ColorMatrix VVd, VQ, QVd, QQVd, VQQ, VQVd, PVd, RVd, SVd, Vd;
  QLA_Real sigma[3];
  size_t nflops = 0;

  // No SVD initially
  perform_svd=0;

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//  printf("Enter QOPPC(u3_un_der_analytic) in hisq_reunit_utilities_p.c\n");
#endif

  if(QOP_hisq_links.reunit_allow_svd){
    // get determinant for future comparison
    det = QOPPC(su3_mat_det)( V );
    // det_check = |det|^2
    det_check = det.real*det.real + det.imag*det.imag;
  }

  /* adjoint */
  QLA_M_eq_Ma( &Vd, V );

  /* Hermitian matrix: Q=V^+ V */
  QLA_M_eq_Ma_times_M( &Q, V, V );
  nflops += 198;

  /* Q^2 */
  QLA_M_eq_M_times_M( &Q2, &Q, &Q );
  nflops += 198;

  /* Q^3 */
  QLA_M_eq_M_times_M( &Q3, &Q2, &Q );
  nflops += 198;

  /* (real) traces */
  QLA_R_eq_re_trace_M(&c0, &Q);
  QLA_R_eq_re_trace_M(&c1, &Q2);
  c1 /= 2;
  QLA_R_eq_re_trace_M(&c2, &Q3);
  c2 /= 3;

  S = c1/3 - c0 * (c0/18);
  nflops += 24;

  if( fabs(S)<QOP_U3_UNIT_ANALYTIC_EPS ) {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
    nflops += 3;
  }
  else {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    //R2 = R*R;
    //CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S*S*S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;

    nflops += 14;

    if( !( fabs(RoS)<1.0 ) ) {
      if( R>0 ) {
        theta = 0.0;
      }
      else {
        theta = QOP_PI;
      }
    } 
    else {
      theta = acos( RoS );
      if(isnan(theta)){
        printf("Hit NaN in QOPPC(u3_un_der_analytic)\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for( i=0; i<3; i++ ) {
          for( j=0; j<3; j++ ) {
            printf( "%24.18g %24.18g\n", QLA_real(QLA_elem_M(*V,i,j)), QLA_imag(QLA_elem_M(*V,i,j)) );
          }
	  nflops += 1;
        }
        QDP_abort(0);
      }
    }

    /* eigenvalues of Q */
    theta3 = theta/3;
    pi23 = (2*QOP_PI) / 3;
    g0 = c0/3 + 2 * S * cos( theta3 );
    g1 = c0/3 + 2 * S * cos( theta3 + pi23 );
    g2 = c0/3 + 2 * S * cos( theta3 + 2*pi23 );

    nflops += 20;
  }

  if(QOP_hisq_links.reunit_allow_svd){

    if(!QOP_hisq_links.reunit_svd_only){

      /* conditions to call SVD */
      if(det_check!=0) {
	if( fabs(det_check-g0*g1*g2)/fabs(det_check) > 
	    QOP_hisq_links.reunit_svd_rel_error ) {
	  perform_svd=1;
	  nflops += 4;
	}
      }
      if( det_check < QOP_hisq_links.reunit_svd_abs_error ) {
	perform_svd=1;
      }
    } else {

      /* exclusively use SVD for finding eigenvalues and reunitarization,
	 this is slow since Q, Q^2 and Q^3 are calculated anyway;
	 this option is for testing: under normal circumstances SVD is
	 rarely used, which makes it harder to test, therefore one can
	 SVD with this switch */
      perform_svd=1;
    }
  } /* reunit_allow_svd */

  if( 0!=perform_svd ) {
    
    /* call SVD */
    QOPPC(svd3x3)( info, V, sigma, &Uleft, &Vright);
    
    (*svd_calls)++;
    
    if(!QOP_hisq_links.reunit_svd_only && QOP_hisq_links.svd_values_info){
      printf("*** QOPQDP reunitarization while calculating fat links ***\n");
      printf("*** Resort to using svd (force) ***\n");
      printf("*** printing from node %d ***\n",QDP_this_node);
      printf("Eigenvalues from cubic equation:\n");
      printf( "  g0 = %28.18f\n", g0 );
      printf( "  g1 = %28.18f\n", g1 );
      printf( "  g2 = %28.18f\n", g2 );
      printf("Eigenvalues from singular value decomposition:\n");
      printf( "  g0 = %28.18f\n", sigma[0]*sigma[0] );
      printf( "  g1 = %28.18f\n", sigma[1]*sigma[1] );
      printf( "  g2 = %28.18f\n", sigma[2]*sigma[2] );
    }
    
    g0=sigma[0]*sigma[0];
    g1=sigma[1]*sigma[1];
    g2=sigma[2]*sigma[2];
    
    nflops += 3;

  }

  if(QOP_hisq_ff.force_filter > 0){
    QLA_Real gmin,g_epsilon;
    gmin=g0;
    if(g1<gmin) gmin=g1;
    if(g2<gmin) gmin=g2;
    g_epsilon = QOP_hisq_ff.force_filter;
    if(gmin<g_epsilon) {
      (*ff_counter)++;
      g0 += g_epsilon;
      g1 += g_epsilon;
      g2 += g_epsilon;
      // +3 flops
      
      // modify also Q and Q2 matrices
      for(i=0;i<3;i++) {
	//Q.e[i][i].real+=g_epsilon;
	QLA_c_peq_r(QLA_elem_M(Q,i,i), g_epsilon);
	nflops += 1;
      }
      QLA_M_eq_M_times_M( &Q2, &Q, &Q );
      nflops += 3+198;
      //QOP_info_hisq_force_filter_counter(info)++;
    }
  }

  /* roots of eigenvalues */
  g0sq = sqrt( g0 );
  g1sq = sqrt( g1 );
  g2sq = sqrt( g2 );
  // +3 flops (sqrt counted as 1 flop)

  /* symmetric combinations */
  us = g1sq + g2sq;
  ws = g1sq * g2sq;
  vs = g0sq * us + ws;
  us += g0sq; 
  ws *= g0sq;
  // +6 flops

  if( ws < QOP_U3_UNIT_ANALYTIC_EPS ) {
    printf( "WARNING: QOPPC(u3_un_der_analytic): ws is too small!\n" );
  }

  denom = ws * ( us*vs - ws );
  // +3 flops

  /* constants in inverse root expression */
  f0 = ( us*vs*vs - ws*(us*us+vs) ) / denom;
  f1 = ( 2*us*vs - ws - us*us*us ) / denom;
  f2 = us / denom;
  // +15 flops

  nflops += 12+15;

  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  QLA_M_eq_r_times_M( &S1, &f2, &Q2 );
  QLA_M_eq_r_times_M_plus_M( &Q12, &f1, &Q, &S1 );
  QLA_c_peq_r(QLA_elem_M(Q12,0,0), f0);
  QLA_c_peq_r(QLA_elem_M(Q12,1,1), f0);
  QLA_c_peq_r(QLA_elem_M(Q12,2,2), f0);

  nflops += 18+36+3;

  /* W = V*S2 */
  QLA_M_eq_M_times_M( &W, V, &Q12 );

  denom3 = 2*denom*denom*denom; // +3 flops
  nflops += 198 + 3;

  /* derivatives of coefficients: B_ij=df_i/dc_j */
  u2 = us*us;
  u3 = u2*us;
  u4 = u3*us;
  u5 = u4*us;
  u6 = u5*us;
  u7 = u6*us;
  u8 = u7*us;
  v2 = vs*vs;
  v3 = v2*vs;
  v4 = v3*vs;
  v5 = v4*vs;
  v6 = v5*vs;
  w2 = ws*ws;
  w3 = w2*ws;
  w4 = w3*ws;
  w5 = w4*ws; // +16 flops
  b00 = -w3*u6 +3*u4*( vs*w3 +v4*ws )
        -u3*( v6 +4*w4 +12*v3*w2 ) + u2*( 16*v2*w3 +3*v5*ws )
        -us*( 8*vs*w4 +3*v4*w2) +w5 +v3*w3;
  b00 /= denom3; // + 33 flops
  b01 = -w2*u7 -v2*ws*u6 +u5*( v4 +6*vs*w2 ) -u4*( 5*w3 +v3*ws )
        -u3*( 2*v5 +6*v2*w2 ) +u2*( 10*vs*w3 +6*v4*ws )
        -us*( 3*w4 +6*v3*w2 ) +2*v2*w3;
  b01 /= denom3; // +38 flops
  b02 = w2*u5 +v2*ws*u4 -u3*( v4 +4*vs*w2 )
        +u2*( 4*w3 +3*v3*ws ) -3*v2*w2*us +vs*w3;
  b02 /= denom3; // +22 flops
  b11 = -ws*u8 -v2*u7 +7*vs*ws*u6 +u5*( 4*v3 -5*w2 ) -16*v2*ws*u4
        -u3*( 4*v4 -16*vs*w2 ) -u2*( 3*w3 -12*v3*ws ) -12*v2*w2*us +3*vs*w3;
  b11 /= denom3; // +37 flops
  b12 = ws*u6 +v2*u5 -5*vs*ws*u4 -u3*( 2*v3 -4*w2 ) +6*v2*ws*u2 -6*vs*w2*us+w3;
  b12 /= denom3; // +22 flops
  b22 = -ws*u4 -v2*u3 +3*vs*ws*u2 -3*w2*us;
  b22 /= denom3; // +12 flops

  nflops += 16+33+38+22+37+22+12;

  /* ** create several building blocks for derivative ** */
  QLA_M_eq_M_times_M( &VQ, V, &Q );
  QLA_M_eq_M_times_Ma( &QVd, &Q, V );
  QLA_M_eq_M_times_Ma( &VVd, V, V );
  QLA_M_eq_M_times_M( &VQQ, V, &Q2 );
  QLA_M_eq_M_times_Ma( &QQVd, &Q2, V );
  QLA_M_eq_M_times_Ma( &VQVd, &VQ, V );
  QLA_M_eq_r_times_M( &S1, &b01, &QVd);
  QLA_M_eq_r_times_M_plus_M( &S2, &b02, &QQVd, &S1 );
  QLA_M_eq_r_times_M_plus_M( &PVd, &b00, &Vd, &S2 );
  QLA_M_eq_r_times_M( &S1, &b11, &QVd );
  QLA_M_eq_r_times_M_plus_M( &S2, &b12, &QQVd, &S1 );
  QLA_M_eq_r_times_M_plus_M( &RVd, &b01, &Vd, &S2 );
  QLA_M_eq_r_times_M( &S1, &b12, &QVd );
  QLA_M_eq_r_times_M_plus_M( &S2, &b22, &QQVd, &S1 );
  QLA_M_eq_r_times_M_plus_M( &SVd, &b02, &Vd, &S2 );

  nflops += 198*6+4*18+5*36;

  /* assemble the derivative rank 4 tensor */
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
	  QLA_c_eq_r(der, 0);
          /* dW/dV */
          if( i==m ) {
	    QLA_c_peq_c(der, QLA_elem_M(Q12, n, j));
	    nflops += 2;
          }
          if( j==n ) {
	    QLA_c_peq_r_times_c(der, f1, QLA_elem_M(VVd,i,m));
	    QLA_c_peq_r_times_c(der, f2, QLA_elem_M(VQVd,i,m));
	    nflops += 4;
          }
          QLA_c_eq_c_times_c( ctmp, QLA_elem_M(*V,i,j), QLA_elem_M(PVd,n,m) );
          QLA_c_peq_c( der, ctmp );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(VQ,i,j),QLA_elem_M(RVd,n,m) );
          QLA_c_peq_c( der, ctmp );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(VQQ,i,j),QLA_elem_M(SVd,n,m) );
          QLA_c_peq_c( der, ctmp );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(VVd,i,m),QLA_elem_M(Q,n,j) );
	  QLA_c_peq_r_times_c(der, f2, ctmp);
          dwdv->t4[i][m][n][j].real = QLA_real(der);
          dwdv->t4[i][m][n][j].imag = QLA_imag(der);
          /* dW^+/dV */
          QLA_c_eq_c_times_c( der,QLA_elem_M(Vd,i,j),QLA_elem_M(PVd,n,m) );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(QVd,i,j),QLA_elem_M(RVd,n,m) );
          QLA_c_peq_c( der, ctmp );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(Vd,i,m),QLA_elem_M(Vd,n,j) );
	  QLA_c_peq_r_times_c(der, f1, ctmp);
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(QQVd,i,j),QLA_elem_M(SVd,n,m) );
          QLA_c_peq_c( der, ctmp );
          QLA_c_eq_c_times_c( ctmp,QLA_elem_M(Vd,i,m),QLA_elem_M(QVd,n,j) );
          QLA_c_eq_c_times_c( ctmp2,QLA_elem_M(Vd,n,j),QLA_elem_M(QVd,i,m) );
	  QLA_c_peq_r_times_c(der, f2, ctmp);
	  QLA_c_peq_r_times_c(der, f2, ctmp2);
          dwdagdv->t4[i][m][n][j].real = QLA_real(der);
          dwdagdv->t4[i][m][n][j].imag = QLA_imag(der);
	  nflops += 10*6+5*2+4*4;
        }
      }
    }
  }

  info->final_flop += nflops;

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//  printf("Exit  QOPPC(u3_un_der_analytic) in hisq_reunit_utilities_p.c\n");
#endif

}





/* **************************************************
   SVD stuff
   ************************************************** */
/* Singular value decomposition (SVD) for
   3x3 complex matrix
   A.Bazavov, Feb 20 2009


   Algorithm sketch:

   SVD is performed in two steps:
     1) 3x3 complex matrix is reduced to real bidiagonal form
        with Householder transformations,
        3x3 matrix requires three left and two right such
        transformations
     2) bidiagonal matrix has the form
        [ b00 b01   0 ]
        [   0 b11 b12 ]
        [   0   0 b22 ]
        It is iteratively diagonalized with QR algorithm with shifts
        (it constructs Given rotations).
        There are many special cases (such as b00==0, b11==0, etc.)
        that are handled separately. If b01==0 then an auxiliary
        routine svd2x2bidiag is used to decompose the lower 2x2 block,
        if b12==0 the same is done for upper block.
        svd2x2bidiag is a separate routine because there are special
        cases for 2x2 superdiagonal matrix that need to be handled.

   This routine needs to be stable for singular matrices. Therefore, most
   of the operations are done to avoid underflow/overflow, for example,
   if norm=sqrt(a^2+b^2) then the calculation proceeds as:
     min=min(a,b), max=max(a,b)
     norm=max*sqrt(1+(min/max)^2)
   and so on.
*/

/* debugging define: prints a lot(!) */
/*#define QOP_SVD3x3_DEBUG*/

/* define precision for chopping small values */
#if ( QOP_PrecisionInt==2 )
#define QOP_SVD3x3_PREC 5e-16
#else
#define QOP_SVD3x3_PREC 5e-7
#endif

/* defines that allow to remap input arrays easily */
#define b00 P[0][0][0]
#define b01 P[0][1][0]
#define b02 P[0][2][0]
#define b10 P[1][0][0]
#define b11 P[1][1][0]
#define b12 P[1][2][0]
#define b20 P[2][0][0]
#define b21 P[2][1][0]
#define b22 P[2][2][0]


/* Input: A -- 3x3 complex matrix,
   Output: sigma[3] -- singular values,
           U,V -- U(3) matrices such, that
           A=U Sigma V^+ */
static int QOPPC(svd3x3)(QOP_info_t *info, QLA_ColorMatrix *A, QLA_Real *sigma, 
			 QLA_ColorMatrix *U, QLA_ColorMatrix *V) {
  QLA_Real Ad[3][3][2], P[3][3][2], Q[3][3][2];
  QLA_Real U1[3][3][2], U2[3][3][2], U3[3][3][2], V1[3][3][2], V2[3][3][2];
  QLA_Real UO3[3][3], VO3[3][3], v[3][2];
  QLA_Real UO2[2][2], VO2[2][2];
  register QLA_Real a, b, c, d, factor, norm, min, max, taure, tauim, beta;
  register QLA_Real m11, m12, m22, dm, lambdamax, cosphi, sinphi, tanphi, cotphi;
  register int i, j, iter;
  size_t nflops = 0;


  /* format of external matrices A, U and V can be arbitrary,
     therefore this routine uses defines (above) to access them
     and never reads A, U and V directly */

  /* original matrix can be in single precision,
     so copy it into double */
#if 0
  Ad[0][0][0]=A00re; Ad[0][0][1]=A00im;
  Ad[0][1][0]=A01re; Ad[0][1][1]=A01im;
  Ad[0][2][0]=A02re; Ad[0][2][1]=A02im;
  Ad[1][0][0]=A10re; Ad[1][0][1]=A10im;
  Ad[1][1][0]=A11re; Ad[1][1][1]=A11im;
  Ad[1][2][0]=A12re; Ad[1][2][1]=A12im;
  Ad[2][0][0]=A20re; Ad[2][0][1]=A20im;
  Ad[2][1][0]=A21re; Ad[2][1][1]=A21im;
  Ad[2][2][0]=A22re; Ad[2][2][1]=A22im;
#endif
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      Ad[i][j][0] = QLA_real(QLA_elem_M(*A,i,j));
      Ad[i][j][1] = QLA_imag(QLA_elem_M(*A,i,j));
    }
  }


  i=0; j=0;

  /* *** Step 1: build first left reflector v,
                 calculate first left rotation U1,
                 apply to the original matrix A *** */
  /* calculate norm of ( A[10] )
                       ( A[20] ) vector
     with minimal loss of accuracy (similar to BLAS) */
  c = 1.;
  factor = fabs( Ad[1][0][0] );
  a = fabs( Ad[1][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + (factor/a)*(factor/a);
      factor = a;
    }
    else {
      c = 1 + (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][0] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  norm = factor*sqrt(c);

  nflops += 15;

  if( norm==0 && Ad[0][0][1]==0 ) { /* no rotation needed */
#ifdef QOP_SVD3x3_DEBUG
    printf("Step 1: no rotation needed\n");
#endif /* QOP_SVD3x3_DEBUG */
    U1[0][0][0]=1.; U1[0][0][1]=0.;
    U1[0][1][0]=0.; U1[0][1][1]=0.;
    U1[0][2][0]=0.; U1[0][2][1]=0.;
    U1[1][0][0]=0.; U1[1][0][1]=0.;
    U1[1][1][0]=1.; U1[1][1][1]=0.;
    U1[1][2][0]=0.; U1[1][2][1]=0.;
    U1[2][0][0]=0.; U1[2][0][1]=0.;
    U1[2][1][0]=0.; U1[2][1][1]=0.;
    U1[2][2][0]=1.; U1[2][2][1]=0.;
    P[0][0][0]=Ad[0][0][0]; P[0][0][1]=Ad[0][0][1];
    P[1][0][0]=Ad[1][0][0]; P[1][0][1]=Ad[1][0][1];
    P[2][0][0]=Ad[2][0][0]; P[2][0][1]=Ad[2][0][1];
    P[0][1][0]=Ad[0][1][0]; P[0][1][1]=Ad[0][1][1];
    P[1][1][0]=Ad[1][1][0]; P[1][1][1]=Ad[1][1][1];
    P[2][1][0]=Ad[2][1][0]; P[2][1][1]=Ad[2][1][1];
    P[0][2][0]=Ad[0][2][0]; P[0][2][1]=Ad[0][2][1];
    P[1][2][0]=Ad[1][2][0]; P[1][2][1]=Ad[1][2][1];
    P[2][2][0]=Ad[2][2][0]; P[2][2][1]=Ad[2][2][1];
  }
  else {


    /* get the norm of full first column of A matrix */
    c=1.;
    factor = norm;
    a = fabs( Ad[0][0][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Ad[0][0][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of first column */
    if( Ad[0][0][0]>0 ) {
      beta = -beta;
    }

#ifdef QOP_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QOP_SVD3x3_DEBUG */


    /* a=Re(A_00-beta), b=Im(A_00-beta) */
    a=Ad[0][0][0]-beta; b=Ad[0][0][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation)
       v_0=1 */
    v[0][0]=1.; v[0][1]=0.;
    /* v_1=A_10/(A_00-beta)=A_10/(a+ib)=(A_10*(a-ib))/norm^2=(A_10/norm)*((a-ib)/norm)
          =(A_10/norm)*(c-id)=|a=Re(A_10)/norm,b=Im(A_10)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[1][0][0]/norm; b=Ad[1][0][1]/norm;
    v[1][0]=a*c+b*d;
    v[1][1]=b*c-a*d;
    /* v_2=A_20/(A_00-beta)=A_20/(a+ib)=(A_20*(a-ib))/norm^2=(A_20/norm)*((a-ib)/norm)
          =(A_20/norm)*(c-id)=|a=Re(A_20)/norm,b=Im(A_20)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[2][0][0]/norm; b=Ad[2][0][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;
#ifdef QOP_SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* QOP_SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Ad[0][0][0])/beta;
    tauim=(Ad[0][0][1])/beta;


    /* assemble left unitary matrix U1=I-tau^+*v*v^+ (store in U1[3][3][2])
       U1_00=A_00/beta */
    U1[0][0][0]=(Ad[0][0][0])/beta;
    U1[0][0][1]=(Ad[0][0][1])/beta;
    /* U1_10=A_10/beta */
    U1[1][0][0]=(Ad[1][0][0])/beta;
    U1[1][0][1]=(Ad[1][0][1])/beta;
    /* U1_20=A_20/beta */
    U1[2][0][0]=(Ad[2][0][0])/beta;
    U1[2][0][1]=(Ad[2][0][1])/beta;
    /* U1_01=-tau^+*v_1^+=-(tau*v_1)^+ */
    U1[0][1][0]=-(taure*v[1][0]-tauim*v[1][1]);
    U1[0][1][1]=taure*v[1][1]+tauim*v[1][0];
    /* U1_11=1-tau^+*v_1*v_1^+ */
    a=v[1][0]*v[1][0]+v[1][1]*v[1][1];
    U1[1][1][0]=1-taure*a;
    U1[1][1][1]=tauim*a;
    /* U1_21=-tau^+*v_2*v_1^+ */
      /* v_2*v_1^+ */
      a=v[2][0]*v[1][0]+v[2][1]*v[1][1];
      b=-v[2][0]*v[1][1]+v[2][1]*v[1][0];
    U1[2][1][0]=-(taure*a+tauim*b);
    U1[2][1][1]=-(taure*b-tauim*a);
    /* U1_02=-tau^+*v_2^+=-(tau*v_2)^+ */
    U1[0][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U1[0][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U1_12=-tau^+*v_1*v_2^+ */
      /* v_1*v_2^+ */
      a=v[1][0]*v[2][0]+v[1][1]*v[2][1];
      b=-v[1][0]*v[2][1]+v[1][1]*v[2][0];
    U1[1][2][0]=-(taure*a+tauim*b);
    U1[1][2][1]=-(taure*b-tauim*a);
    /* U1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U1[2][2][0]=1-taure*a;
    U1[2][2][1]=tauim*a;

    nflops += 91;


    /* apply the transformation to A matrix and store the result in P
       P=U^+A */
    P[0][0][0]=beta;
    P[0][0][1]=0;
    P[1][0][0]=0;
    P[1][0][1]=0;
    P[2][0][0]=0;
    P[2][0][1]=0;
    /* P_01=U1_00^+*A_01+U1_10^+*A_11+U1_20^+*A_21 */
    P[0][1][0]=U1[0][0][0]*Ad[0][1][0]+U1[0][0][1]*Ad[0][1][1]
              +U1[1][0][0]*Ad[1][1][0]+U1[1][0][1]*Ad[1][1][1]
              +U1[2][0][0]*Ad[2][1][0]+U1[2][0][1]*Ad[2][1][1];
    P[0][1][1]=U1[0][0][0]*Ad[0][1][1]-U1[0][0][1]*Ad[0][1][0]
              +U1[1][0][0]*Ad[1][1][1]-U1[1][0][1]*Ad[1][1][0]
              +U1[2][0][0]*Ad[2][1][1]-U1[2][0][1]*Ad[2][1][0];
    /* P_02=U1_00^+*A_02+U1_10^+*A_12+U1_20^+*A_22 */
    P[0][2][0]=U1[0][0][0]*Ad[0][2][0]+U1[0][0][1]*Ad[0][2][1]
              +U1[1][0][0]*Ad[1][2][0]+U1[1][0][1]*Ad[1][2][1]
              +U1[2][0][0]*Ad[2][2][0]+U1[2][0][1]*Ad[2][2][1];
    P[0][2][1]=U1[0][0][0]*Ad[0][2][1]-U1[0][0][1]*Ad[0][2][0]
              +U1[1][0][0]*Ad[1][2][1]-U1[1][0][1]*Ad[1][2][0]
              +U1[2][0][0]*Ad[2][2][1]-U1[2][0][1]*Ad[2][2][0];
    /* P_11=U1_01^+*A_01+U1_11^+*A_11+U1_21^+*A_21 */
    P[1][1][0]=U1[0][1][0]*Ad[0][1][0]+U1[0][1][1]*Ad[0][1][1]
              +U1[1][1][0]*Ad[1][1][0]+U1[1][1][1]*Ad[1][1][1]
              +U1[2][1][0]*Ad[2][1][0]+U1[2][1][1]*Ad[2][1][1];
    P[1][1][1]=U1[0][1][0]*Ad[0][1][1]-U1[0][1][1]*Ad[0][1][0]
              +U1[1][1][0]*Ad[1][1][1]-U1[1][1][1]*Ad[1][1][0]
              +U1[2][1][0]*Ad[2][1][1]-U1[2][1][1]*Ad[2][1][0];
    /* P_12=U1_01^+*A_02+U1_11^+*A_12+U1_21^+*A_22 */
    P[1][2][0]=U1[0][1][0]*Ad[0][2][0]+U1[0][1][1]*Ad[0][2][1]
              +U1[1][1][0]*Ad[1][2][0]+U1[1][1][1]*Ad[1][2][1]
              +U1[2][1][0]*Ad[2][2][0]+U1[2][1][1]*Ad[2][2][1];
    P[1][2][1]=U1[0][1][0]*Ad[0][2][1]-U1[0][1][1]*Ad[0][2][0]
              +U1[1][1][0]*Ad[1][2][1]-U1[1][1][1]*Ad[1][2][0]
              +U1[2][1][0]*Ad[2][2][1]-U1[2][1][1]*Ad[2][2][0];
    /* P_21=U1_02^+*A_01+U1_12^+*A_11+U1_22^+*A_21 */
    P[2][1][0]=U1[0][2][0]*Ad[0][1][0]+U1[0][2][1]*Ad[0][1][1]
              +U1[1][2][0]*Ad[1][1][0]+U1[1][2][1]*Ad[1][1][1]
              +U1[2][2][0]*Ad[2][1][0]+U1[2][2][1]*Ad[2][1][1];
    P[2][1][1]=U1[0][2][0]*Ad[0][1][1]-U1[0][2][1]*Ad[0][1][0]
              +U1[1][2][0]*Ad[1][1][1]-U1[1][2][1]*Ad[1][1][0]
              +U1[2][2][0]*Ad[2][1][1]-U1[2][2][1]*Ad[2][1][0];
    /* P_22=U1_02^+*A_02+U1_12^+*A_12+U1_22^+*A_22 */
    P[2][2][0]=U1[0][2][0]*Ad[0][2][0]+U1[0][2][1]*Ad[0][2][1]
              +U1[1][2][0]*Ad[1][2][0]+U1[1][2][1]*Ad[1][2][1]
              +U1[2][2][0]*Ad[2][2][0]+U1[2][2][1]*Ad[2][2][1];
    P[2][2][1]=U1[0][2][0]*Ad[0][2][1]-U1[0][2][1]*Ad[0][2][0]
              +U1[1][2][0]*Ad[1][2][1]-U1[1][2][1]*Ad[1][2][0]
              +U1[2][2][0]*Ad[2][2][1]-U1[2][2][1]*Ad[2][2][0];

    nflops += 9*12;

  }
#ifdef QOP_SVD3x3_DEBUG
printf("Left unitary matrix U1:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U1[%d][%d].re=%26.18e  U1[%d][%d].im=%26.18e\n",
              i, j, U1[i][j][0], i, j, U1[i][j][1] );
    }
#endif /* QOP_SVD3x3_DEBUG */



  /* *** Step 2: build first right reflector v,
                 calculate first right rotation V1,
                 apply to the matrix P from step 1 *** */
  /* calculate norm of ( P[02] )
     with minimal loss of accuracy */
  a=fabs( P[0][2][0] ); b=fabs( P[0][2][1] );
  /* norm=sqrt(a^2+b^2) */
  if( a>b ) {
    max=a; min=b;
  }
  else {
    max=b; min=a;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && P[0][1][1]==0 ) { /* no rotation needed */
#ifdef QOP_SVD3x3_DEBUG
    printf("Step 2: no rotation needed\n");
#endif /* QOP_SVD3x3_DEBUG */
    V1[0][0][0]=1.; V1[0][0][1]=0.;
    V1[0][1][0]=0.; V1[0][1][1]=0.;
    V1[0][2][0]=0.; V1[0][2][1]=0.;
    V1[1][0][0]=0.; V1[1][0][1]=0.;
    V1[1][1][0]=1.; V1[1][1][1]=0.;
    V1[1][2][0]=0.; V1[1][2][1]=0.;
    V1[2][0][0]=0.; V1[2][0][1]=0.;
    V1[2][1][0]=0.; V1[2][1][1]=0.;
    V1[2][2][0]=1.; V1[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* get the norm of (P_01 P_02) row vector */
    c=1.;
    factor = norm;
    a = fabs( P[0][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( P[0][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (P_01 P_02) row vector */
    if( P[0][1][0]>0 ) {
      beta = -beta;
    }

#ifdef QOP_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QOP_SVD3x3_DEBUG */


    /* a=Re(P_01^+-beta), b=Im(P_01^+-beta) */
    a=P[0][1][0]-beta; b=-P[0][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=P_02^+/(P_01^+-beta)=P_02^+/(a+ib)=(P_02^+*(a-ib))/norm^2=(P_02^+/norm)*((a-ib)/norm)
          =(P_02^+/norm)*(c-id)=|a=Re(P_02^+)/norm,b=Im(P_02^+)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=P[0][2][0]/norm; b=-P[0][2][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;

#ifdef QOP_SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* QOP_SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-P[0][1][0])/beta;
    tauim=-P[0][1][1]/beta;

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V1[0][0][0]=1.;
    V1[0][0][1]=0.;
    V1[1][0][0]=0.;
    V1[1][0][1]=0.;
    V1[2][0][0]=0.;
    V1[2][0][1]=0.;
    V1[0][1][0]=0.;
    V1[0][1][1]=0.;
    V1[0][2][0]=0.;
    V1[0][2][1]=0.;
    /* V1_11=P_01^+/beta */
    V1[1][1][0]=P[0][1][0]/beta;
    V1[1][1][1]=-P[0][1][1]/beta;
    /* V1_21=P_02^+/beta */
    V1[2][1][0]=P[0][2][0]/beta;
    V1[2][1][1]=-P[0][2][1]/beta;
    /* V1_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    V1[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    V1[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* V1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    V1[2][2][0]=1-taure*a;
    V1[2][2][1]=tauim*a;


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=beta;
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    /* Q_11=P_11*V1_11+P_12*V_21 */
    Q[1][1][0]=P[1][1][0]*V1[1][1][0]-P[1][1][1]*V1[1][1][1]
              +P[1][2][0]*V1[2][1][0]-P[1][2][1]*V1[2][1][1];
    Q[1][1][1]=P[1][1][0]*V1[1][1][1]+P[1][1][1]*V1[1][1][0]
              +P[1][2][0]*V1[2][1][1]+P[1][2][1]*V1[2][1][0];
    /* Q_12=P_11*V1_12+P_12*V_22 */
    Q[1][2][0]=P[1][1][0]*V1[1][2][0]-P[1][1][1]*V1[1][2][1]
              +P[1][2][0]*V1[2][2][0]-P[1][2][1]*V1[2][2][1];
    Q[1][2][1]=P[1][1][0]*V1[1][2][1]+P[1][1][1]*V1[1][2][0]
              +P[1][2][0]*V1[2][2][1]+P[1][2][1]*V1[2][2][0];
    /* Q_21=P_21*V1_11+P_22*V_21 */
    Q[2][1][0]=P[2][1][0]*V1[1][1][0]-P[2][1][1]*V1[1][1][1]
              +P[2][2][0]*V1[2][1][0]-P[2][2][1]*V1[2][1][1];
    Q[2][1][1]=P[2][1][0]*V1[1][1][1]+P[2][1][1]*V1[1][1][0]
              +P[2][2][0]*V1[2][1][1]+P[2][2][1]*V1[2][1][0];
    /* Q_22=P_21*V1_12+P_22*V_22 */
    Q[2][2][0]=P[2][1][0]*V1[1][2][0]-P[2][1][1]*V1[1][2][1]
              +P[2][2][0]*V1[2][2][0]-P[2][2][1]*V1[2][2][1];
    Q[2][2][1]=P[2][1][0]*V1[1][2][1]+P[2][1][1]*V1[1][2][0]
              +P[2][2][0]*V1[2][2][1]+P[2][2][1]*V1[2][2][0];

    nflops += 15 + 7*8;
  }
#ifdef QOP_SVD3x3_DEBUG
printf("Right unitary matrix V1:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "V1[%d][%d].re=%26.18e  V1[%d][%d].im=%26.18e\n",
              i, j, V1[i][j][0], i, j, V1[i][j][1] );
    }
#endif /* QOP_SVD3x3_DEBUG */



  /* *** Step 3: build second left reflector v,
                 calculate second left rotation U2,
                 apply to the matrix Q *** */
  /* calculate norm of ( Q[21] )
     with minimal loss of accuracy (similar to BLAS) */
  c=fabs(Q[2][1][0]); d=fabs(Q[2][1][1]);
  if( c>d ) {
    max=c; min=d;
  }
  else {
    max=d; min=c;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && Q[1][1][1]==0 ) { /* no rotation needed */
#ifdef QOP_SVD3x3_DEBUG
    printf("Step 3: no rotation needed\n");
#endif /* QOP_SVD3x3_DEBUG */
    U2[0][0][0]=1.; U2[0][0][1]=0.;
    U2[0][1][0]=0.; U2[0][1][1]=0.;
    U2[0][2][0]=0.; U2[0][2][1]=0.;
    U2[1][0][0]=0.; U2[1][0][1]=0.;
    U2[1][1][0]=1.; U2[1][1][1]=0.;
    U2[1][2][0]=0.; U2[1][2][1]=0.;
    U2[2][0][0]=0.; U2[2][0][1]=0.;
    U2[2][1][0]=0.; U2[2][1][1]=0.;
    U2[2][2][0]=1.; U2[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* get the norm of (Q_11 Q_21) column vector */
    c=1.;
    factor = norm;
    a = fabs( Q[1][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Q[1][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (Q_11 Q_21) column vector */
    if( Q[1][1][0]>0 ) {
      beta = -beta;
    }

#ifdef QOP_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QOP_SVD3x3_DEBUG */


    /* a=Re(Q_11-beta), b=Im(Q_11-beta) */
    a=Q[1][1][0]-beta; b=Q[1][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;

    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=Q_21/(Q_11-beta)=Q_21/(a+ib)=(Q_21*(a-ib))/norm^2=(Q_21/norm)*((a-ib)/norm)
          =(Q_21/norm)*(c-id)=|a=Re(Q_21)/norm,b=Im(Q_21)/norm|=(a+ib)*(c-id)
          =(a*c+b*d)+i(b*c-a*d) */
    a=Q[2][1][0]/norm; b=Q[2][1][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;
#ifdef QOP_SVD3x3_DEBUG
for(i=0;i<3;i++) {
  printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
}
#endif /* QOP_SVD3x3_DEBUG */


    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Q[1][1][0])/beta;
    tauim=Q[1][1][1]/beta;


    /* assemble right unitary matrix U2=I-tau^+*v*v^+ (store in U2[3][3][2]) */
    U2[0][0][0]=1.;
    U2[0][0][1]=0.;
    U2[1][0][0]=0.;
    U2[1][0][1]=0.;
    U2[2][0][0]=0.;
    U2[2][0][1]=0.;
    U2[0][1][0]=0.;
    U2[0][1][1]=0.;
    U2[0][2][0]=0.;
    U2[0][2][1]=0.;
    /* U2_11=Q_11/beta */
    U2[1][1][0]=Q[1][1][0]/beta;
    U2[1][1][1]=Q[1][1][1]/beta;
    /* U2_21=Q_21/beta */
    U2[2][1][0]=Q[2][1][0]/beta;
    U2[2][1][1]=Q[2][1][1]/beta;
    /* U2_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    U2[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U2[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U2_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U2[2][2][0]=1-taure*a;
    U2[2][2][1]=tauim*a;
#ifdef QOP_SVD3x3_DEBUG
printf("Left unitary matrix U2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U2[%d][%d].re=%26.18e  U2[%d][%d].im=%26.18e\n",
              i, j, U2[i][j][0], i, j, U2[i][j][1] );
    }
#endif /* QOP_SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=beta;
    P[1][1][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    /* P_12=U2_11^+*Q_12+U2_21^+*Q_22 */
    P[1][2][0]=U2[1][1][0]*Q[1][2][0]+U2[1][1][1]*Q[1][2][1]
              +U2[2][1][0]*Q[2][2][0]+U2[2][1][1]*Q[2][2][1];
    P[1][2][1]=U2[1][1][0]*Q[1][2][1]-U2[1][1][1]*Q[1][2][0]
              +U2[2][1][0]*Q[2][2][1]-U2[2][1][1]*Q[2][2][0];
    /* P_22=U2_12^+*Q_12+U2_22^+*Q_22 */
    P[2][2][0]=U2[1][2][0]*Q[1][2][0]+U2[1][2][1]*Q[1][2][1]
              +U2[2][2][0]*Q[2][2][0]+U2[2][2][1]*Q[2][2][1];
    P[2][2][1]=U2[1][2][0]*Q[1][2][1]-U2[1][2][1]*Q[1][2][0]
              +U2[2][2][0]*Q[2][2][1]-U2[2][2][1]*Q[2][2][0];

    nflops += 15 + 7*8;

  }



  /* *** Step 4: build second right reflector v,
                 calculate second right rotation V2,
                 apply to the matrix P *** */
  if( P[1][2][1]==0 ) { /* no rotation needed */
#ifdef QOP_SVD3x3_DEBUG
    printf("Step 4: no rotation needed\n");
#endif /* QOP_SVD3x3_DEBUG */
    V2[0][0][0]=1.; V2[0][0][1]=0.;
    V2[0][1][0]=0.; V2[0][1][1]=0.;
    V2[0][2][0]=0.; V2[0][2][1]=0.;
    V2[1][0][0]=0.; V2[1][0][1]=0.;
    V2[1][1][0]=1.; V2[1][1][1]=0.;
    V2[1][2][0]=0.; V2[1][2][1]=0.;
    V2[2][0][0]=0.; V2[2][0][1]=0.;
    V2[2][1][0]=0.; V2[2][1][1]=0.;
    V2[2][2][0]=1.; V2[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* calculate norm of ( P[12] ) */
    c=fabs(P[1][2][0]); d=fabs(P[1][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( P[1][2][0]>0 ) {
      beta = -beta;
    }

#ifdef QOP_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QOP_SVD3x3_DEBUG */

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V2[0][0][0]=1.;
    V2[0][0][1]=0.;
    V2[1][0][0]=0.;
    V2[1][0][1]=0.;
    V2[2][0][0]=0.;
    V2[2][0][1]=0.;
    V2[0][1][0]=0.;
    V2[0][1][1]=0.;
    V2[0][2][0]=0.;
    V2[0][2][1]=0.;
    V2[1][1][0]=1.;
    V2[1][1][1]=0.;
    V2[2][1][0]=0.;
    V2[2][1][1]=0.;
    V2[1][2][0]=0.;
    V2[1][2][1]=0.;
    /* V2_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    V2[2][2][0]=P[1][2][0]/beta;
    V2[2][2][1]=-P[1][2][1]/beta;
#ifdef QOP_SVD3x3_DEBUG
printf("Right unitary matrix V2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "V2[%d][%d].re=%26.18e  V2[%d][%d].im=%26.18e\n",
              i, j, V2[i][j][0], i, j, V2[i][j][1] );
    }
#endif /* QOP_SVD3x3_DEBUG */


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=P[0][1][0];
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    Q[1][1][0]=P[1][1][0];
    Q[1][1][1]=0.;
    Q[1][2][0]=beta;
    Q[1][2][1]=0.;
    Q[2][1][0]=0.;
    Q[2][1][1]=0.;
    /* Q_22=P_22*V2_22 */
    Q[2][2][0]=P[2][2][0]*V2[2][2][0]-P[2][2][1]*V2[2][2][1];
    Q[2][2][1]=P[2][2][0]*V2[2][2][1]+P[2][2][1]*V2[2][2][0];

    nflops += 12;
  }



  /* *** Step 5: build third left reflector v,
                 calculate third left rotation U3,
                 apply to the matrix P *** */
  if( Q[2][2][1]==0 ) { /* no rotation needed */
#ifdef QOP_SVD3x3_DEBUG
    printf("Step 5: no rotation needed\n");
#endif /* QOP_SVD3x3_DEBUG */
    U3[0][0][0]=1.; U3[0][0][1]=0.;
    U3[0][1][0]=0.; U3[0][1][1]=0.;
    U3[0][2][0]=0.; U3[0][2][1]=0.;
    U3[1][0][0]=0.; U3[1][0][1]=0.;
    U3[1][1][0]=1.; U3[1][1][1]=0.;
    U3[1][2][0]=0.; U3[1][2][1]=0.;
    U3[2][0][0]=0.; U3[2][0][1]=0.;
    U3[2][1][0]=0.; U3[2][1][1]=0.;
    U3[2][2][0]=1.; U3[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* calculate norm of ( Q[22] ) */
    c=fabs(Q[2][2][0]); d=fabs(Q[2][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( Q[2][2][0]>0 ) {
      beta = -beta;
    }

#ifdef QOP_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QOP_SVD3x3_DEBUG */

    /* assemble left unitary matrix U3=I-tau^+*v*v^+ (store in U3[3][3][2]) */
    U3[0][0][0]=1.;
    U3[0][0][1]=0.;
    U3[1][0][0]=0.;
    U3[1][0][1]=0.;
    U3[2][0][0]=0.;
    U3[2][0][1]=0.;
    U3[0][1][0]=0.;
    U3[0][1][1]=0.;
    U3[0][2][0]=0.;
    U3[0][2][1]=0.;
    U3[1][1][0]=1.;
    U3[1][1][1]=0.;
    U3[2][1][0]=0.;
    U3[2][1][1]=0.;
    U3[1][2][0]=0.;
    U3[1][2][1]=0.;
    /* U3_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    U3[2][2][0]=Q[2][2][0]/beta;
    U3[2][2][1]=Q[2][2][1]/beta;
#ifdef QOP_SVD3x3_DEBUG
printf("Left unitary matrix U3:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U3[%d][%d].re=%26.18e  U3[%d][%d].im=%26.18e\n",
              i, j, U3[i][j][0], i, j, U3[i][j][1] );
    }
#endif /* QOP_SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=Q[1][1][0];
    P[1][1][1]=0.;
    P[1][2][0]=Q[1][2][0];
    P[1][2][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    P[2][2][0]=beta;
    P[2][2][1]=0.;

    nflops += 6;

  }




  /* *** This part starts with a bidiagonal matrix and uses
         QR algorithm with shifts to eliminate the superdiagonal *** */
  /* prepare left and right real orthogonal matrices that
     accumulate Givens rotations from QR algorithm */
  UO3[0][0]=1.; UO3[0][1]=0.; UO3[0][2]=0.;
  UO3[1][0]=0.; UO3[1][1]=1.; UO3[1][2]=0.;
  UO3[2][0]=0.; UO3[2][1]=0.; UO3[2][2]=1.;
  VO3[0][0]=1.; VO3[0][1]=0.; VO3[0][2]=0.;
  VO3[1][0]=0.; VO3[1][1]=1.; VO3[1][2]=0.;
  VO3[2][0]=0.; VO3[2][1]=0.; VO3[2][2]=1.;

  iter=0;

#ifdef QOP_SVD3x3_DEBUG
printf( "QR iteration: %d\n", iter );
printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* QOP_SVD3x3_DEBUG */

  do {

    iter++;
    if(iter>300) return 1;

    /* chop small superdiagonal elements */
    if( fabs(b01) < QOP_SVD3x3_PREC*(fabs(b00)+fabs(b11)) ) {
      b01=0;
    }
    if( fabs(b12) < QOP_SVD3x3_PREC*(fabs(b00)+fabs(b22)) ) {
      b12=0;
    }

    nflops += 4;

    /* Cases:
       b01=b12=0 -- matrix is already diagonalized,
       b01=0 -- need to work with 2x2 lower block,
       b12=0 -- need to work with 2x2 upper block,
       else -- normal iteration */
    if( !(b01==0 && b12==0) ) {
      if( b01==0 ) {
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case b01==0\n" );
#endif /* QOP_SVD3x3_DEBUG */
        /* need to diagonalize 2x2 lower block */
       QOPPC(svd2x2bidiag)(info, &b11, &b12, &b22, UO2, VO2 );

        /* multiply left UO3 matrix */
        for(i=0;i<3;i++) {
          a=UO3[i][1]; b=UO3[i][2];
          UO3[i][1]=a*UO2[0][0]+b*UO2[1][0];
          UO3[i][2]=a*UO2[0][1]+b*UO2[1][1];
        }
        /* multiply right VO3 matrix */
        for(i=0;i<3;i++) {
          a=VO3[i][1]; b=VO3[i][2];
          VO3[i][1]=a*VO2[0][0]+b*VO2[1][0];
          VO3[i][2]=a*VO2[0][1]+b*VO2[1][1];
        }

	nflops += 36;

      }
      else {
        if( b12==0 ) {
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case b12==0\n" );
#endif /* QOP_SVD3x3_DEBUG */
          /* need to diagonalize 2x2 upper block */
         QOPPC(svd2x2bidiag)(info, &b00, &b01, &b11, UO2, VO2 );

          /* multiply left UO3 matrix */
          for(i=0;i<3;i++) {
            a=UO3[i][0]; b=UO3[i][1];
            UO3[i][0]=a*UO2[0][0]+b*UO2[1][0];
            UO3[i][1]=a*UO2[0][1]+b*UO2[1][1];
          }
          /* multiply right VO3 matrix */
          for(i=0;i<3;i++) {
            a=VO3[i][0]; b=VO3[i][1];
            VO3[i][0]=a*VO2[0][0]+b*VO2[1][0];
            VO3[i][1]=a*VO2[0][1]+b*VO2[1][1];
          }

	  nflops += 36;
        }
        else {
          /* normal 3x3 iteration */

          /* QR shift does not work if there are zeros
             on the diagonal, therefore first check
             for special cases: b00==0 or b11==0 or b22==0 */

          if( b00==0 ) {
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case b00==0\n" );
#endif /* QOP_SVD3x3_DEBUG */
            /* b01 can be rotated away to create b02,
               and then b02 can be rotated away
               (both are left rotations) */
            if( fabs(b01)>fabs(b11) ) {
              cotphi=b11/b01;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b01/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b11=b01*sinphi+b11*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b01=0.;
            if( fabs(b02)>fabs(b22) ) {
              cotphi=b22/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][2];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b02*sinphi+b22*cosphi;
            b02=0.;

	    nflops += 56;
          }
          else if( b11==0 ) {
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case b11==0\n" );
#endif /* QOP_SVD3x3_DEBUG */
            /* b12 is rotated away with left rotation */
            if( fabs(b12)>fabs(b22) ) {
              cotphi=b22/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b12*sinphi+b22*cosphi;
            b12=0.;

	    nflops += 27;
          }
          else if( b22==0 ) {
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case b22==0\n" );
#endif /* QOP_SVD3x3_DEBUG */
            /* b12 is rotated away and b02 appears,
               then b02 is rotated away, both are
               right rotations */
            if( fabs(b12)>fabs(b11) ) {
              cotphi=b11/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b02=-b01*sinphi;
            b01=b01*cosphi;
            b11=b11*cosphi+b12*sinphi;
            b12=0.;
            /* second rotation removes b02 */
            if( fabs(b02)>fabs(b00) ) {
              cotphi=b00/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b00;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][2];
              VO3[i][0]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b00=b00*cosphi+b02*sinphi;
            b02=0.;
	    
	    nflops += 64;
          }
          else {
            /* full iteration with QR shift */
#ifdef QOP_SVD3x3_DEBUG
printf( "Entering case of normal QR iteration\n" );
#endif /* QOP_SVD3x3_DEBUG */

            /* find max eigenvalue of bottom 2x2 minor */
            m11=b11*b11+b01*b01;
            m22=b22*b22+b12*b12;
            m12=b11*b12;
            dm=(m11-m22)/2;

            /* safely calculate sqrt */
            c=fabs(dm); d=fabs(m12);
            if( c>d ) {
              max=c; min=d;
            }
            else {
              max=d; min=c;
            }
            if( min==0 ) {
              norm = max;
            }
            else {
              c = min/max;
              norm = max*sqrt(1+c*c);
            }

            if( dm>=0 ) {
              lambdamax=m22-(m12*m12)/(dm+norm);
            }
            else {
              lambdamax=m22+(m12*m12)/(norm-dm);
            }

            /* calculate first Givens rotation (on the right) */
            a=b00*b00-lambdamax;
            b=b00*b01;
            if( 0==b ) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b)>fabs(a) ) {
                cotphi=-a/b;
                sinphi=1./sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b/a;
                cosphi=1./sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }
	      nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][1];
              VO3[i][0]=a*cosphi-b*sinphi;
              VO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generate b10 */
            a=b00; b=b01;
            b00=a*cosphi-b*sinphi;
            b01=a*sinphi+b*cosphi;
            b10=-b11*sinphi;
            b11=b11*cosphi; 

            /* calculate second Givens rotation (on the left) */
            if(0==b10) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b10)>fabs(b00) ) {
                cotphi=-b00/b10;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b10/b00;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]= a*cosphi-b*sinphi;
              UO3[i][1]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b00=b00*cosphi-b10*sinphi;
            a=b01; b=b11;
            b01=a*cosphi-b*sinphi;
            b11=a*sinphi+b*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b10=0.;

            /* calculate third Givens rotation (on the right) */
            if(0==b02) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b02)>fabs(b01) ) {
                cotphi=-b01/b02;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b02/b01;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]=a*cosphi-b*sinphi;
              VO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b21 */
            b01=b01*cosphi-b02*sinphi;
            a=b11; b=b12;
            b11=a*cosphi-b*sinphi;
            b12=a*sinphi+b*cosphi;
            b21=-b22*sinphi;
            b22=b22*cosphi;
            b02=0.;

            /* calculate fourth Givens rotation (on the left) */
            if(0==b21) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b21)>fabs(b11) ) {
                cotphi=-b11/b21;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b21/b11;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

	      nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]= a*cosphi-b*sinphi;
              UO3[i][2]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this eliminates b21 */
            b11=b11*cosphi-b21*sinphi;
            a=b12; b=b22;
            b12=a*cosphi-b*sinphi;
            b22=a*sinphi+b*cosphi;
            b21=0.;

	    nflops += 127;
          }
        } /* end of normal 3x3 iteration */
      }
    }
#ifdef QOP_SVD3x3_DEBUG
printf( "QR iteration: %d\n", iter );
printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* QOP_SVD3x3_DEBUG */
  }
  while( b01!=0 || b12!=0 );


  /* make singular values positive */
  if(b00<0) {
    b00=-b00;
    VO3[0][0]=-VO3[0][0];
    VO3[1][0]=-VO3[1][0];
    VO3[2][0]=-VO3[2][0];
  }
  if(b11<0) {
    b11=-b11;
    VO3[0][1]=-VO3[0][1];
    VO3[1][1]=-VO3[1][1];
    VO3[2][1]=-VO3[2][1];
  }
  if(b22<0) {
    b22=-b22;
    VO3[0][2]=-VO3[0][2];
    VO3[1][2]=-VO3[1][2];
    VO3[2][2]=-VO3[2][2];
  }



  /* Q=U1*U2 (U2 is block diagonal with U2_00=1) */
  Q[0][0][0]=U1[0][0][0]; Q[0][0][1]=U1[0][0][1];
  Q[1][0][0]=U1[1][0][0]; Q[1][0][1]=U1[1][0][1];
  Q[2][0][0]=U1[2][0][0]; Q[2][0][1]=U1[2][0][1];
  /* Q_01=U1_01*U2_11+U1_02*U2_21 */
  Q[0][1][0]=U1[0][1][0]*U2[1][1][0]-U1[0][1][1]*U2[1][1][1]
            +U1[0][2][0]*U2[2][1][0]-U1[0][2][1]*U2[2][1][1];
  Q[0][1][1]=U1[0][1][0]*U2[1][1][1]+U1[0][1][1]*U2[1][1][0]
            +U1[0][2][0]*U2[2][1][1]+U1[0][2][1]*U2[2][1][0];
  /* Q_02=U1_01*U2_12+U1_02*U2_22 */
  Q[0][2][0]=U1[0][1][0]*U2[1][2][0]-U1[0][1][1]*U2[1][2][1]
            +U1[0][2][0]*U2[2][2][0]-U1[0][2][1]*U2[2][2][1];
  Q[0][2][1]=U1[0][1][0]*U2[1][2][1]+U1[0][1][1]*U2[1][2][0]
            +U1[0][2][0]*U2[2][2][1]+U1[0][2][1]*U2[2][2][0];
  /* Q_11=U1_11*U2_11+U1_12*U2_21 */
  Q[1][1][0]=U1[1][1][0]*U2[1][1][0]-U1[1][1][1]*U2[1][1][1]
            +U1[1][2][0]*U2[2][1][0]-U1[1][2][1]*U2[2][1][1];
  Q[1][1][1]=U1[1][1][0]*U2[1][1][1]+U1[1][1][1]*U2[1][1][0]
            +U1[1][2][0]*U2[2][1][1]+U1[1][2][1]*U2[2][1][0];
  /* Q_12=U1_11*U2_12+U1_12*U2_22 */
  Q[1][2][0]=U1[1][1][0]*U2[1][2][0]-U1[1][1][1]*U2[1][2][1]
            +U1[1][2][0]*U2[2][2][0]-U1[1][2][1]*U2[2][2][1];
  Q[1][2][1]=U1[1][1][0]*U2[1][2][1]+U1[1][1][1]*U2[1][2][0]
            +U1[1][2][0]*U2[2][2][1]+U1[1][2][1]*U2[2][2][0];
  /* Q_21=U1_21*U2_11+U1_22*U2_21 */
  Q[2][1][0]=U1[2][1][0]*U2[1][1][0]-U1[2][1][1]*U2[1][1][1]
            +U1[2][2][0]*U2[2][1][0]-U1[2][2][1]*U2[2][1][1];
  Q[2][1][1]=U1[2][1][0]*U2[1][1][1]+U1[2][1][1]*U2[1][1][0]
            +U1[2][2][0]*U2[2][1][1]+U1[2][2][1]*U2[2][1][0];
  /* Q_22=U1_21*U2_12+U1_22*U2_22 */
  Q[2][2][0]=U1[2][1][0]*U2[1][2][0]-U1[2][1][1]*U2[1][2][1]
            +U1[2][2][0]*U2[2][2][0]-U1[2][2][1]*U2[2][2][1];
  Q[2][2][1]=U1[2][1][0]*U2[1][2][1]+U1[2][1][1]*U2[1][2][0]
            +U1[2][2][0]*U2[2][2][1]+U1[2][2][1]*U2[2][2][0];

  /* Q=Q*U3 (U3 is block diagonal with U3_00=1, U3_11=1)
     (this changes only third column of Q */
  a=Q[0][2][0]*U3[2][2][0]-Q[0][2][1]*U3[2][2][1];
  b=Q[0][2][0]*U3[2][2][1]+Q[0][2][1]*U3[2][2][0];
  Q[0][2][0]=a; Q[0][2][1]=b;
  a=Q[1][2][0]*U3[2][2][0]-Q[1][2][1]*U3[2][2][1];
  b=Q[1][2][0]*U3[2][2][1]+Q[1][2][1]*U3[2][2][0];
  Q[1][2][0]=a; Q[1][2][1]=b;
  a=Q[2][2][0]*U3[2][2][0]-Q[2][2][1]*U3[2][2][1];
  b=Q[2][2][0]*U3[2][2][1]+Q[2][2][1]*U3[2][2][0];
  Q[2][2][0]=a; Q[2][2][1]=b;

  nflops += 102;

  /* final U=Q*UO3
     (unitary times orthogonal that accumulated Givens rotations) */
#if 0
  U00re=Q[0][0][0]*UO3[0][0]+Q[0][1][0]*UO3[1][0]+Q[0][2][0]*UO3[2][0];
  U00im=Q[0][0][1]*UO3[0][0]+Q[0][1][1]*UO3[1][0]+Q[0][2][1]*UO3[2][0];
  U01re=Q[0][0][0]*UO3[0][1]+Q[0][1][0]*UO3[1][1]+Q[0][2][0]*UO3[2][1];
  U01im=Q[0][0][1]*UO3[0][1]+Q[0][1][1]*UO3[1][1]+Q[0][2][1]*UO3[2][1];
  U02re=Q[0][0][0]*UO3[0][2]+Q[0][1][0]*UO3[1][2]+Q[0][2][0]*UO3[2][2];
  U02im=Q[0][0][1]*UO3[0][2]+Q[0][1][1]*UO3[1][2]+Q[0][2][1]*UO3[2][2];
  U10re=Q[1][0][0]*UO3[0][0]+Q[1][1][0]*UO3[1][0]+Q[1][2][0]*UO3[2][0];
  U10im=Q[1][0][1]*UO3[0][0]+Q[1][1][1]*UO3[1][0]+Q[1][2][1]*UO3[2][0];
  U11re=Q[1][0][0]*UO3[0][1]+Q[1][1][0]*UO3[1][1]+Q[1][2][0]*UO3[2][1];
  U11im=Q[1][0][1]*UO3[0][1]+Q[1][1][1]*UO3[1][1]+Q[1][2][1]*UO3[2][1];
  U12re=Q[1][0][0]*UO3[0][2]+Q[1][1][0]*UO3[1][2]+Q[1][2][0]*UO3[2][2];
  U12im=Q[1][0][1]*UO3[0][2]+Q[1][1][1]*UO3[1][2]+Q[1][2][1]*UO3[2][2];
  U20re=Q[2][0][0]*UO3[0][0]+Q[2][1][0]*UO3[1][0]+Q[2][2][0]*UO3[2][0];
  U20im=Q[2][0][1]*UO3[0][0]+Q[2][1][1]*UO3[1][0]+Q[2][2][1]*UO3[2][0];
  U21re=Q[2][0][0]*UO3[0][1]+Q[2][1][0]*UO3[1][1]+Q[2][2][0]*UO3[2][1];
  U21im=Q[2][0][1]*UO3[0][1]+Q[2][1][1]*UO3[1][1]+Q[2][2][1]*UO3[2][1];
  U22re=Q[2][0][0]*UO3[0][2]+Q[2][1][0]*UO3[1][2]+Q[2][2][0]*UO3[2][2];
  U22im=Q[2][0][1]*UO3[0][2]+Q[2][1][1]*UO3[1][2]+Q[2][2][1]*UO3[2][2];
#endif
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      QLA_Real tr,ti;
      tr = Q[i][0][0]*UO3[0][j]+Q[i][1][0]*UO3[1][j]+Q[i][2][0]*UO3[2][j];
      ti = Q[i][0][1]*UO3[0][j]+Q[i][1][1]*UO3[1][j]+Q[i][2][1]*UO3[2][j];
      QLA_c_eq_r_plus_ir(QLA_elem_M(*U,i,j), tr, ti);
    }
  }

  nflops += 90;

  /* Q=V1*V2 (V1 is block diagonal with V2_11=1,
              V2 is block diagonal with V2_11=1, V2_22=1) */
  Q[0][0][0]=V1[0][0][0]; Q[0][0][1]=V1[0][0][1];
  Q[1][0][0]=V1[1][0][0]; Q[1][0][1]=V1[1][0][1];
  Q[2][0][0]=V1[2][0][0]; Q[2][0][1]=V1[2][0][1];
  Q[0][1][0]=V1[0][1][0]; Q[0][1][1]=V1[0][1][1];
  Q[0][2][0]=V1[0][2][0]; Q[0][2][1]=V1[0][2][1];
  Q[1][1][0]=V1[1][1][0]; Q[1][1][1]=V1[1][1][1];
  Q[2][1][0]=V1[2][1][0]; Q[2][1][1]=V1[2][1][1];
  Q[1][2][0]=V1[1][2][0]*V2[2][2][0]-V1[1][2][1]*V2[2][2][1];
  Q[1][2][1]=V1[1][2][0]*V2[2][2][1]+V1[1][2][1]*V2[2][2][0];
  Q[2][2][0]=V1[2][2][0]*V2[2][2][0]-V1[2][2][1]*V2[2][2][1];
  Q[2][2][1]=V1[2][2][0]*V2[2][2][1]+V1[2][2][1]*V2[2][2][0];

  /* final V=Q*VO3
     (unitary times orthogonal that accumulated Givens rotations) */
#if 0
  V00re=Q[0][0][0]*VO3[0][0]+Q[0][1][0]*VO3[1][0]+Q[0][2][0]*VO3[2][0];
  V00im=Q[0][0][1]*VO3[0][0]+Q[0][1][1]*VO3[1][0]+Q[0][2][1]*VO3[2][0];
  V01re=Q[0][0][0]*VO3[0][1]+Q[0][1][0]*VO3[1][1]+Q[0][2][0]*VO3[2][1];
  V01im=Q[0][0][1]*VO3[0][1]+Q[0][1][1]*VO3[1][1]+Q[0][2][1]*VO3[2][1];
  V02re=Q[0][0][0]*VO3[0][2]+Q[0][1][0]*VO3[1][2]+Q[0][2][0]*VO3[2][2];
  V02im=Q[0][0][1]*VO3[0][2]+Q[0][1][1]*VO3[1][2]+Q[0][2][1]*VO3[2][2];
  V10re=Q[1][0][0]*VO3[0][0]+Q[1][1][0]*VO3[1][0]+Q[1][2][0]*VO3[2][0];
  V10im=Q[1][0][1]*VO3[0][0]+Q[1][1][1]*VO3[1][0]+Q[1][2][1]*VO3[2][0];
  V11re=Q[1][0][0]*VO3[0][1]+Q[1][1][0]*VO3[1][1]+Q[1][2][0]*VO3[2][1];
  V11im=Q[1][0][1]*VO3[0][1]+Q[1][1][1]*VO3[1][1]+Q[1][2][1]*VO3[2][1];
  V12re=Q[1][0][0]*VO3[0][2]+Q[1][1][0]*VO3[1][2]+Q[1][2][0]*VO3[2][2];
  V12im=Q[1][0][1]*VO3[0][2]+Q[1][1][1]*VO3[1][2]+Q[1][2][1]*VO3[2][2];
  V20re=Q[2][0][0]*VO3[0][0]+Q[2][1][0]*VO3[1][0]+Q[2][2][0]*VO3[2][0];
  V20im=Q[2][0][1]*VO3[0][0]+Q[2][1][1]*VO3[1][0]+Q[2][2][1]*VO3[2][0];
  V21re=Q[2][0][0]*VO3[0][1]+Q[2][1][0]*VO3[1][1]+Q[2][2][0]*VO3[2][1];
  V21im=Q[2][0][1]*VO3[0][1]+Q[2][1][1]*VO3[1][1]+Q[2][2][1]*VO3[2][1];
  V22re=Q[2][0][0]*VO3[0][2]+Q[2][1][0]*VO3[1][2]+Q[2][2][0]*VO3[2][2];
  V22im=Q[2][0][1]*VO3[0][2]+Q[2][1][1]*VO3[1][2]+Q[2][2][1]*VO3[2][2];
#endif
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      QLA_Real tr,ti;
      tr = Q[i][0][0]*VO3[0][j]+Q[i][1][0]*VO3[1][j]+Q[i][2][0]*VO3[2][j];
      ti = Q[i][0][1]*VO3[0][j]+Q[i][1][1]*VO3[1][j]+Q[i][2][1]*VO3[2][j];
      QLA_c_eq_r_plus_ir(QLA_elem_M(*V,i,j), tr, ti);
    }
  }

  nflops += 102;

  /* singular values */
  sigma[0]=b00; sigma[1]=b11; sigma[2]=b22;

  info->final_flop += nflops;

  return 0;
}

/* SVD of 2x2 real matrix brought to the form:
    [ a00 a01]
    [   0 a11]
   This routine eliminates off-diagonal element, handling special cases */
static int QOPPC(svd2x2bidiag)(QOP_info_t *info, QLA_Real *a00, QLA_Real *a01, 
			       QLA_Real *a11, QLA_Real U2[2][2], QLA_Real V2[2][2]) {
  register double sinphi, cosphi, tanphi, cotphi;
  register double a, b, min, max, abs00, abs01, abs11;
  register double lna01a11, lna00, ln_num, tau, t;
  register double P00, P01, P10, P11;
  register int isign;
  size_t nflops = 0;

  U2[0][0]=1.; U2[0][1]=0.;
  U2[1][0]=0.; U2[1][1]=1.;
  V2[0][0]=1.; V2[0][1]=0.;
  V2[1][0]=0.; V2[1][1]=1.;

  if( *a00==0 ) {
    if( *a11==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a11)>fabs(*a01) ) {
        cotphi=-(*a01)/(*a11);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a11)/(*a01);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 6;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a01)-sinphi*(*a11);
    (*a01)=0.; (*a11)=0.;
    /* exchange columns in matrix V */
    V2[0][0]=0.; V2[0][1]=1.;
    V2[1][0]=1.; V2[1][1]=0.;
    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 3;
  }
  else if( *a11==0 ) {
    if( *a01==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a01)>fabs(*a00) ) {
        cotphi=-(*a00)/(*a01);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a01)/(*a00);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a00)-sinphi*(*a01);
    (*a01)=0.; (*a11)=0.;
    /* V is just Givens rotation */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;
    nflops += 3;
  }
  else if( *a01==0 ){ /* nothing to be done */
    ;
  }
  else {
    /* need to calculate ( a11^2+a01^2-a00^2 )/( 2*a00*a01 )
       avoiding overflow/underflow,
       use logarithmic coding */
    abs01=fabs(*a01); abs11=fabs(*a11);
    if(abs01>abs11) {
      min=abs11; max=abs01;
    }
    else {
      min=abs01; max=abs11;
    }
    a=min/max;
    lna01a11=2*log(max)+log(1+a*a);

    abs00=fabs(*a00);
    lna00=2*log(abs00);
    if( lna01a11>lna00 ) {
      /* subtract smaller from larger, overall "+" */
      isign=1;
      ln_num=lna01a11+log(1.-exp(lna00-lna01a11));
    }
    else {
      /* subtract larger from smaller, need to change order, overall "-" */
      isign=-1;
      ln_num=lna00+log(1.-exp(lna01a11-lna00));
    }
    a=ln_num-log(2)-log(abs00)-log(abs01);
    tau=exp(a);
    tau*=isign;
    if(*a00<0.)
      {
        tau*=-1;
      }
    if(*a01<0.)
      {
        tau*=-1;
      }

    /* calculate b=sqrt(1+tau^2) */
    a=fabs(tau);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    if(tau>=0.) {
      t = 1.0/(tau + b);
    }
    else {
      t = 1.0/(tau - b);
    }

    /* calculate b=sqrt(1+t^2) */
    a=fabs(t);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    cosphi=1./b;
    sinphi=t*cosphi;

    /* transform matrix A so it has othogonal columns */
    P00= cosphi*(*a00)-sinphi*(*a01);
    P10=-sinphi*(*a11);
    P01= sinphi*(*a00)+cosphi*(*a01);
    P11= cosphi*(*a11);

    /* prepare V  */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;

    /* make column with the largest norm first column */
    if( sqrt(P00*P00+P10*P10)<sqrt(P01*P01+P11*P11) ) {
      a=P00; P00=P01; P01=a;
      a=P10; P10=P11; P11=a;
      /* exchange columns in matrix V2 */
      a=V2[0][0]; V2[0][0]=V2[0][1]; V2[0][1]=a;
      a=V2[1][0]; V2[1][0]=V2[1][1]; V2[1][1]=a;
    }

    /* calculate left Givens rotation and diagonalize */
    if( P10==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(P10)>fabs(P00) ) {
        cotphi=-P00/P10;
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-P10/P00;
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    *a00=P00*cosphi-P10*sinphi;
    *a01=0.;
    *a11=P01*sinphi+P11*cosphi;

    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 56;
  }

  info->final_flop += nflops;

  return 0;
}


