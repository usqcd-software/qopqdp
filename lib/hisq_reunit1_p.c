/**************** hisq_reunit_p.c ******************************/
/*
 * AG: HISQ reunitarization routines. Originally in serial 
 * Fortran code developed by Kit Wong, see
 * A Practical Guide to Lattice Simulations Involving Reunitarized Links
 * K. Y. Wong and R. M. Woloshyn
 * (HPQCD Collaboration)
 * April 2007
 * A.Bazavov 02/2010 Added analytic reunitarization similar to MILC,
 *                   split U(3) and SU(3) projection
 *                   (SU(3) is not supposed to be used so far)
 */

#include <qop_internal.h>
#include <qla_d3.h>
#include <qla_df3.h>

#if 0
static void
u3reunit_site_d(QLA_D3_ColorMatrix *Ur, QLA_D3_ColorMatrix *U)
{
  QLA_D3_ColorMatrix M1, M2;
  QLA_D3_M_eq_Ma_times_M(&M1, U, U);
  QLA_D3_M_eq_sqrt_M(&M2, &M1);
  QLA_D3_M_eq_inverse_M(&M1, &M2);
  QLA_D3_M_eq_M_times_M(Ur, U, &M1);
}

static void
u3reunit_site(QLA_ColorMatrix *Ur, int i, void *args)
{
  QLA_ColorMatrix *U = &((QLA_ColorMatrix *)args)[i];
#if (QLA_Precision=='F')
    QLA_D3_ColorMatrix Ud, Urd;
    QLA_DF3_M_eq_M(&Ud, U);
    u3reunit_site_d(&Urd, &Ud);
    QLA_FD3_M_eq_M(Ur, &Urd);
#else
    u3reunit_site_d(Ur, U);
#endif
}

// requires QLA >= 1.7.0 and QDP >= 1.9.0
void
QOPPC(u3reunit)(QOP_into_t *info, QDP_ColorMatrix *U, QDP_ColorMatrix *Ur)
{
  QLA_ColorMatrix *Uq = QDP_expose_M(U);
  QDP_M_eq_funcia(Ur, u3reunit_site, (void *)Uq, QDP_all);
  QDP_reset_M(U);
}
#endif

//AB output for debugging
// static void PrintT4( QLA_ColorTensor4 dwdv ) {
//   int i,j,k,l;
//   i=j=k=l=0;
// 
//   for(i=0;i<3;i++) {
//     for(j=0;j<3;j++) {
//       for(k=0;k<3;k++) {
//         for(l=0;l<3;l++) {
//           printf("T4[%d][%d][%d][%d].real=%lf  T4[%d][%d][%d][%d].imag=%lf\n", 
//             i,j,k,l,dwdv.t4[i][j][k][l].real,
//             i,j,k,l,dwdv.t4[i][j][k][l].imag );
//         }
//       }
//     }
//   }
// }
// 
// static void 
// PrintMone(QLA_ColorMatrix m){
//   printf("Field: %e %e %e\n", m.e[0][0].real,
// 	 m.e[0][1].real, m.e[0][2].real);
//   printf("       %e %e %e\n", m.e[1][0].real,
// 	 m.e[1][1].real, m.e[1][2].real);
//   printf("       %e %e %e\n\n", m.e[2][0].real,
// 	 m.e[2][1].real, m.e[2][2].real);
// }

//#define QOPQDP_REUNIT_FOLLOW_PRECISION

//AB Define precision for reunitarization
#ifdef QOPQDP_REUNIT_FOLLOW_PRECISION
#  define custom_QLA_ColorMatrix QLA_ColorMatrix
#  define custom_QLA_Real QLA_Real
#else
#  define custom_QLA_ColorMatrix QLA_D3_ColorMatrix
#  define custom_QLA_Real QLA_D_Real
#endif

//AB Reunitarization to U(3), similar to the MILC code,
//   this is single link operation, it can be done outside
//   QDP, it is not field based operation.
void
QOP_u3reunit(QOP_info_t *info, QDP_ColorMatrix *V, QDP_ColorMatrix *W)
{
  int svd_calls = 0;
  double dtime = -QOP_time();
  // stop QDP operations and expose QDP links
  QLA_ColorMatrix *Vlinks = QDP_expose_M(V);
  QLA_ColorMatrix *Wlinks = QDP_expose_M(W);
  info->final_flop = 0.0;

#pragma omp parallel
  {
    int svd_callst = 0;
    QOP_info_t infot;
    infot.final_flop = 0.0;
#pragma omp for
    for(int i=0; i<QDP_sites_on_node; i++) {
      svd_callst += QOP_u3_un_analytic(&infot, &(Vlinks[i]), &(Wlinks[i]));
    }
#pragma omp critical
    {
      info->final_flop += infot.final_flop;
      svd_calls += svd_callst;
    }
  }

  // Tally across all nodes and update global counter
  QMP_sum_int(&svd_calls);
  QOP_info_hisq_svd_counter(info) += svd_calls;

  // resume QDP operations on links
  QDP_reset_M(V);
  QDP_reset_M(W);
  dtime += QOP_time();
  info->final_sec = dtime;
}

// not yet supported
#if 0
//Reuniterise links to SU3 (see eq. 53, 54 in the practical guide) 
void QOPPC(su3reunit)(QOP_info_t *info, QDP_ColorMatrix *U,QDP_ColorMatrix *Ur)
{
  QDP_ColorMatrix *msquare,*mpolytot,*mtemp1,*mtemp2, *Ur_temp;
  QDP_Complex *complex_temp;
  QDP_Complex *det;
  QDP_Real *real_temp, *im_temp, *real_temp1, *im_temp1;
  QDP_Real *poly_temp, *z_norm, *z_arg, *re_det, *im_det;
  QDP_Real *qdp_const;
  QDP_Int *logical_temp;
  QLA_Real qla_const;
  const int npolyterm=14;
  QLA_Real poly_c[npolyterm+1], poly_d[npolyterm+1];
  int j,idiag,ic,jc;

// Definition for the fraction:
  poly_c[0] = 0.085091;
  poly_c[1] = 2.413330975e+00;
  poly_c[2] = 6.257184884e-01;
  poly_c[3] = 2.925707925e-01;
  poly_c[4] = 1.737405612e-01;
  poly_c[5] = 1.166359792e-01;
  poly_c[6] = 8.372555094e-02;
  poly_c[7] = 6.216038074e-02;
  poly_c[8] = 4.652496186e-02;
  poly_c[9] = 3.423610040e-02;
  poly_c[10] = 2.404754621e-02;
  poly_c[11] = 1.545550091e-02;
  poly_c[12] = 8.436481876e-03;
  poly_c[13] = 3.419245947e-03;
  poly_c[14] = 1.138166539e-03;

  poly_d[1] = 1.361747338e+01;
  poly_d[2] = 3.135687028e+00;
  poly_d[3] = 1.213113539e+00;
  poly_d[4] = 5.596349298e-01;
  poly_d[5] = 2.752627333e-01;
  poly_d[6] = 1.364115846e-01;
  poly_d[7] = 6.543005714e-02;
  poly_d[8] = 2.923946484e-02;
  poly_d[9] = 1.164228894e-02;
  poly_d[10] = 3.887745892e-03;
  poly_d[11] = 9.937321442e-04;
  poly_d[12] = 1.684882417e-04;
  poly_d[13] = 1.585925699e-05;
  poly_d[14] = 5.914114023e-07;

  msquare = QDP_create_M();
  mpolytot = QDP_create_M();
  mtemp1 = QDP_create_M();
  mtemp2 = QDP_create_M();
  Ur_temp = QDP_create_M();
  complex_temp=QDP_create_C();
  det=QDP_create_C();
  real_temp=QDP_create_R();
  im_temp=QDP_create_R();
  real_temp1=QDP_create_R();
  im_temp1=QDP_create_R();
  poly_temp=QDP_create_R();
  re_det=QDP_create_R();
  im_det=QDP_create_R();
  z_norm=QDP_create_R();
  z_arg=QDP_create_R();
  qdp_const=QDP_create_R();
  logical_temp=QDP_create_I();

  //The reuniterisation
  QDP_M_eq_zero(mpolytot,QDP_all);
  
  QDP_M_eq_Ma_times_M(msquare,U,U,QDP_all);
  
  for(j=1;j<=npolyterm;j++) {

    QDP_M_eq_M(mtemp1,msquare,QDP_all);
    
    //convert constant to QDP Real field
    QDP_R_eq_r(poly_temp,&poly_d[j],QDP_all);
    
    //calculate mpolytot
    for(idiag=0;idiag<3;idiag++){
      //extract complex color matrix element
      QDP_C_eq_elem_M(complex_temp,mtemp1,idiag,idiag,QDP_all);
      //extract real part
      QDP_R_eq_re_C(real_temp,complex_temp,QDP_all);
      //extract imaginary part
      QDP_R_eq_im_C(im_temp,complex_temp,QDP_all);
      //update real part
      QDP_R_peq_R(real_temp,poly_temp,QDP_all);
      //insert real,im parts back to complex
      QDP_C_eq_R_plus_i_R(complex_temp,real_temp,im_temp,QDP_all);
      //insert complex back to color matrix
      QDP_M_eq_elem_C(mtemp1,complex_temp,idiag,idiag,QDP_all);  
    }
    
    QOPPC(su3inverse)(mtemp1,mtemp2);
    
    QDP_M_peq_r_times_M(mpolytot,&poly_c[j],mtemp2,QDP_all);
  }
  
  
  
  //convert constant to QDP Real field
  QDP_R_eq_r(poly_temp,&poly_c[0],QDP_all);
  
  for(idiag=0;idiag<3;idiag++){
    //extract complex color matrix element
    QDP_C_eq_elem_M(complex_temp,mpolytot,idiag,idiag,QDP_all);
    //extract real part
    QDP_R_eq_re_C(real_temp,complex_temp,QDP_all);
    //extract imaginary part
    QDP_R_eq_im_C(im_temp,complex_temp,QDP_all);
    //update real part
    QDP_R_peq_R(real_temp,poly_temp,QDP_all);
    //insert real,im parts back to complex
    QDP_C_eq_R_plus_i_R(complex_temp,real_temp,im_temp,QDP_all);
    //insert complex back to color matrix
    QDP_M_eq_elem_C(mpolytot,complex_temp,idiag,idiag,QDP_all);  
  }
  

  //calculate U times mpolytot and det of this
  QDP_M_eq_M_times_M(mtemp1,U,mpolytot,QDP_all);
  QOPPC(su3det)(mtemp1,det);
  
  //extract re and im parts of det
  QDP_R_eq_re_C(re_det,det,QDP_all);
  QDP_R_eq_im_C(im_det,det,QDP_all);


//AB QUICK FIX TO REMOVE PROJECTION TO SU(3)
  qla_const=1.;
  QDP_R_eq_r(re_det,&qla_const,QDP_all);
  qla_const=0.;
  QDP_R_eq_r(im_det,&qla_const,QDP_all);
  
  //normalise det
  qla_const=0.166666667;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_R_eq_R_times_R(real_temp,re_det,re_det,QDP_all);
  QDP_R_peq_R_times_R(real_temp,im_det,im_det,QDP_all);
  QDP_R_eq_R_pow_R(z_norm,real_temp,qdp_const,QDP_all);
  
  QDP_R_eq_R_divide_R(real_temp,im_det,re_det,QDP_all);
  QDP_R_eq_atan_R(z_arg,real_temp,QDP_all);
  
  
  qla_const=3.141592654;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  QDP_R_eq_R_plus_R(real_temp,z_arg,qdp_const,QDP_all);
  
  
  qla_const=0.;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_I_eq_R_lt_R(logical_temp,re_det,qdp_const,QDP_all);     
  QDP_R_eq_R_mask_I(z_arg,real_temp,logical_temp,QDP_all);
  
  
  qla_const=3.;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_R_eq_R_divide_R(real_temp,z_arg,qdp_const,QDP_all);
  QDP_R_eq_R(z_arg, real_temp,QDP_all);
  
  QDP_R_eq_cos_R(real_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(re_det,z_norm,real_temp,QDP_all);
  
  QDP_R_eq_sin_R(im_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(im_det,z_norm,im_temp,QDP_all);
  
  
  //calculate z_norm
  qla_const=1.;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_R_eq_R_times_R(real_temp,re_det,re_det,QDP_all);
  QDP_R_peq_R_times_R(real_temp,im_det,im_det,QDP_all);
  QDP_R_eq_R_divide_R(z_norm,qdp_const,real_temp,QDP_all);
  
  
  //do det multiplication and normalise by z_norm
  for (ic=0;ic<3;ic++){
    for (jc=0;jc<3;jc++){
      
      //extract complex color matrix elements of mtemp1
      QDP_C_eq_elem_M(complex_temp,mtemp1,ic,jc,QDP_all);
      //extract real part of mtemp1
      QDP_R_eq_re_C(real_temp,complex_temp,QDP_all);
      //extract im part of mtemp1
      QDP_R_eq_im_C(im_temp,complex_temp,QDP_all);
      //calculate real part
      QDP_R_eq_R_times_R(real_temp1,real_temp,re_det,QDP_all);
      QDP_R_peq_R_times_R(real_temp1,im_temp,im_det,QDP_all);
      //calculate imaginary part
      QDP_R_eq_R_times_R(im_temp1,im_temp,re_det,QDP_all);
      QDP_R_meq_R_times_R(im_temp1,real_temp,im_det,QDP_all);
      //multiply by z_norm
      QDP_R_eq_R_times_R(real_temp,real_temp1,z_norm,QDP_all);
      QDP_R_eq_R_times_R(im_temp,im_temp1,z_norm,QDP_all);
      //insert real,im parts back to complex
      QDP_C_eq_R_plus_i_R(complex_temp,real_temp,im_temp,QDP_all);
      //insert complex to color matrix
      QDP_M_eq_elem_C(Ur,complex_temp,ic,jc,QDP_all); 

    }
  }
 
  QDP_destroy_M(msquare);
  QDP_destroy_M(mpolytot);
  QDP_destroy_M(mtemp1);
  QDP_destroy_M(mtemp2);
  QDP_destroy_M(Ur_temp);
  QDP_destroy_C(complex_temp);
  QDP_destroy_C(det);
  QDP_destroy_R(real_temp);
  QDP_destroy_R(im_temp);
  QDP_destroy_R(real_temp1);
  QDP_destroy_R(im_temp1);
  QDP_destroy_R(poly_temp);
  QDP_destroy_R(re_det);
  QDP_destroy_R(im_det);
  QDP_destroy_R(z_norm);
  QDP_destroy_R(z_arg);
  QDP_destroy_R(qdp_const);
  QDP_destroy_I(logical_temp);

  return;
}
#endif

void
QOP_hisq_force_multi_reunit(QOP_info_t *info,
			    QDP_ColorMatrix *V[4],
			    QDP_ColorMatrix *Force[4],
			    QDP_ColorMatrix *Force_old[4])
{
  int svd_calls = 0;
  int ff_counter = 0;
  info->final_flop = 0.0;
  int nd = QDP_ndim();

  for(int i=0; i<nd; i++ ) {
    // stop QDP operations and expose QDP links
    QLA_ColorMatrix *Vlinks = QDP_expose_M( V[i] );
    QLA_ColorMatrix *ff_new = QDP_expose_M( Force[i] );
    QLA_ColorMatrix *ff_old = QDP_expose_M( Force_old[i] );
    //AB NO NEED TO REPHASE V LINKS???
#pragma omp parallel
    {
      int svd_callst = 0;
      int ff_countert = 0;
      QOP_info_t infot;
      infot.final_flop = 0.0;
#pragma omp for
      for(int j=0; j<QDP_sites_on_node; j++ ) {
	QLA_ColorTensor4 dwdv, dwdagdv;
	// derivative with respect to V and V^+
	QOPPC(u3_un_der_analytic)( &infot, &( Vlinks[j] ),
				   &dwdv, &dwdagdv, &svd_calls, &ff_countert );
	// adjoint piece of force from the previous level
	QLA_ColorMatrix ff_old_adj;
	QLA_M_eq_Ma( &ff_old_adj, &( ff_old[j] ) );
	// see LONG COMMENT in fermion_force_fn_multi_hisq
	QLA_M_eq_zero( &( ff_new[j] ) );
	for(int m=0; m<3; m++) {
	  for(int n=0; n<3; n++) {
	    for(int k=0; k<3; k++) {
	      for(int l=0; l<3; l++) {
		// direct part
		QLA_Complex ftmp;
		QLA_c_eq_c_times_c( ftmp, dwdv.t4[k][m][n][l], QLA_elem_M(ff_old[j],l,k) );
		QLA_c_peq_c( QLA_elem_M(ff_new[j],n,m), ftmp );
		// adjoint part
		QLA_c_eq_c_times_c( ftmp, dwdagdv.t4[k][m][n][l], QLA_elem_M(ff_old_adj,l,k) );
		QLA_c_peq_c( QLA_elem_M(ff_new[j],n,m), ftmp );
	      }
	    }
	  }
	}
#pragma omp critical
	{
	  info->final_flop += infot.final_flop;
	  svd_calls += svd_callst;
	  ff_counter += ff_countert;
	}
      }
    }
    // resume QDP operations on links
    QDP_reset_M( V[i] );
    QDP_reset_M( Force[i] );
    QDP_reset_M( Force_old[i] );
  }

  // Tally across all nodes and update global counter
  QMP_sum_int(&svd_calls);
  QMP_sum_int(&ff_counter);
  QOP_info_hisq_svd_counter(info) += svd_calls;
  QOP_info_hisq_force_filter_counter(info) += ff_counter;
  // not sure why, but the accumulated flop count seems too large
  info->final_flop = ((double)(198 + 81*16))*QDP_sites_on_node;
}
