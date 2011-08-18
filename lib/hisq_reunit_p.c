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

// This alternate version requires some thought before adoption.
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
QOPPC(u3reunit)(QDP_ColorMatrix *U, QDP_ColorMatrix *Ur)
{
  QLA_ColorMatrix *Uq = QDP_expose_M(U);
  QDP_M_eq_funcia(Ur, u3reunit_site, (void *)Uq, QDP_all);
  QDP_reset_M(U);
}
#endif

//utilities
static void QOPPC(su3det)(QDP_ColorMatrix *mat,QDP_Complex *det);
static void QOPPC(su3inverse)(QDP_ColorMatrix *mat,QDP_ColorMatrix *inv);
//static void QOPPC(su3adjoint)(QDP_ColorMatrix *mat,QDP_ColorMatrix *adj);



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

#define custom_QLA_ColorMatrix QLA_ColorMatrix
#define custom_QLA_Real QLA_Real

#else

#define custom_QLA_ColorMatrix QLA_D3_ColorMatrix
#define custom_QLA_Real QLA_D_Real

#endif

//AB Reunitarization to U(3), similar to the MILC code,
//   this is single link operation, it can be done outside
//   QDP, it is not field based operation.
//   Depending on QOPQDP_REUNIT_FOLLOW_PRECISION
//   reunitarization can be done either in dominant precision
//   or always in double precision
void QOPPC(u3reunit)(QDP_ColorMatrix *V, QDP_ColorMatrix *W)
{
  QLA_ColorMatrix *Vlinks;
  QLA_ColorMatrix *Wlinks;
  
  //custom_QLA_ColorMatrix A;
  
  int i;
  
//  printf("sizeof(custom_QLA_ColorMatrix)=%d\n",sizeof(custom_QLA_ColorMatrix));

  // stop QDP operations and expose QDP links
  Vlinks = QDP_expose_M(V);
  Wlinks = QDP_expose_M(W);

  
  for( i=0; i<QDP_sites_on_node; i++ ) {
    QOPPC(su3_un_analytic)( &(Vlinks[i]), &(Wlinks[i]) );
  }
  
  // resume QDP operations on links
  QDP_reset_M( V );
  QDP_reset_M( W );
}


//Reuniterise links to SU3 (see eq. 53, 54 in the practical guide) 
void QOPPC(su3reunit)(QDP_ColorMatrix *U,QDP_ColorMatrix *Ur)
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
  
  
  for(j=1;j<=npolyterm;j++){
    
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




void
QOPPC(hisq_force_multi_reunit)(QOP_info_t *info,
			       QDP_ColorMatrix *V[4],
			       QDP_ColorMatrix *Force[4],
			       QDP_ColorMatrix *Force_old[4])
{
  int i, j, k, l, m, n, nd;
  QLA_ColorMatrix *Vlinks;
  QLA_ColorMatrix *ff_new, *ff_old, ff_old_adj;
  QLA_Complex ftmp;
  QLA_ColorTensor4 dwdv, dwdagdv;

  nd = QDP_ndim();

  for( i=0; i<nd; i++ ) {
    // stop QDP operations and expose QDP links
    Vlinks = QDP_expose_M( V[i] );
    ff_new = QDP_expose_M( Force[i] );
    ff_old = QDP_expose_M( Force_old[i] );

    //AB NO NEED TO REPHASE V LINKS???
  
  
    for( j=0; j<QDP_sites_on_node; j++ ) {
      // derivative with respect to V and V^+
      QOPPC(su3_un_der_analytic)( &( Vlinks[j] ),
            &dwdv, &dwdagdv );
    
//      if( j==0 ) {
//        printf("FORCE REUNIT V link at site 0\n");
//        PrintMone( Vlinks[j] );
//        printf("FORCE REUNIT derivative at site 0\n");
//        PrintT4(dwdv);
//        PrintT4(dwdagdv);
//      }

      // adjoint piece of force from the previous level
      QLA_M_eq_Ma( &ff_old_adj, &( ff_old[j] ) );

      // see LONG COMMENT in fermion_force_fn_multi_hisq
      QLA_M_eq_zero( &( ff_new[j] ) );
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          for( k=0; k<3; k++) {
            for( l=0; l<3; l++) {
              // direct part
              QLA_c_eq_c_times_c( ftmp, dwdv.t4[k][m][n][l], QLA_elem_M(ff_old[j],l,k) );
              QLA_c_peq_c( QLA_elem_M(ff_new[j],n,m), ftmp );
              // adjoint part
              QLA_c_eq_c_times_c( ftmp, dwdagdv.t4[k][m][n][l], QLA_elem_M(ff_old_adj,l,k) );
              QLA_c_peq_c( QLA_elem_M(ff_new[j],n,m), ftmp );
            }
          }
        }
      }
//      if( j==0 ) {
//        printf("FORCE REUNIT old force\n");
//        PrintMone( ff_old[j] );
//        printf("FORCE REUNIT new force\n");
//        PrintMone( ff_new[j] );
//      }
    }


    // resume QDP operations on links
    QDP_reset_M( V[i] );
    QDP_reset_M( Force[i] );
    QDP_reset_M( Force_old[i] );
  }
}


#ifdef ALAN_GRAY_REUNIT_FORCE
//compute derivatives (see eq. 55, 56 in the practical guide) 
void
QOPPC(hisq_force_multi_reunit)(QOP_info_t *info,
			       QDP_ColorMatrix *u[4],
			       QDP_ColorMatrix *f1[4],
			       QDP_ColorMatrix *f0[4])
{

  //variable declarations

  int j, id, irow, icol, jrow, jcol, idiag,ipoly;

  const int npolyterm=14;
  QLA_Real poly_c[npolyterm+1], poly_d[npolyterm+1];

  //fields
  QDP_ColorMatrix *msquare;

  QDP_ColorMatrix *sumterm;
  QDP_Real *re_sumterm,*im_sumterm;

  QDP_ColorMatrix *mtemp1;
  QDP_Real *re_mtemp1,*im_mtemp1;

  QDP_ColorMatrix *term[npolyterm+1];
  QDP_Real *re_term,*im_term;

  QDP_ColorMatrix *Sterm[npolyterm+1];
  QDP_Real *re_Sterm1,*im_Sterm1;
  QDP_Real *re_Sterm2,*im_Sterm2;

  QDP_ColorMatrix *StermSd[npolyterm+1];
  QDP_Real *re_StermSd,*im_StermSd;

  QDP_ColorMatrix *mW; 
  QDP_Real *re_mW,*im_mW;
 
  QDP_ColorMatrix *mWadj; 

  QDP_ColorMatrix *dWdsigma; 
  QDP_Real *re_dWdsigma,*im_dWdsigma;

  QDP_ColorMatrix *dWdsigmad; 
  QDP_Real *re_dWdsigmad,*im_dWdsigmad;
 
  QDP_ColorMatrix *durdsigma; 
  QDP_Real *re_durdsigma,*im_durdsigma;

  QDP_ColorMatrix *durdsigmad; 
  QDP_Real *re_durdsigmad,*im_durdsigmad;

  QDP_ColorMatrix *Tlocal; 
  QDP_Real *re_Tlocal,*im_Tlocal;

  QDP_ColorMatrix *Tplocal; 
  QDP_Real *re_Tplocal,*im_Tplocal;

  QDP_Real *re_f0,*im_f0;
 
  QDP_Real *real_temp, *real_temp1, *real_temp2;
  QDP_Real *poly_temp, *z_norm, *z_arg;
  QDP_Real *re_ctemp,*im_ctemp;
  QDP_Real *re_det, *im_det,*re_det13, *im_det13, *re_det43, *im_det43;

  QDP_Complex *complex_temp, *det;
  
  QLA_Real qla_const;
  QDP_Real *qdp_const;

  QDP_Int *logical_temp;

  
  //Allocate fields
  msquare = QDP_create_M();
  mWadj = QDP_create_M();

  mtemp1 = QDP_create_M();
  re_mtemp1=QDP_create_R(); im_mtemp1=QDP_create_R();

  dWdsigma = QDP_create_M();
  re_dWdsigma=QDP_create_R(); im_dWdsigma=QDP_create_R();

  dWdsigmad = QDP_create_M();
  re_dWdsigmad=QDP_create_R(); im_dWdsigmad=QDP_create_R();

  mW = QDP_create_M();
  re_mW=QDP_create_R(); im_mW=QDP_create_R();

  durdsigma = QDP_create_M();
  re_durdsigma=QDP_create_R(); im_durdsigma=QDP_create_R();

  durdsigmad = QDP_create_M();
  re_durdsigmad=QDP_create_R(); im_durdsigmad=QDP_create_R();

  Tlocal = QDP_create_M();
  re_Tlocal=QDP_create_R(); im_Tlocal=QDP_create_R();

  Tplocal = QDP_create_M();
  re_Tplocal=QDP_create_R(); im_Tplocal=QDP_create_R();

  sumterm = QDP_create_M();
  re_sumterm=QDP_create_R(); im_sumterm=QDP_create_R();

  re_f0=QDP_create_R(); im_f0=QDP_create_R();

  for(j=0;j<=npolyterm;j++) Sterm[j] = QDP_create_M();
  re_Sterm1=QDP_create_R(); im_Sterm1=QDP_create_R();
  re_Sterm2=QDP_create_R(); im_Sterm2=QDP_create_R();

  for(j=0;j<=npolyterm;j++) StermSd[j] = QDP_create_M();
  re_StermSd=QDP_create_R(); im_StermSd=QDP_create_R();

  for(j=0;j<=npolyterm;j++) term[j] = QDP_create_M();
  re_term=QDP_create_R(); im_term=QDP_create_R();

  logical_temp=QDP_create_I();
  complex_temp=QDP_create_C();

  det=QDP_create_C();
  re_det=QDP_create_R();
  im_det=QDP_create_R();
  re_det13=QDP_create_R();
  im_det13=QDP_create_R();
  re_det43=QDP_create_R();
  im_det43=QDP_create_R();
  real_temp=QDP_create_R();
  real_temp1=QDP_create_R();
  real_temp2=QDP_create_R();
  poly_temp=QDP_create_R();
  z_norm=QDP_create_R();
  z_arg=QDP_create_R();
  qdp_const=QDP_create_R();
  re_ctemp=QDP_create_R();
  im_ctemp=QDP_create_R();


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


  //main body of routine
 for(id=0;id<4;id++)
    {

 // Preparation:
 QDP_M_eq_zero(sumterm,QDP_all);
  
 QDP_M_eq_Ma_times_M(msquare,u[id],u[id],QDP_all);

 for(j=1;j<=npolyterm;j++){
    
    QDP_M_eq_M(mtemp1,msquare,QDP_all);
    
    //convert constant to QDP Real field
    QDP_R_eq_r(poly_temp,&poly_d[j],QDP_all);
    
    //calculate mtemp1
    for(idiag=0;idiag<3;idiag++){
      //extract complex color matrix element
      QDP_C_eq_elem_M(complex_temp,mtemp1,idiag,idiag,QDP_all);
      //extract real part
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      //extract imaginary part
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      //update real part
      QDP_R_peq_R(re_mtemp1,poly_temp,QDP_all);
      //insert real,im parts back to complex
      QDP_C_eq_R_plus_i_R(complex_temp,re_mtemp1,im_mtemp1,QDP_all);
      //insert complex back to color matrix
      QDP_M_eq_elem_C(mtemp1,complex_temp,idiag,idiag,QDP_all);  
    }
    
    QOPPC(su3inverse)(mtemp1,term[j]);
    
    QDP_M_peq_r_times_M(sumterm,&poly_c[j],term[j],QDP_all);

    QDP_M_eq_M_times_M(Sterm[j],u[id],term[j],QDP_all);
    QDP_M_eq_M_times_Ma(StermSd[j],Sterm[j],u[id],QDP_all);
    
  }


 // We also need W itself and its determinants:

    QDP_M_eq_M(mtemp1,sumterm,QDP_all);

   //convert constant to QDP Real field
    QDP_R_eq_r(poly_temp,&poly_c[0],QDP_all);
    
  
  for(idiag=0;idiag<3;idiag++){
      //extract complex color matrix element
      QDP_C_eq_elem_M(complex_temp,mtemp1,idiag,idiag,QDP_all);
      //extract real part
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      //extract imaginary part
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      //update real part
      QDP_R_peq_R(re_mtemp1,poly_temp,QDP_all);
      //insert real,im parts back to complex
      QDP_C_eq_R_plus_i_R(complex_temp,re_mtemp1,im_mtemp1,QDP_all);
      //insert complex back to color matrix
      QDP_M_eq_elem_C(mtemp1,complex_temp,idiag,idiag,QDP_all);  
    }

  QDP_M_eq_M_times_M(mW,u[id],mtemp1,QDP_all);

  QOPPC(su3adjoint)(mW,mWadj);

  //calculate determinant
  QOPPC(su3det)(mW,det);
  
  //extract re and im parts of det
  QDP_R_eq_re_C(re_det,det,QDP_all);
  QDP_R_eq_im_C(im_det,det,QDP_all);
  
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
  

  //det 13
  qla_const=3.;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_R_eq_R_divide_R(real_temp,z_arg,qdp_const,QDP_all);
  QDP_R_eq_R(z_arg, real_temp,QDP_all);

  QDP_R_eq_cos_R(real_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(re_det13,z_norm,real_temp,QDP_all);
  
  QDP_R_eq_sin_R(real_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(im_det13,z_norm,real_temp,QDP_all);


  //det43
  qla_const=4.;
  QDP_R_eq_r(qdp_const,&qla_const,QDP_all);
  
  QDP_R_eq_R_times_R(real_temp,z_arg,qdp_const,QDP_all);
  QDP_R_eq_R(z_arg, real_temp,QDP_all);


  ///z_norm=z_norm**4
  QDP_R_eq_R_times_R(real_temp,z_norm,z_norm,QDP_all);
  QDP_R_eq_R_times_R(z_norm,real_temp,real_temp,QDP_all);
  
  QDP_R_eq_cos_R(real_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(re_det43,z_norm,real_temp,QDP_all);
  
  QDP_R_eq_sin_R(real_temp,z_arg,QDP_all);
  QDP_R_eq_R_times_R(im_det43,z_norm,real_temp,QDP_all);
 

 //!*==================
 //!* Obtaining T & Tp
 //!*==================

 QDP_M_eq_zero(Tlocal,QDP_all);
 QDP_M_eq_zero(Tplocal,QDP_all);

  for(icol=0;icol<3;icol++){
    for(irow=0;irow<3;irow++){
      
 
      //!*-----------
      //!* durdsigma
      //!*-----------
      //!* poly_c(0) term:

      QDP_M_eq_zero(dWdsigma,QDP_all);

      qla_const=0.;
      QDP_R_eq_r(re_dWdsigma,&poly_c[0],QDP_all);
      QDP_R_eq_r(im_dWdsigma,&qla_const,QDP_all);
      //insert real,im parts back to complex
      QDP_C_eq_R_plus_i_R(complex_temp,re_dWdsigma,im_dWdsigma,QDP_all);
      //insert complex back to color matrix
      QDP_M_eq_elem_C(dWdsigma,complex_temp,irow,icol,QDP_all);


      //dS/dS term in the product rule:
      for (jcol=0;jcol<3;jcol++)
	{

	  QDP_C_eq_elem_M(complex_temp,dWdsigma,irow,jcol,QDP_all);
	  QDP_R_eq_re_C(re_dWdsigma,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_dWdsigma,complex_temp,QDP_all);
	  
	  QDP_C_eq_elem_M(complex_temp,sumterm,icol,jcol,QDP_all);
	  QDP_R_eq_re_C(re_sumterm,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_sumterm,complex_temp,QDP_all);
	  
	  QDP_R_peq_R(re_dWdsigma,re_sumterm,QDP_all);
	  QDP_R_peq_R(im_dWdsigma,im_sumterm,QDP_all);
	  
	  QDP_C_eq_R_plus_i_R(complex_temp,re_dWdsigma,im_dWdsigma,QDP_all);
	  QDP_M_eq_elem_C(dWdsigma,complex_temp,irow,jcol,QDP_all);  

	}

      //!* The messy term in the product rule:

      for (ipoly=1;ipoly<=npolyterm;ipoly++)
	{
	QDP_R_eq_r(poly_temp,&poly_c[ipoly],QDP_all);

	for (jcol=0;jcol<3;jcol++)
	  {
	  for (jrow=0;jrow<3;jrow++)
		{

	  QDP_C_eq_elem_M(complex_temp,dWdsigma,jrow,jcol,QDP_all);
	  QDP_R_eq_re_C(re_dWdsigma,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_dWdsigma,complex_temp,QDP_all);

	  QDP_C_eq_elem_M(complex_temp,StermSd[ipoly],jrow,irow,QDP_all);
	  	  QDP_R_eq_re_C(re_StermSd,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_StermSd,complex_temp,QDP_all);

	  QDP_C_eq_elem_M(complex_temp,term[ipoly],icol,jcol,QDP_all);
	  QDP_R_eq_re_C(re_term,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_term,complex_temp,QDP_all);

	  //real part
	  QDP_R_eq_R_times_R(real_temp,re_StermSd,re_term,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_StermSd,im_term,QDP_all);
	  QDP_R_eq_R_minus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_meq_R_times_R(re_dWdsigma,poly_temp,real_temp2,QDP_all);
			     
	  //im part
	  QDP_R_eq_R_times_R(real_temp,re_StermSd,im_term,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_StermSd,re_term,QDP_all);
	  QDP_R_eq_R_plus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_meq_R_times_R(im_dWdsigma,poly_temp,real_temp2,QDP_all);

	  QDP_C_eq_R_plus_i_R(complex_temp,re_dWdsigma,im_dWdsigma,QDP_all);
	  QDP_M_eq_elem_C(dWdsigma,complex_temp,jrow,jcol,QDP_all);  

		}
	      }
	    }


      //!* Finally durdsigma:

      //calculate ctemp
      QDP_M_eq_M_times_M(mtemp1,mWadj,dWdsigma,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,0,0,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_eq_R(im_ctemp,im_mtemp1,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,1,1,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_peq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_peq_R(im_ctemp,im_mtemp1,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,2,2,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_peq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_peq_R(im_ctemp,im_mtemp1,QDP_all);

      //update ctemp real part
      QDP_R_eq_R_times_R(real_temp1,re_ctemp,re_det43,QDP_all);
      QDP_R_eq_R_times_R(real_temp2,im_ctemp,im_det43,QDP_all);
      QDP_R_eq_R_plus_R(re_ctemp,real_temp1,real_temp2,QDP_all);
      
      //update ctemp im part
      QDP_R_eq_R_times_R(real_temp1,im_ctemp,re_det43,QDP_all);
      QDP_R_eq_R_times_R(real_temp2,re_ctemp,im_det43,QDP_all);
      QDP_R_eq_R_minus_R(im_ctemp,real_temp1,real_temp2,QDP_all);

      //znorm
      QDP_R_eq_R_times_R(z_norm,re_det43,re_det43,QDP_all);
      QDP_R_peq_R_times_R(z_norm,im_det43,im_det43,QDP_all);

      //update ctemp
      qla_const=3.;
      QDP_R_eq_r_times_R(real_temp,&qla_const,z_norm,QDP_all);

      QDP_R_eq_R_divide_R(real_temp2,re_ctemp,real_temp,QDP_all);
      QDP_R_eq_R(re_ctemp,real_temp2,QDP_all);

      QDP_R_eq_R_divide_R(real_temp2,im_ctemp,real_temp,QDP_all);
      QDP_R_eq_R(im_ctemp,real_temp2,QDP_all);

      for (jcol=0;jcol<3;jcol++)
	{
	  for (jrow=0;jrow<3;jrow++)
	    {
	  
	      QDP_C_eq_elem_M(complex_temp,dWdsigma,jrow,jcol,QDP_all);
	      QDP_R_eq_re_C(re_dWdsigma,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_dWdsigma,complex_temp,QDP_all);

	      QDP_R_eq_R_times_R(real_temp1,re_dWdsigma,re_det13,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,im_dWdsigma,im_det13,QDP_all);
	      QDP_R_eq_R_plus_R(re_durdsigma,real_temp1,real_temp2,QDP_all);
	      
	      QDP_R_eq_R_times_R(real_temp1,im_dWdsigma,re_det13,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,re_dWdsigma,im_det13,QDP_all);
	      QDP_R_eq_R_minus_R(im_durdsigma,real_temp1,real_temp2,QDP_all);
	      QDP_C_eq_R_plus_i_R(complex_temp,re_durdsigma,im_durdsigma,
				  QDP_all);
	      QDP_M_eq_elem_C(durdsigma,complex_temp,jrow,jcol,QDP_all);

		}
	      }

      //znorm
      QDP_R_eq_R_times_R(real_temp,re_det13,re_det13,QDP_all);
      QDP_R_peq_R_times_R(real_temp,im_det13,im_det13,QDP_all);
      
      qla_const=1.;
      QDP_R_eq_r(real_temp1,&qla_const,QDP_all);
      QDP_R_eq_R_divide_R(z_norm,real_temp1,real_temp,QDP_all);
      
      //update durdsigma
      for (jcol=0;jcol<3;jcol++)
	{
	  for (jrow=0;jrow<3;jrow++)
	    {
	      QDP_C_eq_elem_M(complex_temp,durdsigma,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_durdsigma,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_durdsigma,complex_temp,QDP_all);
	      
	      QDP_R_eq_R_times_R(real_temp1,re_durdsigma,z_norm,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,im_durdsigma,z_norm,QDP_all);
	      
	      QDP_C_eq_R_plus_i_R(complex_temp,real_temp1,real_temp2,QDP_all);
	      QDP_M_eq_elem_C(durdsigma,complex_temp,jcol,jrow,QDP_all);
	    }
	}


      for (jcol=0;jcol<3;jcol++)
	{
	  for (jrow=0;jrow<3;jrow++)
	    {
	      
	      QDP_C_eq_elem_M(complex_temp,mW,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_mW,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_mW,complex_temp,QDP_all);
	      
	      QDP_C_eq_elem_M(complex_temp,durdsigma,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_durdsigma,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_durdsigma,complex_temp,QDP_all);
	      
	      QDP_R_meq_R_times_R(re_durdsigma,re_mW,re_ctemp,QDP_all);
	      QDP_R_peq_R_times_R(re_durdsigma,im_mW,im_ctemp,QDP_all);
	      
	      QDP_R_meq_R_times_R(im_durdsigma,re_mW,im_ctemp,QDP_all);
	      QDP_R_meq_R_times_R(im_durdsigma,im_mW,re_ctemp,QDP_all);
	      
	      QDP_C_eq_R_plus_i_R(complex_temp,re_durdsigma,
				  im_durdsigma,QDP_all);
	      QDP_M_eq_elem_C(durdsigma,complex_temp,jcol,jrow,QDP_all);
	      
	    }
	}
      
      //!*------------
      //!* durdsigmad
      //!*------------
      //!* Only has the messy term:
      
      QDP_M_eq_zero(dWdsigmad,QDP_all);
      for (ipoly=1;ipoly<=npolyterm;ipoly++)
	{
	  QDP_R_eq_r(poly_temp,&poly_c[ipoly],QDP_all);
	  
	  for (jcol=0;jcol<3;jcol++)
	    {
	      for (jrow=0;jrow<3;jrow++)
		{

	  QDP_C_eq_elem_M(complex_temp,dWdsigmad,jrow,jcol,QDP_all);
	  QDP_R_eq_re_C(re_dWdsigmad,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_dWdsigmad,complex_temp,QDP_all);
	  
	  QDP_C_eq_elem_M(complex_temp,Sterm[ipoly],jrow,irow,QDP_all);
	  QDP_R_eq_re_C(re_Sterm1,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_Sterm1,complex_temp,QDP_all);

	  QDP_C_eq_elem_M(complex_temp,Sterm[ipoly],icol,jcol,QDP_all);
	  QDP_R_eq_re_C(re_Sterm2,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_Sterm2,complex_temp,QDP_all);

	  //real part
	  QDP_R_eq_R_times_R(real_temp,re_Sterm1,re_Sterm2,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_Sterm1,im_Sterm2,QDP_all);
	  QDP_R_eq_R_minus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_meq_R_times_R(re_dWdsigmad,poly_temp,real_temp2,QDP_all);

	  //im part
	  QDP_R_eq_R_times_R(real_temp,re_Sterm1,im_Sterm2,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_Sterm1,re_Sterm2,QDP_all);
	  QDP_R_eq_R_plus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_meq_R_times_R(im_dWdsigmad,poly_temp,real_temp2,QDP_all);

	  QDP_C_eq_R_plus_i_R(complex_temp,re_dWdsigmad,im_dWdsigmad,QDP_all);
	  QDP_M_eq_elem_C(dWdsigmad,complex_temp,jrow,jcol,QDP_all);  

		}
	      }
	}


      //!* Finally durdsigmad:

      //calculate ctemp
      QDP_M_eq_M_times_M(mtemp1,mWadj,dWdsigmad,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,0,0,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_eq_R(im_ctemp,im_mtemp1,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,1,1,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_peq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_peq_R(im_ctemp,im_mtemp1,QDP_all);

      QDP_C_eq_elem_M(complex_temp,mtemp1,2,2,QDP_all);
      QDP_R_eq_re_C(re_mtemp1,complex_temp,QDP_all);
      QDP_R_eq_im_C(im_mtemp1,complex_temp,QDP_all);
      QDP_R_peq_R(re_ctemp,re_mtemp1,QDP_all);
      QDP_R_peq_R(im_ctemp,im_mtemp1,QDP_all);

      //update ctemp real part
      QDP_R_eq_R_times_R(real_temp1,re_ctemp,re_det43,QDP_all);
      QDP_R_eq_R_times_R(real_temp2,im_ctemp,im_det43,QDP_all);
      QDP_R_eq_R_plus_R(re_ctemp,real_temp1,real_temp2,QDP_all);
      
      //update ctemp im part
      QDP_R_eq_R_times_R(real_temp1,im_ctemp,re_det43,QDP_all);
      QDP_R_eq_R_times_R(real_temp2,re_ctemp,im_det43,QDP_all);
      QDP_R_eq_R_minus_R(im_ctemp,real_temp1,real_temp2,QDP_all);
      
      //znorm
      QDP_R_eq_R_times_R(z_norm,re_det43,re_det43,QDP_all);
      QDP_R_peq_R_times_R(z_norm,im_det43,im_det43,QDP_all);

      //update ctemp
      qla_const=3.;
      QDP_R_eq_r_times_R(real_temp,&qla_const,z_norm,QDP_all);

      QDP_R_eq_R_divide_R(real_temp2,re_ctemp,real_temp,QDP_all);
      QDP_R_eq_R(re_ctemp,real_temp2,QDP_all);

      QDP_R_eq_R_divide_R(real_temp2,im_ctemp,real_temp,QDP_all);
      QDP_R_eq_R(im_ctemp,real_temp2,QDP_all);


	for (jcol=0;jcol<3;jcol++)
	  {
	  for (jrow=0;jrow<3;jrow++)
	    {
	      
	      QDP_C_eq_elem_M(complex_temp,dWdsigmad,jrow,jcol,QDP_all);
	      QDP_R_eq_re_C(re_dWdsigmad,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_dWdsigmad,complex_temp,QDP_all);
	      
	      
	      QDP_R_eq_R_times_R(real_temp1,re_dWdsigmad,re_det13,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,im_dWdsigmad,im_det13,QDP_all);
	      QDP_R_eq_R_plus_R(re_durdsigmad,real_temp1,real_temp2,QDP_all);
	      
	      QDP_R_eq_R_times_R(real_temp1,im_dWdsigmad,re_det13,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,re_dWdsigmad,im_det13,QDP_all);
	      QDP_R_eq_R_minus_R(im_durdsigmad,real_temp1,real_temp2,QDP_all);
	      QDP_C_eq_R_plus_i_R(complex_temp,re_durdsigmad,im_durdsigmad,
				  QDP_all);
	      QDP_M_eq_elem_C(durdsigmad,complex_temp,jrow,jcol,QDP_all);


	    }
	  }


      //znorm
      QDP_R_eq_R_times_R(real_temp,re_det13,re_det13,QDP_all);
      QDP_R_peq_R_times_R(real_temp,im_det13,im_det13,QDP_all);
      
      qla_const=1.;
      QDP_R_eq_r(real_temp1,&qla_const,QDP_all);
      QDP_R_eq_R_divide_R(z_norm,real_temp1,real_temp,QDP_all);
      
      //update durdsigmad
      for (jcol=0;jcol<3;jcol++)
	{
	  for (jrow=0;jrow<3;jrow++)
	    {
	      QDP_C_eq_elem_M(complex_temp,durdsigmad,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_durdsigmad,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_durdsigmad,complex_temp,QDP_all);
	      
	      QDP_R_eq_R_times_R(real_temp1,re_durdsigmad,z_norm,QDP_all);
	      QDP_R_eq_R_times_R(real_temp2,im_durdsigmad,z_norm,QDP_all);
	      
	      QDP_C_eq_R_plus_i_R(complex_temp,real_temp1,real_temp2,QDP_all);
	      QDP_M_eq_elem_C(durdsigmad,complex_temp,jcol,jrow,QDP_all);
	    }
	}
      

      
      for (jcol=0;jcol<3;jcol++)
	{
	  for (jrow=0;jrow<3;jrow++)
	    {
	      
	      QDP_C_eq_elem_M(complex_temp,mW,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_mW,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_mW,complex_temp,QDP_all);
	      
	      QDP_C_eq_elem_M(complex_temp,durdsigmad,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_durdsigmad,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_durdsigmad,complex_temp,QDP_all);
	      
	      QDP_R_meq_R_times_R(re_durdsigmad,re_mW,re_ctemp,QDP_all);
	      QDP_R_peq_R_times_R(re_durdsigmad,im_mW,im_ctemp,QDP_all);
	      
	      QDP_R_meq_R_times_R(im_durdsigmad,re_mW,im_ctemp,QDP_all);
	      QDP_R_meq_R_times_R(im_durdsigmad,im_mW,re_ctemp,QDP_all);

	      QDP_C_eq_R_plus_i_R(complex_temp,re_durdsigmad,
				  im_durdsigmad,QDP_all);
	      QDP_M_eq_elem_C(durdsigmad,complex_temp,jcol,jrow,QDP_all);

	    }
	}
      

      //!* Take T=Tr[fer*durdsigma] and Tp=Tr[fer*durdsigmad]:

      for (jrow=0;jrow<3;jrow++)
	{
	  for (jcol=0;jcol<3;jcol++)
	    {
	      
	      
	      
	      
	      QDP_C_eq_elem_M(complex_temp,Tlocal,icol,irow,QDP_all);
	      QDP_R_eq_re_C(re_Tlocal,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_Tlocal,complex_temp,QDP_all);
	      
	      
	      
	      QDP_C_eq_elem_M(complex_temp,f0[id],jrow,jcol,QDP_all);
	      QDP_R_eq_re_C(re_f0,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_f0,complex_temp,QDP_all);
	      
	      
	      
	      QDP_C_eq_elem_M(complex_temp,durdsigma,jcol,jrow,QDP_all);
	      QDP_R_eq_re_C(re_durdsigma,complex_temp,QDP_all);
	      QDP_R_eq_im_C(im_durdsigma,complex_temp,QDP_all);
	      
	      
	      
	      //real part
	      QDP_R_eq_R_times_R(real_temp,re_f0,re_durdsigma,QDP_all);
	      QDP_R_eq_R_times_R(real_temp1,im_f0,im_durdsigma,QDP_all);
	      QDP_R_eq_R_minus_R(real_temp2,real_temp,real_temp1,QDP_all);
	      QDP_R_peq_R(re_Tlocal,real_temp2,QDP_all);
	      
	      //im part
	      QDP_R_eq_R_times_R(real_temp,re_f0,im_durdsigma,QDP_all);
	      QDP_R_eq_R_times_R(real_temp1,im_f0,re_durdsigma,QDP_all);
	      QDP_R_eq_R_plus_R(real_temp2,real_temp,real_temp1,QDP_all);
	      QDP_R_peq_R(im_Tlocal,real_temp2,QDP_all);
	      
	      
	      QDP_C_eq_R_plus_i_R(complex_temp,re_Tlocal,im_Tlocal,QDP_all);
	      QDP_M_eq_elem_C(Tlocal,complex_temp,icol,irow,QDP_all);  
	  
	      
	    }
	}



	for (jrow=0;jrow<3;jrow++)
	  {
	  for (jcol=0;jcol<3;jcol++)
		{
	  
	  QDP_C_eq_elem_M(complex_temp,Tplocal,icol,irow,QDP_all);
	  QDP_R_eq_re_C(re_Tplocal,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_Tplocal,complex_temp,QDP_all);
	  
	  QDP_C_eq_elem_M(complex_temp,f0[id],jrow,jcol,QDP_all);
	  	  QDP_R_eq_re_C(re_f0,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_f0,complex_temp,QDP_all);

	  QDP_C_eq_elem_M(complex_temp,durdsigmad,jcol,jrow,QDP_all);
	  QDP_R_eq_re_C(re_durdsigmad,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_durdsigmad,complex_temp,QDP_all);

	  //real part
	  QDP_R_eq_R_times_R(real_temp,re_f0,re_durdsigmad,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_f0,im_durdsigmad,QDP_all);
	  QDP_R_eq_R_minus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_peq_R(re_Tplocal,real_temp2,QDP_all);
			     
	  //im part
	  QDP_R_eq_R_times_R(real_temp,re_f0,im_durdsigmad,QDP_all);
	  QDP_R_eq_R_times_R(real_temp1,im_f0,re_durdsigmad,QDP_all);
	  QDP_R_eq_R_plus_R(real_temp2,real_temp,real_temp1,QDP_all);
	  QDP_R_peq_R(im_Tplocal,real_temp2,QDP_all);

	  QDP_C_eq_R_plus_i_R(complex_temp,re_Tplocal,im_Tplocal,QDP_all);
	  QDP_M_eq_elem_C(Tplocal,complex_temp,icol,irow,QDP_all);  
	  

		}
	      }

	//!*=================================
	//!* End of loop over T's components
	//!*=================================
	
    }//end loop over irow
  }//end loop over icol


  //Finally compute T:

  for (jcol=0;jcol<3;jcol++)
    {
      for (jrow=0;jrow<3;jrow++)
	{
	  
	  QDP_C_eq_elem_M(complex_temp,Tlocal,jrow,jcol,QDP_all);
	  QDP_R_eq_re_C(re_Tlocal,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_Tlocal,complex_temp,QDP_all);
	  
	  QDP_C_eq_elem_M(complex_temp,Tplocal,jcol,jrow,QDP_all);
	  QDP_R_eq_re_C(re_Tplocal,complex_temp,QDP_all);
	  QDP_R_eq_im_C(im_Tplocal,complex_temp,QDP_all);
	  
	  QDP_R_eq_R_plus_R(real_temp,re_Tlocal,re_Tplocal,QDP_all);
	  QDP_R_eq_R_minus_R(real_temp2,im_Tlocal,im_Tplocal,QDP_all);
	  
	  QDP_C_eq_R_plus_i_R(complex_temp,real_temp,real_temp2,
			      QDP_all);
	  QDP_M_eq_elem_C(f1[id],complex_temp,jrow,jcol,QDP_all);  
	  
	}
      
    }
  
  
    }//end loop over id
 
 
 //Destroy fields
 QDP_destroy_M(msquare);
 QDP_destroy_M(mWadj);
 
 QDP_destroy_M (mtemp1);
 QDP_destroy_R(re_mtemp1); QDP_destroy_R(im_mtemp1);
 
 QDP_destroy_M(dWdsigma);
 QDP_destroy_R(re_dWdsigma); QDP_destroy_R(im_dWdsigma);
 
 QDP_destroy_M(dWdsigmad);
 QDP_destroy_R(re_dWdsigmad); QDP_destroy_R(im_dWdsigmad);
 
 QDP_destroy_M(mW);
 QDP_destroy_R(re_mW); QDP_destroy_R(im_mW);
 
 QDP_destroy_M(durdsigma);
 QDP_destroy_R(re_durdsigma); QDP_destroy_R(im_durdsigma);
 
 QDP_destroy_M(durdsigmad);
 QDP_destroy_R(re_durdsigmad);QDP_destroy_R(im_durdsigmad);
 
 QDP_destroy_M(Tlocal);
 QDP_destroy_R(re_Tlocal); QDP_destroy_R(im_Tlocal);
 
 QDP_destroy_M(Tplocal);
 QDP_destroy_R(re_Tplocal); QDP_destroy_R(im_Tplocal);
 
 QDP_destroy_M(sumterm);
 QDP_destroy_R(re_sumterm); QDP_destroy_R(im_sumterm);
 
 QDP_destroy_R(re_f0); QDP_destroy_R(im_f0);
 
 for(j=0;j<=npolyterm;j++) QDP_destroy_M(Sterm[j]);
 QDP_destroy_R(re_Sterm1); QDP_destroy_R(im_Sterm1);
 QDP_destroy_R(re_Sterm2); QDP_destroy_R(im_Sterm2);
 
 for(j=0;j<=npolyterm;j++) QDP_destroy_M(StermSd[j]);
 QDP_destroy_R(re_StermSd); QDP_destroy_R(im_StermSd);
 
 for(j=0;j<=npolyterm;j++) QDP_destroy_M(term[j]);
 QDP_destroy_R(re_term); QDP_destroy_R(im_term);
 
 QDP_destroy_I(logical_temp);
 QDP_destroy_C(complex_temp);
 
 QDP_destroy_C(det);
 QDP_destroy_R(re_det);
 QDP_destroy_R(im_det);
 QDP_destroy_R(re_det13);
 QDP_destroy_R(im_det13);
 QDP_destroy_R(re_det43);
 QDP_destroy_R(im_det43);
 QDP_destroy_R(real_temp);
 QDP_destroy_R(real_temp1);
 QDP_destroy_R(real_temp2);
 QDP_destroy_R(poly_temp);
 QDP_destroy_R(z_norm);
 QDP_destroy_R(z_arg);
 QDP_destroy_R(qdp_const);
 QDP_destroy_R(re_ctemp);
 QDP_destroy_R(im_ctemp);

 return;


}
#endif /* ALAN_GRAY_REUNIT_FORCE */



//utilities
//find inverse of SU3 matrix
static void QOPPC(su3inverse)(QDP_ColorMatrix *mat,QDP_ColorMatrix *inv)
{

  int i,j,i1=0,j1=0,i2=0,j2=0;

  QDP_Real *re_cofactor[3][3];
  QDP_Real *im_cofactor[3][3];
  QDP_Real *re_mat[3][3];
  QDP_Real *im_mat[3][3];
  QDP_Real *re_inv;
  QDP_Real *im_inv;
  QDP_Real *re_det;
  QDP_Real *im_det;
  QDP_Complex *complextemp;
  QDP_Real *re_temp;
  QDP_Real *im_temp;
  QDP_Real *QDPone;
  QDP_Real *absdet;
  QLA_Real minusone=-1.;
  QLA_Real one=1.;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      re_mat[i][j]=QDP_create_R();
      im_mat[i][j]=QDP_create_R();
      re_cofactor[i][j]=QDP_create_R();
      im_cofactor[i][j]=QDP_create_R();
    }
  }
  re_det=QDP_create_R();
  im_det=QDP_create_R();
  re_inv=QDP_create_R();
  im_inv=QDP_create_R();
  complextemp=QDP_create_C();
  re_temp=QDP_create_R();
  im_temp=QDP_create_R();
  absdet=QDP_create_R();
  QDPone=QDP_create_R();

  //extract matrix elements
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      QDP_C_eq_elem_M(complextemp,mat,i,j,QDP_all);
      QDP_R_eq_re_C(re_mat[i][j],complextemp,QDP_all);
      QDP_R_eq_im_C(im_mat[i][j],complextemp,QDP_all);
    }
  }

  // calculate cofactors

  for(i=0;i<3;i++){

      switch (i)
	{
	case 0:
	  i1=1;i2=2;
	  break;
	case 1:
	  i1=0;i2=2;
	  break;
	case 2:
	  i1=0;i2=1;
	  break;
	}

    for(j=0;j<3;j++){

      switch (j)
	{
	case 0:
	  j1=1;j2=2;
	  break;
	case 1:
	  j1=0;j2=2;
	  break;
	case 2:
	  j1=0;j2=1;
	  break;
	}

      //cofactor real part
      QDP_R_eq_R_times_R(re_cofactor[i][j],re_mat[i1][j1],re_mat[i2][j2],
			 QDP_all);
      QDP_R_meq_R_times_R(re_cofactor[i][j],im_mat[i1][j1],im_mat[i2][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(re_cofactor[i][j],re_mat[i2][j1],re_mat[i1][j2],
			  QDP_all);
      QDP_R_peq_R_times_R(re_cofactor[i][j],im_mat[i2][j1],im_mat[i1][j2],
			  QDP_all);
      
      //cofactor im part
      QDP_R_eq_R_times_R(im_cofactor[i][j],re_mat[i1][j1],im_mat[i2][j2],
			 QDP_all);
      QDP_R_peq_R_times_R(im_cofactor[i][j],im_mat[i1][j1],re_mat[i2][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(im_cofactor[i][j],re_mat[i2][j1],im_mat[i1][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(im_cofactor[i][j],im_mat[i2][j1],re_mat[i1][j2],
			  QDP_all);
    }
  }

  //negate as necessary
  QDP_R_eq_R(re_temp,re_cofactor[0][1],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[0][1],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[0][1],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[0][1],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[1][0],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[1][0],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[1][0],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[1][0],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[1][2],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[1][2],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[1][2],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[1][2],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[2][1],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[2][1],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[2][1],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[2][1],&minusone,re_temp,QDP_all);


  //calculate det
  QDP_R_eq_zero(re_det,QDP_all);
  QDP_R_eq_zero(im_det,QDP_all);

  for(i=0;i<3;i++){
 
    QDP_R_peq_R_times_R(re_det,re_mat[i][0],re_cofactor[i][0],QDP_all);
    QDP_R_meq_R_times_R(re_det,im_mat[i][0],im_cofactor[i][0],QDP_all);
    QDP_R_peq_R_times_R(im_det,re_mat[i][0],im_cofactor[i][0],QDP_all);
    QDP_R_peq_R_times_R(im_det,im_mat[i][0],re_cofactor[i][0],QDP_all);
  }


  //calculate inverse

  QDP_R_eq_R_times_R(re_temp,re_det,re_det,QDP_all);
  QDP_R_peq_R_times_R(re_temp,im_det,im_det,QDP_all);
  QDP_R_eq_r(QDPone,&one,QDP_all);
  QDP_R_eq_R_divide_R(absdet,QDPone,re_temp,QDP_all);


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
 
    QDP_R_eq_R_times_R(re_inv,re_cofactor[i][j],re_det,QDP_all);
    QDP_R_peq_R_times_R(re_inv,im_cofactor[i][j],im_det,QDP_all);

    QDP_R_eq_R_times_R(im_inv,im_cofactor[i][j],re_det,QDP_all);
    QDP_R_meq_R_times_R(im_inv,re_cofactor[i][j],im_det,QDP_all);


    QDP_R_eq_R_times_R(re_temp,absdet,re_inv,QDP_all);
    QDP_R_eq_R_times_R(im_temp,absdet,im_inv,QDP_all);

    //insert real,im parts back to complex
    QDP_C_eq_R_plus_i_R(complextemp,re_temp,im_temp,QDP_all);
    //insert complex to color matrix
    QDP_M_eq_elem_C(inv,complextemp,j,i,QDP_all);

    }
  }


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      QDP_destroy_R(re_mat[i][j]);
      QDP_destroy_R(im_mat[i][j]);
      QDP_destroy_R(re_cofactor[i][j]);
      QDP_destroy_R(im_cofactor[i][j]);
    }
  }
  QDP_destroy_R(re_det);
  QDP_destroy_R(im_det);
  QDP_destroy_R(re_inv);
  QDP_destroy_R(im_inv);
  QDP_destroy_C(complextemp);
  QDP_destroy_R(re_temp);
  QDP_destroy_R(im_temp);
  QDP_destroy_R(absdet);
  QDP_destroy_R(QDPone);

  return;

}


//find determinant of SU3 matrix
static void QOPPC(su3det)(QDP_ColorMatrix *mat,QDP_Complex *det)
{

  int i,j;
  int j1=0,j2=0;

  QDP_Real *re_cofactor[3];
  QDP_Real *im_cofactor[3];
  QDP_Real *re_mat[3][3];
  QDP_Real *im_mat[3][3];
  QDP_Real *re_det;
  QDP_Real *im_det;
  QDP_Complex *complextemp;
  QDP_Real *rtemp;
  QLA_Real minusone=-1.;

  for(i=0;i<3;i++){
      re_cofactor[i]=QDP_create_R();
      im_cofactor[i]=QDP_create_R();
    for(j=0;j<3;j++){
      re_mat[i][j]=QDP_create_R();
      im_mat[i][j]=QDP_create_R();
    }
  }
  re_det=QDP_create_R();
  im_det=QDP_create_R();
  complextemp=QDP_create_C();
  rtemp=QDP_create_R();

  //extract 9 retemp and 9 imtemp QDP reals
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      QDP_C_eq_elem_M(complextemp,mat,i,j,QDP_all);
      QDP_R_eq_re_C(re_mat[i][j],complextemp,QDP_all);
      QDP_R_eq_im_C(im_mat[i][j],complextemp,QDP_all);
    }
  }

  // calculate cofactors
  // Expand about the first row:

  for(i=0;i<3;i++){
    if (i==0){
      j1=1; j2=2; 
    }
    if (i==1){
      j1=0; j2=2; 
    }
    if (i==2){
      j1=0; j2=1; 
    }

    //cofactor real part
    QDP_R_eq_R_times_R(re_cofactor[i],re_mat[1][j1],re_mat[2][j2],QDP_all);
    QDP_R_meq_R_times_R(re_cofactor[i],im_mat[1][j1],im_mat[2][j2],QDP_all);
    QDP_R_meq_R_times_R(re_cofactor[i],re_mat[1][j2],re_mat[2][j1],QDP_all);
    QDP_R_peq_R_times_R(re_cofactor[i],im_mat[1][j2],im_mat[2][j1],QDP_all);
    
    //cofactor im part
    QDP_R_eq_R_times_R(im_cofactor[i],re_mat[1][j1],im_mat[2][j2],QDP_all);
    QDP_R_peq_R_times_R(im_cofactor[i],im_mat[1][j1],re_mat[2][j2],QDP_all);
    QDP_R_meq_R_times_R(im_cofactor[i],re_mat[1][j2],im_mat[2][j1],QDP_all);
    QDP_R_meq_R_times_R(im_cofactor[i],im_mat[1][j2],re_mat[2][j1],QDP_all);
  }
  
  //negate re_cofactor[1], im_cofactor[1]
  QDP_R_eq_R(rtemp,re_cofactor[1],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[1],&minusone,rtemp,QDP_all);

  QDP_R_eq_R(rtemp,im_cofactor[1],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[1],&minusone,rtemp,QDP_all);
  
  
  //calculate det
  QDP_R_eq_zero(re_det,QDP_all);
  QDP_R_eq_zero(im_det,QDP_all);

  for(i=0;i<3;i++){

 
    QDP_R_peq_R_times_R(re_det,re_mat[0][i],re_cofactor[i],QDP_all);
    QDP_R_meq_R_times_R(re_det,im_mat[0][i],im_cofactor[i],QDP_all);
    QDP_R_peq_R_times_R(im_det,re_mat[0][i],im_cofactor[i],QDP_all);
    QDP_R_peq_R_times_R(im_det,im_mat[0][i],re_cofactor[i],QDP_all);
  }

  QDP_C_eq_R_plus_i_R(det,re_det,im_det,QDP_all);

 

    for(i=0;i<3;i++){
      QDP_destroy_R(re_cofactor[i]);
      QDP_destroy_R(im_cofactor[i]);
    for(j=0;j<3;j++){
      QDP_destroy_R(re_mat[i][j]);
      QDP_destroy_R(im_mat[i][j]);
    }
  }
  QDP_destroy_R(re_det);
  QDP_destroy_R(im_det);
  QDP_destroy_C(complextemp);
  QDP_destroy_R(rtemp);
  
  return;


}


#if 0
//find adjoint of SU3 matrix
static void QOPPC(su3adjoint)(QDP_ColorMatrix *mat,QDP_ColorMatrix *adj)
{

  int i,j,i1=0,j1=0,i2=0,j2=0;

  QDP_Real *re_cofactor[3][3];
  QDP_Real *im_cofactor[3][3];
  QDP_Real *re_mat[3][3];
  QDP_Real *im_mat[3][3];
  QDP_Real *re_inv;
  QDP_Real *im_inv;
  QDP_Real *re_det;
  QDP_Real *im_det;
  QDP_Complex *complextemp;
  QDP_Real *re_temp;
  QDP_Real *im_temp;
  QDP_Real *QDPone;
  QDP_Real *absdet;
  QLA_Real minusone=-1.;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      re_mat[i][j]=QDP_create_R();
      im_mat[i][j]=QDP_create_R();
      re_cofactor[i][j]=QDP_create_R();
      im_cofactor[i][j]=QDP_create_R();
    }
  }
  re_det=QDP_create_R();
  im_det=QDP_create_R();
  re_inv=QDP_create_R();
  im_inv=QDP_create_R();
  complextemp=QDP_create_C();
  re_temp=QDP_create_R();
  im_temp=QDP_create_R();
  absdet=QDP_create_R();
  QDPone=QDP_create_R();

  //extract matrix elements
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      QDP_C_eq_elem_M(complextemp,mat,i,j,QDP_all);
      QDP_R_eq_re_C(re_mat[i][j],complextemp,QDP_all);
      QDP_R_eq_im_C(im_mat[i][j],complextemp,QDP_all);
    }
  }

  // calculate cofactors

  for(i=0;i<3;i++){

      switch (i)
	{
	case 0:
	  i1=1;i2=2;
	  break;
	case 1:
	  i1=0;i2=2;
	  break;
	case 2:
	  i1=0;i2=1;
	  break;
	}

    for(j=0;j<3;j++){

      switch (j)
	{
	case 0:
	  j1=1;j2=2;
	  break;
	case 1:
	  j1=0;j2=2;
	  break;
	case 2:
	  j1=0;j2=1;
	  break;
	}

      //cofactor real part
      QDP_R_eq_R_times_R(re_cofactor[i][j],re_mat[i1][j1],re_mat[i2][j2],
			 QDP_all);
      QDP_R_meq_R_times_R(re_cofactor[i][j],im_mat[i1][j1],im_mat[i2][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(re_cofactor[i][j],re_mat[i2][j1],re_mat[i1][j2],
			  QDP_all);
      QDP_R_peq_R_times_R(re_cofactor[i][j],im_mat[i2][j1],im_mat[i1][j2],
			  QDP_all);
      
      //cofactor im part
      QDP_R_eq_R_times_R(im_cofactor[i][j],re_mat[i1][j1],im_mat[i2][j2],
			 QDP_all);
      QDP_R_peq_R_times_R(im_cofactor[i][j],im_mat[i1][j1],re_mat[i2][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(im_cofactor[i][j],re_mat[i2][j1],im_mat[i1][j2],
			  QDP_all);
      QDP_R_meq_R_times_R(im_cofactor[i][j],im_mat[i2][j1],re_mat[i1][j2],
			  QDP_all);
    }
  }

  //negate as necessary
  QDP_R_eq_R(re_temp,re_cofactor[0][1],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[0][1],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[0][1],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[0][1],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[1][0],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[1][0],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[1][0],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[1][0],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[1][2],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[1][2],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[1][2],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[1][2],&minusone,re_temp,QDP_all);

  QDP_R_eq_R(re_temp,re_cofactor[2][1],QDP_all);
  QDP_R_eq_r_times_R(re_cofactor[2][1],&minusone,re_temp,QDP_all);
  
  QDP_R_eq_R(re_temp,im_cofactor[2][1],QDP_all);
  QDP_R_eq_r_times_R(im_cofactor[2][1],&minusone,re_temp,QDP_all);


  


  //calculate adjoint


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
 
 
    //insert real,im parts back to complex

    QDP_C_eq_R_plus_i_R(complextemp,re_cofactor[i][j],im_cofactor[i][j],QDP_all);

   

    //insert complex to color matrix
    QDP_M_eq_elem_C(adj,complextemp,j,i,QDP_all);

    }
  }


  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      QDP_destroy_R(re_mat[i][j]);
      QDP_destroy_R(im_mat[i][j]);
      QDP_destroy_R(re_cofactor[i][j]);
      QDP_destroy_R(im_cofactor[i][j]);
    }
  }
  QDP_destroy_R(re_det);
  QDP_destroy_R(im_det);
  QDP_destroy_R(re_inv);
  QDP_destroy_R(im_inv);
  QDP_destroy_C(complextemp);
  QDP_destroy_R(re_temp);
  QDP_destroy_R(im_temp);
  QDP_destroy_R(absdet);
  QDP_destroy_R(QDPone);

  return;

}

#endif
