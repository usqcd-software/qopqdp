
#include <qop_internal.h>

static void QOPPC(su3det)(QDP_ColorMatrix *mat,QDP_Complex *det);
static void QOPPC(su3inverse)(QDP_ColorMatrix *mat,QDP_ColorMatrix *inv);
static void QOPPC(su3reunit)(QDP_ColorMatrix *U,QDP_ColorMatrix *Ur);


//create HISQ links from QOPPC(GaugeField)
QOPPC(FermionLinksHisq) *
QOPPC(hisq_create_L_from_G)(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOPPC(GaugeField) *gauge)
{
  QOPPC(FermionLinksAsqtad) *fla;
  QOPPC(FermionLinksHisq) *flh;
  QOP_asqtad_coeffs_t acoeffs;
  int i;
  QDP_ColorMatrix *tempM;


  tempM=QDP_create_M();

  // fat7 stage 
  // convert hisq style coeffs to asqtad style 
  acoeffs.one_link = coeffs->fat7_one_link;
  acoeffs.three_staple = coeffs->fat7_three_staple;
  acoeffs.five_staple = coeffs->fat7_five_staple;
  acoeffs.seven_staple = coeffs->fat7_seven_staple;
  // by definition lepage and naik coeffs are 0 for fat7
  acoeffs.lepage = 0.;
  acoeffs.naik = 0.;


  // smear gauge with fat7 
  fla=QOP_asqtad_create_L_from_G(info, &acoeffs, gauge);


  // convert fla back to gauge to enable further smearing 
  // since fat7 contains only fat links, can just copy fla->fatlinks to gauge
  for (i=0;i<4;i++)
    {
      QDP_M_eq_M(gauge->links[i], fla->fatlinks[i], QDP_all);
    }


  // projection stage 
  // project gauge fields to SU3 
  for (i=0;i<4;i++)
    {
      QOPPC(su3reunit)(gauge->links[i],tempM);
      QDP_M_eq_M(gauge->links[i], tempM, QDP_all);
    }


  // asqtad stage 
  // convert hisq style coeffs to asqtad style 
  acoeffs.one_link = coeffs->asqtad_one_link;
  acoeffs.three_staple = coeffs->asqtad_three_staple;
  acoeffs.five_staple = coeffs->asqtad_five_staple;
  acoeffs.seven_staple = coeffs->asqtad_seven_staple;
  acoeffs.lepage = coeffs->asqtad_lepage;
  acoeffs.naik = coeffs->asqtad_naik;

  // smear gauge with asqtad 
  fla=QOP_asqtad_create_L_from_G(info, &acoeffs, gauge);

  // copy fla pointer to flh pointer
  flh= (QOPPC(FermionLinksHisq)*) fla;

  QDP_destroy_M(tempM);

  return flh;

}

//wrapper of asqtad inverter for hisq links
void
QOPPC(hisq_invert)(QOP_info_t *info,
		     QOPPC(FermionLinksHisq) *flh,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     REAL mass,
		     QOPPC(ColorVector) *out,
		     QOPPC(ColorVector) *in)
{

  QOPPC(FermionLinksAsqtad) *fla;
  fla=(QOPPC(FermionLinksAsqtad)*) flh;
  
  QOP_asqtad_invert(info, fla, inv_arg, res_arg, mass, out, in);


}

//wrapper of asqtad inverter for hisq links
void
QOPPC(hisq_invert_multi)(QOP_info_t *info,
			   QOP_FermionLinksHisq *flh,
			   QOP_invert_arg_t *inv_arg,
			   QOP_resid_arg_t **res_arg[],
			   REAL *masses[], int nmass[],
			   QOP_ColorVector **out_pt[],
			   QOP_ColorVector *in_pt[],
			   int nsrc)
{

  QOPPC(FermionLinksAsqtad) *fla;
  fla=(QOPPC(FermionLinksAsqtad)*) flh;
  
  QOP_asqtad_invert_multi(info, fla, inv_arg, res_arg, masses,
			  nmass, out_pt, in_pt, nsrc);

}


//The following reunit, inverse and det routines were translated from 
//Kit Wong's serial fortran HISQ code.

//Reuniterise links to SU3
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

//find inverse of SU3 matrix
void QOPPC(su3inverse)(QDP_ColorMatrix *mat,QDP_ColorMatrix *inv)
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
void QOPPC(su3det)(QDP_ColorMatrix *mat,QDP_Complex *det)
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
