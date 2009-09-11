/******* hisq_invert_p.c ****************/


#include <qop_internal.h>

extern void QOPPC(su3det)(QDP_ColorMatrix *mat,QDP_Complex *det);
extern void QOPPC(su3inverse)(QDP_ColorMatrix *mat,QDP_ColorMatrix *inv);
extern void QOPPC(su3reunit)(QDP_ColorMatrix *U,QDP_ColorMatrix *Ur);


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
  QDP_ColorMatrix *tempM, *Vgauge[4], *Wgauge[4];

  QOPPC(GaugeField) *qopgf_tmp;

  tempM=QDP_create_M();

  for (i=0;i<4;i++){
    Vgauge[i]=QDP_create_M();
    Wgauge[i]=QDP_create_M();
  }

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
      QLA_Real two = 2.0; // correct for factor of 1/2 in links
      //QDP_M_eq_M(gauge->links[i], fla->fatlinks[i], QDP_all);

      //EPCC correction: dont overwrite original gauge field
      //QDP_M_eq_r_times_M(gauge->links[i], &two, fla->fatlinks[i], QDP_all);

      QDP_M_eq_r_times_M(Vgauge[i], &two, fla->fatlinks[i], QDP_all);

    }
  QOP_asqtad_destroy_L(fla);

  // projection stage 
  // project gauge fields to SU3 
  for (i=0;i<4;i++)
    {
      //EPCC correction: dont overwrite original gauge field
      //QOPPC(su3reunit)(gauge->links[i],tempM);
      //QDP_M_eq_M(gauge->links[i], tempM, QDP_all);

      QOPPC(su3reunit)(Vgauge[i],Wgauge[i]);

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

  //EPCC correction: dont overwrite original gauge field
  //fla=QOP_asqtad_create_L_from_G(info, &acoeffs, gauge);

  qopgf_tmp = QOP_convert_G_from_qdp(Wgauge);
  fla=QOP_asqtad_create_L_from_G(info, &acoeffs, qopgf_tmp);

  // copy fla pointer to flh pointer
  flh= (QOPPC(FermionLinksHisq)*) fla;

  QDP_destroy_M(tempM);

  for (i=0;i<4;i++){
    QDP_destroy_M(Vgauge[i]);
    QDP_destroy_M(Wgauge[i]);
}

  QOP_destroy_G(qopgf_tmp);

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

