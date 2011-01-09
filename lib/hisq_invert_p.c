/******* hisq_invert_p.c ****************/
//#define DO_TRACE

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
  QDP_ColorMatrix *Vgauge[4], *Wgauge[4];
  QOPPC(GaugeField) *qopgf_tmp;

  TRACE;
  for(int i=0; i<4; i++) {
    Vgauge[i] = QDP_create_M();
    Wgauge[i] = QDP_create_M();
  }

  TRACE;
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
  fla = QOP_asqtad_create_L_from_G(info, &acoeffs, gauge);

  TRACE;
  // convert fla back to gauge to enable further smearing 
  // since fat7 contains only fat links, can just copy fla->fatlinks to gauge
  for(int i=0; i<4; i++) {
    QLA_Real two = 2.0; // correct for factor of 1/2 in links
    QDP_M_eq_r_times_M(Vgauge[i], &two, fla->fatlinks[i], QDP_all);
  }
  TRACE;
  QOP_asqtad_destroy_L(fla);

  TRACE;
  // projection stage 
  // project gauge fields to SU3 
  for(int i=0; i<4; i++) {
    QOPPC(su3reunit)(Vgauge[i], Wgauge[i]);
  }

  TRACE;
  // asqtad stage 
  // convert hisq style coeffs to asqtad style 
  acoeffs.one_link = coeffs->asqtad_one_link;
  acoeffs.three_staple = coeffs->asqtad_three_staple;
  acoeffs.five_staple = coeffs->asqtad_five_staple;
  acoeffs.seven_staple = coeffs->asqtad_seven_staple;
  acoeffs.lepage = coeffs->asqtad_lepage;
  acoeffs.naik = coeffs->asqtad_naik;

  // smear gauge with asqtad 
  qopgf_tmp = QOP_convert_G_from_qdp(Wgauge);
  fla = QOP_asqtad_create_L_from_G(info, &acoeffs, qopgf_tmp);

  TRACE;
  // copy fla pointer to flh pointer
  flh = (QOPPC(FermionLinksHisq)*) fla;

  for(int i=0; i<4; i++){
    QDP_destroy_M(Vgauge[i]);
    //QDP_destroy_M(Wgauge[i]);  don't destroy converted links
  }
  QOP_destroy_G(qopgf_tmp);

  TRACE;
  return flh;
}

// wrappers of asqtad dslash and inverter for hisq links

void
QOPPC(hisq_dslash_qdp)(QOP_info_t *info,
		       QOPPC(FermionLinksHisq) *flh,
		       REAL mass,
		       QDPPC(ColorVector) *out,
		       QDPPC(ColorVector) *in,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in)
{
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOPPC(asqtad_dslash_qdp)(info, fla, mass, out, in, eo_out, eo_in);
}

void
QOPPC(hisq_diaginv_qdp)(QOP_info_t *info,
			QOPPC(FermionLinksHisq) *flh,
			REAL mass,
			QDPPC(ColorVector) *out,
			QDPPC(ColorVector) *in,
			QOP_evenodd_t eo)
{
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOPPC(asqtad_diaginv_qdp)(info, fla, mass, out, in, eo);
}

void
QOPPC(hisq_invert)(QOP_info_t *info,
		   QOPPC(FermionLinksHisq) *flh,
		   QOP_invert_arg_t *inv_arg,
		   QOP_resid_arg_t *res_arg,
		   REAL mass,
		   QOPPC(ColorVector) *out,
		   QOPPC(ColorVector) *in)
{
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOP_asqtad_invert(info, fla, inv_arg, res_arg, mass, out, in);
}

void
QOPPC(hisq_invert_qdp)(QOP_info_t *info,
		       QOPPC(FermionLinksHisq) *flh,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       REAL mass,
		       QDPPC(ColorVector) *out,
		       QDPPC(ColorVector) *in)
{
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOP_asqtad_invert_qdp(info, fla, inv_arg, res_arg, mass, out, in);
}

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
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOP_asqtad_invert_multi(info, fla, inv_arg, res_arg, masses,
			  nmass, out_pt, in_pt, nsrc);
}

void
QOPPC(hisq_invert_multi_qdp)(QOP_info_t *info,
			     QOP_FermionLinksHisq *flh,
			     QOP_invert_arg_t *inv_arg,
			     QOP_resid_arg_t **res_arg[],
			     REAL *masses[], int nmass[],
			     QDP_ColorVector **out_pt[],
			     QDP_ColorVector *in_pt[],
			     int nsrc)
{
  QOPPC(FermionLinksAsqtad) *fla = (QOPPC(FermionLinksAsqtad)*) flh;
  QOP_asqtad_invert_multi_qdp(info, fla, inv_arg, res_arg, masses,
			      nmass, out_pt, in_pt, nsrc);
}
