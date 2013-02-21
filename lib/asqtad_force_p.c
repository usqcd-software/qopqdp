#include <qop_internal.h>

void
QOP_asqtad_deriv_multi_qdp(QOP_info_t *info, QDP_ColorMatrix *links[],
			   QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
			   REAL eps[], QDP_ColorVector *in_pt[], int nsrc)
{
  ASQTAD_FORCE_BEGIN;

  QOP_asqtad_deriv_multi_fnmat_qdp(info, links, force, coef, eps, in_pt, nsrc);

  ASQTAD_FORCE_END;
}

void
QOP_asqtad_force_multi_qdp(QOP_info_t *info, QDP_ColorMatrix *links[],
			   QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
			   REAL eps[], QDP_ColorVector *in_pt[], int nsrc)
{
  ASQTAD_FORCE_BEGIN;

  if( nsrc < QOP_asqtad_ff.fnmat_src_min ) {
    int vl = QOP_asqtad_ff.veclength;
    QOP_info_t tinfo;
    info->final_sec = 0;
    info->final_flop = 0;
    for(int j=0; j<nsrc; j+=vl) {
      int ns = nsrc-j;
      if(ns>vl) ns = vl;
      QOP_asqtad_force_multi_asvec_qdp(&tinfo, links, force, coef, eps+j, in_pt+j, ns);
      info->final_sec += tinfo.final_sec;
      info->final_flop += tinfo.final_flop;
    }
  } else {
    QOP_asqtad_force_multi_fnmat_qdp(info, links, force, coef, eps, in_pt, nsrc);
  }

  ASQTAD_FORCE_END;
}

void
QOP_asqtad_force_multi(QOP_info_t *info, QOP_GaugeField *gauge,
		       QOP_Force *force, QOP_asqtad_coeffs_t *coef,
		       QLA_Real eps[], QOP_ColorVector *in_pt[], int nsrc)
{
  ASQTAD_FORCE_BEGIN;

  QDP_ColorVector *x[nsrc];
  for(int i=0; i<nsrc; i++) x[i] = in_pt[i]->cv;
  QOP_asqtad_force_multi_qdp(info, gauge->links, force->force, coef, eps, x, nsrc);

  ASQTAD_FORCE_END;
}

void
QOP_asqtad_force(QOP_info_t *info, QOP_GaugeField *gauge, QOP_Force *force,
		 QOP_asqtad_coeffs_t *coeffs, QLA_Real eps, QOP_ColorVector *in_pt)
{
  ASQTAD_FORCE_BEGIN;

  QOP_asqtad_force_multi(info, gauge, force, coeffs, &eps, &in_pt, 1);

  ASQTAD_FORCE_END;
}
