#define QOP_Precision 'D'
#include <qop_internal.h>

/* Copy with multiplication from double to single precision */

static QOP_F3_FermionLinksAsqtad *
QOP_FD3_asqtad_create_L_from_r_times_L(QLA_F_Real *s,
				       QOP_D3_FermionLinksAsqtad *fla_src)
{
  QOP_F3_FermionLinksAsqtad *fla;
  QDP_D3_ColorMatrix **fl_src = fla_src->fatlinks;
  QDP_D3_ColorMatrix **ll_src = fla_src->longlinks;
  QDP_F3_ColorMatrix *fl_dst[4], *ll_dst[4];
  int i;

  ASQTAD_DSLASH_BEGIN;

  /* Copy and multiply. */

  for(i=0; i<4; i++) {
    QLA_F_Real two_s = 2.0 * (*s); /* Correct for the 1/2 in source */
    fl_dst[i] = QDP_F3_create_M();
    QDP_FD3_M_eq_M(fl_dst[i], fl_src[i], QDP_all);
    QDP_F3_M_eq_r_times_M(fl_dst[i], &two_s, fl_dst[i], QDP_all);
    if(ll_src) {
      ll_dst[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(ll_dst[i], ll_src[i], QDP_all);
      QDP_F3_M_eq_r_times_M(ll_dst[i], &two_s, ll_dst[i], QDP_all);
    }
  }

  /* Hand over the QDP fields to the links structure */

  if(ll_src)
    fla = QOP_F3_asqtad_convert_L_from_qdp(fl_dst, ll_dst);
  else
    fla = QOP_F3_asqtad_convert_L_from_qdp(fl_dst, NULL);


  ASQTAD_DSLASH_END;
  return fla;
}

/* Copy from double to single precision */

QOP_F3_FermionLinksAsqtad *
QOP_FD3_asqtad_create_L_from_L(QOP_D3_FermionLinksAsqtad *fla_src)
{
  QOP_F3_FermionLinksAsqtad *fla;
  QLA_F_Real one = 1.0;

  ASQTAD_DSLASH_BEGIN;

  /* Copy the source links */
  fla = QOP_FD3_asqtad_create_L_from_r_times_L(&one, fla_src);

  ASQTAD_DSLASH_END;
  return fla;
}

