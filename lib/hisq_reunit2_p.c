#include <qop_internal.h>
#if QOP_Colors == 1
#include <qla_f1.h>
#include <qla_d1.h>
#include <qla_df1.h>
#elif QOP_Colors == 2
#include <qla_f2.h>
#include <qla_d2.h>
#include <qla_df2.h>
#elif QOP_Colors == 3
#include <qla_f3.h>
#include <qla_d3.h>
#include <qla_df3.h>
#else
#include <qla_fn.h>
#include <qla_dn.h>
#include <qla_dfn.h>
#endif

static void
projectU_site_d(NCPROT1 QLA_D_ColorMatrix(*Ur), QLA_D_ColorMatrix(*U))
{
  QLA_D_ColorMatrix(M1);
  QLA_D_ColorMatrix(M2);
  QLA_D_M_eq_Ma_times_M(&M1, U, U);
  //QLA_D_M_eq_invsqrt_M(&M2, &M1);
  QLA_D_M_eq_invsqrtPH_M(&M2, &M1);
  QLA_D_M_eq_M_times_M(Ur, U, &M2);
}

static void
projectU_site(NCPROT1 QLA_ColorMatrix(*Ur), int i, void *args)
{
  QLA_ColorMatrix(*U) = &((QLA_ColorMatrix(*))args)[i];
#if QLA_Precision == 'F'
  QLA_D_ColorMatrix(Ud);
  QLA_D_ColorMatrix(Urd);
  QLA_DF_M_eq_M(&Ud, U);
  projectU_site_d(NCARG1 &Urd, &Ud);
  QLA_FD_M_eq_M(Ur, &Urd);
#else
  projectU_site_d(NCARG1 Ur, U);
#endif
}
#undef NC

// requires QLA >= 1.7.0 and QDP >= 1.9.0
void
QOP_projectU_qdp(QOP_info_t *info, QDP_ColorMatrix *pU, QDP_ColorMatrix *U)
{
#define NC QDP_get_nc(U)
  double dtime = -QOP_time();
  QLA_ColorMatrix(*Uq) = QDP_expose_M(U);
  QDP_M_eq_funcia(pU, projectU_site, (void *)Uq, QDP_all);
  QDP_reset_M(U);
  dtime += QOP_time();
  info->final_sec = dtime;
  info->final_flop = 0.0;
#undef NC
}

void
QOP_projectU_deriv_qdp(QOP_info_t *info,
		       QDP_ColorMatrix *f,
		       QDP_ColorMatrix *pU,
		       QDP_ColorMatrix *U,
		       QDP_ColorMatrix *chain)
{
  // TODO
  info->final_sec = 0;
  info->final_flop = 0;
}
