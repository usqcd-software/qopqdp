#include <qop_internal.h>
#if QOP_Colors == 2
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
u3reunit_site_d(QLA_D_ColorMatrix *Ur, QLA_D_ColorMatrix *U)
{
  QLA_D_ColorMatrix M1, M2;
  QLA_D_M_eq_Ma_times_M(&M1, U, U);
  QLA_D_M_eq_invsqrt_M(&M2, &M1);
  QLA_D_M_eq_M_times_M(Ur, U, &M2);
}

static void
u3reunit_site(QLA_ColorMatrix *Ur, int i, void *args)
{
  QLA_ColorMatrix *U = &((QLA_ColorMatrix *)args)[i];
#if QLA_Precision == 'F'
  QLA_D_ColorMatrix Ud, Urd;
  QLA_DF_M_eq_M(&Ud, U);
  u3reunit_site_d(&Urd, &Ud);
  QLA_FD_M_eq_M(Ur, &Urd);
#else
  u3reunit_site_d(Ur, U);
#endif
}

// requires QLA >= 1.7.0 and QDP >= 1.9.0
void
QOP_u3reunit(QOP_info_t *info, QDP_ColorMatrix *U, QDP_ColorMatrix *Ur)
{
  double dtime = -QOP_time();
  QLA_ColorMatrix *Uq = QDP_expose_M(U);
  QDP_M_eq_funcia(Ur, u3reunit_site, (void *)Uq, QDP_all);
  QDP_reset_M(U);
  dtime += QOP_time();
  info->final_sec = dtime;
  info->final_flop = 0.0;
}

void
QOP_hisq_force_multi_reunit(QOP_info_t *info,
			    QDP_ColorMatrix *V[4],
			    QDP_ColorMatrix *Force[4],
			    QDP_ColorMatrix *Force_old[4])
{
  // TODO
  info->final_sec = 0;
  info->final_flop = 0;
}
