#define QDP_Precision 'F'
#define QDP_Nc 3
#define QOP_Precision 'F'
#define QOP_Nc 3

#if 0
/* Old version did unitarization in single precision */
#include "hisq_reunit_utilities_p.c"
#else
/* Do unitarization in double precision, regardless */

#include <stdlib.h>
#include <math.h>
#include <qop_internal.h>

int
QOP_u3_un_analytic(QOP_info_t *info, QLA_ColorMatrix *V, QLA_ColorMatrix *W )
{
  QLA_D_ColorMatrix VD, WD;
  int svd_calls;
  QLA_DF_M_eq_M(&VD, V);
  svd_calls = QOP_D3_u3_un_analytic(info, &VD, &WD);
  QLA_FD_M_eq_M(W, &WD);
  return svd_calls;
}

/* Convert 4th rank tensor from double to single precision */
static void 
FD3_T4_eq_T4(QLA_ColorTensor4 *t, QLA_D3_ColorTensor4 *tD)
{
  for(int i0 = 0; i0 < 3; i0++)
    for(int i1 = 0; i1 < 3; i1++)
      for(int i2 = 0; i2 < 3; i2++)
	for(int i3 = 0; i3 < 3; i3++)
	  QLA_FD_C_eq_C(&t->t4[i0][i1][i2][i3], &tD->t4[i0][i1][i2][i3]);
}

void
QOP_u3_un_der_analytic(QOP_info_t *info, QLA_ColorMatrix *V, 
		       QLA_ColorTensor4 *dwdv, QLA_ColorTensor4 *dwdagdv,
		       int *svd_calls, int *ff_counter)
{
  QLA_D_ColorMatrix VD;
  QLA_D3_ColorTensor4 dwdvD, dwdagdvD;
  QLA_DF_M_eq_M(&VD, V);
  QOP_D3_u3_un_der_analytic(info, &VD, &dwdvD, &dwdagdvD, svd_calls, ff_counter );
  FD3_T4_eq_T4(dwdv, &dwdvD);
  FD3_T4_eq_T4(dwdagdv, &dwdagdvD);
}

#endif
