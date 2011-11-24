#define QDP_Precision 1
#define QOP_Precision 1


//#include "hisq_reunit_utilities_p.c"  /* Old version did unitarization in single precision */

#include <stdlib.h>
#include <math.h>
#include <qop_internal.h>

/* Do unitarization in double precision, regardless */

int QOP_F3_u3_un_analytic( QOP_info_t *info, QLA_F3_ColorMatrix *V, QLA_F3_ColorMatrix *W ){
  QLA_D3_ColorMatrix VD, WD;
  int svd_calls;

  QLA_DF3_M_eq_M(&VD, V);

  svd_calls = QOP_D3_u3_un_analytic(info, &VD, &WD);

  QLA_FD3_M_eq_M(W, &WD);

  return svd_calls;

}

/* Convert 4th rank tensor from double to single precision */

static void FD3_T4_eq_T4(QLA_F3_ColorTensor4 *t, QLA_D3_ColorTensor4 *tD){
  int i0,i1,i2,i3;

  for(i0 = 0; i0 < 3; i0++)
    for(i1 = 0; i1 < 3; i1++)
      for(i2 = 0; i2 < 3; i2++)
	for(i3 = 0; i3 < 3; i3++)
	  QLA_FD_C_eq_C(&t->t4[i0][i1][i2][i3], &tD->t4[i0][i1][i2][i3]);
}

void QOP_F3_u3_un_der_analytic( QOP_info_t *info, QLA_F3_ColorMatrix *V, 
				QLA_F3_ColorTensor4 *dwdv, QLA_F3_ColorTensor4 *dwdagdv,
				int *svd_calls, int *ff_counter){
  QLA_D3_ColorMatrix VD;
  QLA_D3_ColorTensor4 dwdvD, dwdagdvD;

  QLA_DF3_M_eq_M(&VD, V);

  QOP_D3_u3_un_der_analytic(info, &VD, &dwdvD, &dwdagdvD, svd_calls, ff_counter );

  FD3_T4_eq_T4(dwdv, &dwdvD);
  FD3_T4_eq_T4(dwdagdv, &dwdagdvD);

}

// The following is needed only for testing
// Determinant of 3x3 complex matrix
QLA_F_Complex QOP_F3_su3_mat_det( QLA_F3_ColorMatrix *U) {
  QLA_F_Complex a, b, m0, m1, m2, cdet;
#define Uelem(a,b) QLA_elem_M(*U,a,b)

  // brute-force calculation
  QLA_c_eq_c_times_c( a, Uelem(1,1), Uelem(2,2) );
  QLA_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,1) );
  QLA_c_eq_c_minus_c( m0, a, b );
  
  QLA_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,2) );
  QLA_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,0) );
  QLA_c_eq_c_minus_c( m1, a, b );
  
  QLA_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,1) );
  QLA_c_eq_c_times_c( b, Uelem(1,1), Uelem(2,0) );
  QLA_c_eq_c_minus_c( m2, a, b );

  QLA_c_eq_c_times_c( a, Uelem(0,0), m0 );
  QLA_c_eq_c_times_c( b, Uelem(0,1), m1 );
  QLA_c_eq_c_minus_c( cdet, a, b );
  QLA_c_eq_c_times_c( a, Uelem(0,2), m2 );
  QLA_c_eq_c_plus_c( cdet, cdet, a );
  
  
  return cdet;
}

