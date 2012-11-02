/****** hisq_force_p.c  -- ******************/
/***AG: wrapper for function in hisq_force_fnmat_p.c***/ 

#include <qop_internal.h>

/*

  Multi-term HISQ force.  

  The coef parameter includes a specification of the number of Naik
  terms and their unique Naik epsilon corrections.  The zeroth Naik
  epsilon is required to be zero.

  The input random color vector fields are arranged as a flat array
  with logical groupings specified by n_orders_naik.  The indexing of
  the groups corresponds to the indexing of the Naik epsilon.  The
  zeroth group thus corresponds to epsilon = 0, the group with index
  1, to the first nonzero epsilon, etc.  The integer n_orders_naik[k]
  specifies the number of terms in the in_pt array belonging to the
  kth Naik epsilon.  The terms in in_pt must be arranged in sequence
  according to the indexing k.  The total number of fields in the
  in_pt array must equal the sum of n_orders_naik[k] over all unique
  Naik epsilons.

 */

void
QOP_hisq_deriv_multi_qdp(QOP_info_t *info, 
			 QOP_FermionLinksHisq *flh,
			 QDP_ColorMatrix *deriv[],
			 QOP_hisq_coeffs_t *coef,
			 REAL epsv[], 
			 QDP_ColorVector *in_pt[], 
			 int *n_orders_naik,
			 int doLastScale)
{
  HISQ_FORCE_BEGIN;

  QOP_hisq_deriv_multi_fnmat2_qdp(info, flh, deriv, coef, epsv, in_pt,
				  n_orders_naik/*, doLastScale*/);

  HISQ_FORCE_END;
}

void
QOP_hisq_force_multi_qdp(QOP_info_t *info, 
			 QOP_FermionLinksHisq *flh,
			 QDP_ColorMatrix *force[],
			 QOP_hisq_coeffs_t *coef,
			 REAL epsv[],
			 QDP_ColorVector *in_pt[], 
			 int *n_orders_naik)
{
  HISQ_FORCE_BEGIN;

  QOP_hisq_force_multi_fnmat2_qdp(info, flh, force, coef, epsv, in_pt, n_orders_naik);

  HISQ_FORCE_END;
}

void
QOP_hisq_force_multi(QOP_info_t *info, 
		     QOP_FermionLinksHisq *flh,
		     QOP_Force *force, 
		     QOP_hisq_coeffs_t *coef,
		     REAL epsv[], 
		     QOP_ColorVector *in_pt[], 
		     int *n_orders_naik)
{
  HISQ_FORCE_BEGIN;

  int n_naiks = coef->n_naiks;
  int nterms = 0;
  for(int inaik = 0; inaik < n_naiks; inaik++)
    nterms += n_orders_naik[inaik];

  QDP_ColorVector *x[nterms];
  for(int i=0; i<nterms; i++) x[i] = in_pt[i]->cv;

  QOP_hisq_force_multi_qdp(info, flh, force->force, coef, epsv, x, n_orders_naik);

  HISQ_FORCE_END;
}
