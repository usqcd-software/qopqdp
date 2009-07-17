/****** hisq_force_p.c  -- ******************/
/***AG: wrapper for function in hisq_force_fnmat_p.c***/ 

#include <qop_internal.h>

void 
QOPPC(hisq_force_multi_wrapper_fnmat)(QOP_info_t *info,  
				      QOP_GaugeField *UGauge,
				      QOP_GaugeField *VGauge,
				      QOP_GaugeField *WGauge,
				      QOP_Force *Force, 
				      QOP_hisq_coeffs_t *hisq_coeff,
				      REAL *epsv,
				      QOP_ColorVector *in_pt[], 
				      int nterms,
				      int n_naiks,
				      int n_order_naik_total,
				      int *n_orders_naik,
				      REAL *eps_naik);

void
QOPPC(hisq_force_multi)(QOP_info_t *info, 
			QOP_GaugeField *Ugauge, 
			QOP_GaugeField *Vgauge,
			QOP_GaugeField *Wgauge,
			QOP_Force *force, 
			QOP_hisq_coeffs_t *coef,
			REAL epsv[], 
			QOP_ColorVector *in_pt[], 
			int nsrc,
			int n_naiks,
			int n_order_naik_total,
			int *n_orders_naik,
			REAL *eps_naik)
{

  //HISQ_FORCE_BEGIN;

  QOPPC(hisq_force_multi_wrapper_fnmat)(info, Ugauge, Vgauge, Wgauge, force, 
					coef, epsv, in_pt, nsrc,
					n_naiks, n_order_naik_total,
					n_orders_naik, eps_naik);

  //HISQ_FORCE_END;
}

