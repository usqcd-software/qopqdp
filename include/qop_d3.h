// DO NOT EDIT
// generated from qop_pc.h
#ifndef _QOP_D3_H
#define _QOP_D3_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QOP_D3_ColorVector_struct   QOP_D3_ColorVector;
typedef struct QOP_D3_DiracFermion_struct  QOP_D3_DiracFermion;
typedef struct QOP_D3_GaugeField_struct    QOP_D3_GaugeField;
typedef struct QOP_D3_Force_struct         QOP_D3_Force;

typedef struct QOP_D3_FermionLinksAsqtad_struct  QOP_D3_FermionLinksAsqtad;
typedef struct QOP_D3_FermionLinksHisq_struct    QOP_D3_FermionLinksHisq;
typedef struct QOP_D3_FermionLinksWilson_struct  QOP_D3_FermionLinksWilson;
typedef struct QOP_D3_FermionLinksDW_struct      QOP_D3_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

#define QOP_D3_qla_type_V QLA_D3_ColorVector
#define QOP_D3_qla_type_D QLA_D3_DiracFermion
#define QOP_D3_qla_type_M QLA_D3_ColorMatrix
#define QOP_D3_raw_size(T) (QDP_sites_on_node*sizeof(QOP_D3_qla_type_##T))
#define QOP_D3_raw_size_V(evenodd) QOP_D3_raw_size(V)
#define QOP_D3_raw_size_D(evenodd) QOP_D3_raw_size(D)
#define QOP_D3_raw_size_G(evenodd) QOP_D3_raw_size(M)
#define QOP_D3_raw_size_F(evenodd) QOP_D3_raw_size(M)

#define QOP_D3_elem(T, raw, i, ...) QLA_D3_elem_##T(((QOP_D3_qla_type_##T *)raw)[i], __ARGV__)
#define QOP_D3_set(T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_D3_elem(T, raw, i, __ARGV__), re, im)
#define QOP_D3_raw_set_V(raw, evenodd, i, ic, re, im) QOP_D3_set(V, raw, i, re, im, ic)
#define QOP_D3_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_D3_set(D, raw, i, re, im, ic, is)
#define QOP_D3_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_D3_set(M, raw, i, re, im, ic, jc)
#define QOP_D3_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_D3_set(M, raw, i, re, im, ic, jc)
#define QOP_D3_get(T, re, im, raw, i, ...) {			\
    QLA_D_Complex _c = QOP_D3_elem(T, raw, i, __ARGV__);	\
    re = QLA_real(_c); im = QLA_imag(_c);			\
  }
#define QOP_D3_raw_get_V(re, im, raw, evenodd, i, ic)     QOP_D3_set(V, re, im, raw, i, ic)
#define QOP_D3_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_D3_set(D, re, im, raw, i, ic, is)
#define QOP_D3_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_D3_set(M, re, im, raw, i, ic, jc)
#define QOP_D3_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_D3_set(M, re, im, raw, i, ic, jc)

/* create a QOP field with a copy of the raw source field */
QOP_D3_ColorVector  *QOP_D3_create_V_from_raw(QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_D3_DiracFermion *QOP_D3_create_D_from_raw(QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_D3_GaugeField   *QOP_D3_create_G_from_raw(QOP_D_Real *links[], QOP_evenodd_t evenodd);
QOP_D3_Force        *QOP_D3_create_F_from_raw(QOP_D_Real *force[], QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_D3_extract_V_to_raw(QOP_D_Real *dest, QOP_D3_ColorVector *src, QOP_evenodd_t evenodd);
void QOP_D3_extract_D_to_raw(QOP_D_Real *dest, QOP_D3_DiracFermion *src, QOP_evenodd_t evenodd);
void QOP_D3_extract_G_to_raw(QOP_D_Real *dest[], QOP_D3_GaugeField *src, QOP_evenodd_t evenodd);
void QOP_D3_extract_F_to_raw(QOP_D_Real *dest[], QOP_D3_Force *src, QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_D3_destroy_V(QOP_D3_ColorVector *field);
void QOP_D3_destroy_D(QOP_D3_DiracFermion *field);
void QOP_D3_destroy_G(QOP_D3_GaugeField *field);
void QOP_D3_destroy_F(QOP_D3_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_D3_ColorVector  *QOP_D3_convert_V_from_raw(QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_D3_DiracFermion *QOP_D3_convert_D_from_raw(QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_D3_GaugeField   *QOP_D3_convert_G_from_raw(QOP_D_Real *links[], QOP_evenodd_t evenodd);
QOP_D3_Force        *QOP_D3_convert_F_from_raw(QOP_D_Real *force[], QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
QOP_D_Real  *QOP_D3_convert_V_to_raw(QOP_D3_ColorVector *src, QOP_evenodd_t evenodd);
QOP_D_Real  *QOP_D3_convert_D_to_raw(QOP_D3_DiracFermion *src, QOP_evenodd_t evenodd);
QOP_D_Real **QOP_D3_convert_G_to_raw(QOP_D3_GaugeField *src, QOP_evenodd_t evenodd);
QOP_D_Real **QOP_D3_convert_F_to_raw(QOP_D3_Force *src, QOP_evenodd_t evenodd);

QOP_D3_ColorVector  *QOP_D3_create_V_from_qdp(QDP_D3_ColorVector *src);
QOP_D3_DiracFermion *QOP_D3_create_D_from_qdp(QDP_D3_DiracFermion *src);
QOP_D3_GaugeField   *QOP_D3_create_G_from_qdp(QDP_D3_ColorMatrix *src[]);
QOP_D3_Force        *QOP_D3_create_F_from_qdp(QDP_D3_ColorMatrix *src[]);

void QOP_D3_extract_V_to_qdp(QDP_D3_ColorVector *d, QOP_D3_ColorVector *src);
void QOP_D3_extract_D_to_qdp(QDP_D3_DiracFermion *d, QOP_D3_DiracFermion *src);
void QOP_D3_extract_G_to_qdp(QDP_D3_ColorMatrix *d[], QOP_D3_GaugeField *src);
void QOP_D3_extract_F_to_qdp(QDP_D3_ColorMatrix *d[], QOP_D3_Force *src);

QOP_D3_ColorVector  *QOP_D3_convert_V_from_qdp(QDP_D3_ColorVector *src);
QOP_D3_DiracFermion *QOP_D3_convert_D_from_qdp(QDP_D3_DiracFermion *src);
QOP_D3_GaugeField   *QOP_D3_convert_G_from_qdp(QDP_D3_ColorMatrix *src[]);
QOP_D3_Force        *QOP_D3_convert_F_from_qdp(QDP_D3_ColorMatrix *src[]);

QDP_D3_ColorVector   *QOP_D3_convert_V_to_qdp(QOP_D3_ColorVector *src);
QDP_D3_DiracFermion  *QOP_D3_convert_D_to_qdp(QOP_D3_DiracFermion *src);
QDP_D3_ColorMatrix  **QOP_D3_convert_G_to_qdp(QOP_D3_GaugeField *src);
QDP_D3_ColorMatrix  **QOP_D3_convert_F_to_qdp(QOP_D3_Force *src);


  /********************/
  /*  Gauge routines  */
  /********************/

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */
void QOP_D3_rephase_G(QOP_D3_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_D3_rephase_G_qdp(QDP_D3_ColorMatrix *links[],
			  int *r0,
			  QOP_bc_t *bc,
			  QOP_staggered_sign_t *sign);

void QOP_D3_smear_fat7l_qdp(QOP_info_t *info, QDP_D3_ColorMatrix *sg[],
			    QDP_D3_ColorMatrix *g[],
			    QOP_asqtad_coeffs_t *coeffs);

void QOP_D3_gauge_deriv_multi_qdp(QOP_info_t *info, QDP_D3_ColorMatrix *deriv[],
				  QOP_D3_GaugeField *g[], QDP_D3_ColorMatrix **chain[],
				  int n, int doLastScale);

void QOP_D3_gauge_force_multi_qdp(QOP_info_t *info, QDP_D3_ColorMatrix *f[],
				  QOP_D3_GaugeField *g[], QDP_D3_ColorMatrix **chain[], int n);

void QOP_D3_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_D3_GaugeField *gauge,
					QOP_D_Real *acts, QOP_D_Real *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_D3_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info,
					   QOP_D3_GaugeField *gauge,
					   QDP_D3_ColorMatrix *deriv[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_D_Real eps, int doLastScale);

void QOP_D3_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, 
					   QOP_D3_GaugeField *gauge, 
					   QDP_D3_ColorMatrix *force[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_D_Real eps);

void QOP_D3_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_D3_GaugeField *gauge, 
				       QOP_D3_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       QOP_D_Real eps);


  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_raw(QOP_D_Real *fatlinks[], QOP_D_Real *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_D3_GaugeField *gauge);

void QOP_D3_asqtad_extract_L_to_raw(QOP_D_Real *fatlinks[], QOP_D_Real *longlinks[],
				    QOP_D3_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_D3_asqtad_destroy_L(QOP_D3_FermionLinksAsqtad *field);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_convert_L_from_raw(QOP_D_Real *fatlinks[], QOP_D_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_D3_asqtad_convert_L_to_raw(QOP_D_Real ***fatlinks, QOP_D_Real ***longlinks,
				    QOP_D3_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_D3_asqtad_load_L_from_raw(QOP_D3_FermionLinksAsqtad *asqtad,
				   QOP_D_Real *fatlinks[], QOP_D_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_D3_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_D3_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

void QOP_D3_asqtad_rephase_L(QOP_D3_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_qdp(QDP_D3_ColorMatrix *fatlinks[],
				  QDP_D3_ColorMatrix *longlinks[]);

void QOP_D3_asqtad_extract_L_to_qdp(QDP_D3_ColorMatrix *fatlinks[],
				    QDP_D3_ColorMatrix *longlinks[],
				    QOP_D3_FermionLinksAsqtad *src);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_convert_L_from_qdp(QDP_D3_ColorMatrix *fatlinks[],
				   QDP_D3_ColorMatrix *longlinks[]);

void QOP_D3_asqtad_convert_L_to_qdp(QDP_D3_ColorMatrix ***fatlinks,
				    QDP_D3_ColorMatrix ***longlinks,
				    QOP_D3_FermionLinksAsqtad *src);

void QOP_D3_asqtad_load_L_from_qdp(QOP_D3_FermionLinksAsqtad *asqtad,
				   QDP_D3_ColorMatrix *fatlinks[],
				   QDP_D3_ColorMatrix *longlinks[]);

void QOP_D3_asqtad_rephase_field_L_qdp(QOP_D3_FermionLinksAsqtad *fla,
				       QDP_D_Complex *fatphase[],
				       QDP_D_Complex *longphase[]);


  /* inverter routines */

void QOP_D3_asqtad_dslash(QOP_info_t *info,
			  QOP_D3_FermionLinksAsqtad *asqtad,
			  QOP_D_Real mass,
			  QOP_D3_ColorVector *out,
			  QOP_D3_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_D3_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_D3_ColorVector *out,
			      QOP_D3_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_D3_asqtad_diaginv(QOP_info_t *info,
			   QOP_D3_FermionLinksAsqtad *asqtad,
			   QOP_D_Real mass,
			   QOP_D3_ColorVector *out,
			   QOP_D3_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_D3_asqtad_invert(QOP_info_t *info,
			  QOP_D3_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real mass,
			  QOP_D3_ColorVector *out_pt,
			  QOP_D3_ColorVector *in_pt);

void QOP_D3_asqtad_invert_threaded(QOP_info_t *info,
				   QOP_D3_FermionLinksAsqtad *asqtad,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg,
				   QOP_D_Real mass,
				   QOP_D3_ColorVector *out_pt,
				   QOP_D3_ColorVector *in_pt,
				   int nthreads);

void QOP_D3_asqtad_invert_multi(QOP_info_t *info,
				QOP_D3_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *masses[],
				int nmass[],
				QOP_D3_ColorVector **out_pt[],
				QOP_D3_ColorVector *in_pt[],
				int nsrc);

void QOP_D3_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      QOP_D_Real mass,
			      QDP_D3_ColorVector *out,
			      QDP_D3_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_D3_asqtad_dslash_dir_qdp(QOP_info_t *info,
				  QOP_D3_FermionLinksAsqtad *asqtad,
				  int dir, int fb,
				  double wtfat, double wtlong,
				  QDP_D3_ColorVector *out,
				  QDP_D3_ColorVector *in,
				  QOP_evenodd_t eo_out);

void QOP_D3_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_D3_FermionLinksAsqtad *asqtad,
			       QOP_D_Real mass,
			       QDP_D3_ColorVector *out,
			       QDP_D3_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_D3_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_D_Real mass,
			      QDP_D3_ColorVector *out,
			      QDP_D3_ColorVector *in);

void QOP_D3_asqtad_invert_threaded_qdp(QOP_info_t *info,
				       QOP_D3_FermionLinksAsqtad *asqtad,
				       QOP_invert_arg_t *inv_arg,
				       QOP_resid_arg_t *res_arg,
				       QOP_D_Real mass,
				       QDP_D3_ColorVector *out,
				       QDP_D3_ColorVector *in,
				       int nthreads);

void QOP_D3_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_D3_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_D_Real *masses[],
				    int nmass[],
				    QDP_D3_ColorVector **out[],
				    QDP_D3_ColorVector *in[],
				    int nsrc);

void QOP_D3_asqtad_get_eigcg(QOP_D3_FermionLinksAsqtad *asqtad,
			     QDP_D3_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

  /* fermion force routines */

void QOP_D3_asqtad_deriv(QOP_info_t *info, QDP_D3_ColorMatrix *gauge[],
			 QDP_D3_ColorMatrix *force[],
			 QOP_asqtad_coeffs_t *coef,
			 QDP_D3_ColorMatrix *mid_fat[],
			 QDP_D3_ColorMatrix *mid_naik[]);

void QOP_D3_asqtad_force(QOP_info_t *info,
			 QOP_D3_GaugeField *gauge,
			 QOP_D3_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 QOP_D_Real eps,
			 QOP_D3_ColorVector *in_pt);

void QOP_D3_asqtad_force_multi(QOP_info_t *info,
			       QOP_D3_GaugeField *gauge,
			       QOP_D3_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       QOP_D_Real eps[],
			       QOP_D3_ColorVector *in_pt[],
			       int nsrc);

void QOP_D3_asqtad_force_multi_qdp(QOP_info_t *info,
				   QOP_D3_GaugeField *gauge,
				   QDP_D3_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_D_Real eps[],
				   QDP_D3_ColorVector *in_pt[],
				   int nsrc);

  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* single precision */

QOP_D3_FermionLinksHisq *
  QOP_D3_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_D3_GaugeField *gauge);

void QOP_D3_hisq_destroy_L(QOP_D3_FermionLinksHisq *field);

QOP_D3_FermionLinksAsqtad **
  QOP_D3_get_asqtad_links_from_hisq(QOP_D3_FermionLinksHisq *hl);
  
QOP_D3_FermionLinksAsqtad *
  QOP_D3_get_asqtad_deps_links_from_hisq(QOP_D3_FermionLinksHisq *hl);

  /* fermion force routines */

void QOP_D3_hisq_force_multi(QOP_info_t *info,
			     QOP_D3_FermionLinksHisq *flh,
			     QOP_D3_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     QOP_D_Real eps[],
			     QOP_D3_ColorVector *in_pt[],
			     int *n_orders_naik);

void QOP_D3_hisq_deriv_multi_qdp(QOP_info_t *info,
				 QOP_D3_FermionLinksHisq *flh,
				 QDP_D3_ColorMatrix *deriv[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_D_Real eps[],
				 QDP_D3_ColorVector *in_pt[],
				 int *n_orders_naik,
				 int doLastScale);

void QOP_D3_hisq_force_multi_qdp(QOP_info_t *info,
				 QOP_D3_FermionLinksHisq *flh,
				 QDP_D3_ColorMatrix *force[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_D_Real eps[],
				 QDP_D3_ColorVector *in_pt[],
				 int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_raw(QOP_D_Real *links[], QOP_D_Real *clov,
				  QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_D3_GaugeField *gauge);

void QOP_D3_wilson_extract_L_to_raw(QOP_D_Real *links[], QOP_D_Real *clov,
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_D3_wilson_destroy_L(QOP_D3_FermionLinksWilson *field);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_raw(QOP_D_Real *links[], QOP_D_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_D3_wilson_convert_L_to_raw(QOP_D_Real ***links, QOP_D_Real **clov,
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_wilson_convert_L_to_G(QOP_D3_FermionLinksWilson *links);

void QOP_D3_wilson_load_L_from_raw(QOP_D3_FermionLinksWilson *wilson,
				   QOP_D_Real *links[], QOP_D_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_D3_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_D3_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_qdp(QDP_D3_ColorMatrix *links[],
				  QDP_D3_DiracPropagator *clov);

void QOP_D3_wilson_extract_L_to_qdp(QDP_D3_ColorMatrix *links[],
				    QDP_D3_DiracPropagator *clov,
				    QOP_D3_FermionLinksWilson *src);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_qdp(QDP_D3_ColorMatrix *links[],
				   QDP_D3_DiracPropagator *clov);

void QOP_D3_wilson_convert_L_to_qdp(QDP_D3_ColorMatrix ***links,
				    QDP_D3_DiracPropagator **clov,
				    QOP_D3_FermionLinksWilson *src);

void QOP_D3_wilson_load_L_from_qdp(QOP_D3_FermionLinksWilson *wilson,
				   QDP_D3_ColorMatrix *links[],
				   QDP_D3_DiracPropagator *clov);


  /* inverter routines */

void QOP_D3_wilson_dslash(QOP_info_t *info,
			  QOP_D3_FermionLinksWilson *flw,
			  QOP_D_Real kappa,
			  int sign,
			  QOP_D3_DiracFermion *out,
			  QOP_D3_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_D3_wilson_diaginv(QOP_info_t *info,
			   QOP_D3_FermionLinksWilson *flw,
			   QOP_D_Real kappa,
			   QOP_D3_DiracFermion *out,
			   QOP_D3_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_D3_wilson_invert(QOP_info_t *info,
			  QOP_D3_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real kappa,
			  QOP_D3_DiracFermion *out_pt,
			  QOP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_invert_multi(QOP_info_t *info,
				QOP_D3_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *kappas[],
				int nkappa[],
				QOP_D3_DiracFermion **out_pt[],
				QOP_D3_DiracFermion *in_pt[],
				int nsrc);

void QOP_D3_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksWilson *flw,
			      QOP_D_Real kappa,
			      int sign,
			      QDP_D3_DiracFermion *out,
			      QDP_D3_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_D3_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_D3_FermionLinksWilson *flw,
			       QOP_D_Real kappa,
			       QDP_D3_DiracFermion *out,
			       QDP_D3_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_D3_wilson_invert_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_D_Real kappa,
			      QDP_D3_DiracFermion *out_pt,
			      QDP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_D3_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_D_Real *kappas[],
				    int nkappa[],
				    QDP_D3_DiracFermion **out_pt[],
				    QDP_D3_DiracFermion *in_pt[],
				    int nsrc);

void QOP_D3_wilson_invert_ne_qdp(QOP_info_t *info,
				 QOP_D3_FermionLinksWilson *flw,
				 QOP_invert_arg_t *inv_arg,
				 QOP_resid_arg_t *res_arg,
				 QOP_D_Real kappa,
				 QDP_D3_DiracFermion *out,
				 QDP_D3_DiracFermion *in);

  /* fermion force routines */

void QOP_D3_wilson_force(QOP_info_t *info,
			 QOP_D3_GaugeField *gauge,
			 QOP_D3_Force *force,
			 QOP_wilson_coeffs_t *coeffs,
			 QOP_D_Real eps,
			 QOP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_deriv_multi_qdp(QOP_info_t *info,
				   QOP_D3_FermionLinksWilson *flw,
				   QDP_D3_ColorMatrix *deriv[],
				   QOP_D_Real eps[],
				   QDP_D3_DiracFermion *x[],
				   QDP_D3_DiracFermion *y[],
				   int n);

void QOP_D3_wilson_force_multi(QOP_info_t *info,
			       QOP_D3_GaugeField *gauge,
			       QOP_D3_Force *force,
			       QOP_wilson_coeffs_t *coef,
			       QOP_D_Real eps[],
			       QOP_D3_DiracFermion *in_pt[],
			       int nsrc);

void QOP_D3_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
					QOP_D3_FermionLinksWilson *flw,
					QDP_D3_ColorMatrix *deriv[],
					QOP_D_Real kappa[],
					QOP_D_Real eps[],
					QDP_D3_DiracFermion *x[],
					QDP_D3_DiracFermion *y[],
					int n);

void QOP_D3_wilson_force_prec_multi_qdp(QOP_info_t *info,
					QOP_D3_FermionLinksWilson *flw,
					QDP_D3_ColorMatrix *force[],
					QOP_D_Real kappa[],
					QOP_D_Real eps[],
					QDP_D3_DiracFermion *x[],
					QDP_D3_DiracFermion *y[],
					int n);

  // new fermilab action IFLA -- added by bugra --------------- :

void QOP_D3_wilson_ifla_dslash(QOP_info_t *info,
			       QOP_D3_FermionLinksWilson *flw,
			       QOP_D_Real kappa,
			       int sign,
			       QOP_wilson_ifla_coeffs_t *coeffs,
			       QOP_D3_DiracFermion *out,
			       QOP_D3_DiracFermion *in,
			       QOP_evenodd_t eo_out,
			       QOP_evenodd_t eo_in);

void QOP_D3_wilson_ifla_dslash_qdp(QOP_info_t *info,
				   QOP_D3_FermionLinksWilson *flw,
				   QOP_D_Real kappa,
				   int sign,
				   QOP_wilson_ifla_coeffs_t *coeffs,
				   QDP_D3_DiracFermion *out,
				   QDP_D3_DiracFermion *in,
				   QOP_evenodd_t eo_out,
				   QOP_evenodd_t eo_in);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_raw(QOP_D_Real *links[], QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_D3_GaugeField *gauge);

void QOP_D3_dw_extract_L_to_raw(QOP_D_Real *links[],
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_D3_dw_destroy_L(QOP_D3_FermionLinksDW *field);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_raw(QOP_D_Real *links[], QOP_evenodd_t evenodd);

void QOP_D3_dw_convert_L_to_raw(QOP_D_Real ***links,
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_dw_convert_L_to_G(QOP_D3_FermionLinksDW *links);

void QOP_D3_dw_load_L_from_raw(QOP_D3_FermionLinksDW *dw,
			       QOP_D_Real *links[], QOP_evenodd_t evenodd);

void QOP_D3_dw_load_L_from_G(QOP_info_t *info,
			     QOP_D3_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_D3_GaugeField *gauge);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_qdp(QDP_D3_ColorMatrix *links[]);

void QOP_D3_dw_extract_L_to_qdp(QDP_D3_ColorMatrix *links[],
				QOP_D3_FermionLinksDW *src);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_qdp(QDP_D3_ColorMatrix *links[]);

void QOP_D3_dw_convert_L_to_qdp(QDP_D3_ColorMatrix ***links,
				QOP_D3_FermionLinksDW *src);

void QOP_D3_dw_load_L_from_qdp(QOP_D3_FermionLinksDW *dw,
			       QDP_D3_ColorMatrix *links[]);

  /* inverter routines */

void QOP_D3_dw_dslash(QOP_info_t *info,
		      QOP_D3_FermionLinksDW *links,
		      QOP_D_Real M5,
		      QOP_D_Real m,
		      int sign,
		      QOP_D3_DiracFermion *out_pt[],
		      QOP_D3_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_D3_dw_dslash2(QOP_info_t *info,
		       QOP_D3_FermionLinksDW *links,
		       QOP_D_Real M5,
		       QOP_D_Real m,
		       QOP_D3_DiracFermion *out_pt[],
		       QOP_D3_DiracFermion *in_pt[],
		       int Ls,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in);

void QOP_D3_dw_invert(QOP_info_t *info,
		      QOP_D3_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QOP_D_Real M5,
		      QOP_D_Real m,
		      QOP_D3_DiracFermion *out_pt[],
		      QOP_D3_DiracFermion *in_pt[],
		      int Ls);

void QOP_D3_dw_invert_multi(QOP_info_t *info,
			    QOP_D3_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_D_Real *M5[],
			    QOP_D_Real *m[],
			    int nmass[],
			    QOP_D3_DiracFermion ***out_pt[],
			    QOP_D3_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_D3_dw_dslash_qdp(QOP_info_t *info,
			  QOP_D3_FermionLinksDW *links,
			  QOP_D_Real M5, 
			  QOP_D_Real m,
			  int sign,
			  QDP_D3_DiracFermion *out_pt[],
			  QDP_D3_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_D3_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_D3_FermionLinksDW *links,
			   QOP_D_Real M5,
			   QOP_D_Real m,
			   QDP_D3_DiracFermion *out_pt[],
			   QDP_D3_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);

void QOP_D3_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_D3_FermionLinksDW *fldw,
			   QOP_D_Real M5,
			   QOP_D_Real m,
			   QDP_D3_DiracFermion **out,
			   QDP_D3_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_D3_dw_invert_qdp(QOP_info_t *info,
			  QOP_D3_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real M5,
			  QOP_D_Real m,
			  QDP_D3_DiracFermion *out[],
			  QDP_D3_DiracFermion *in[],
			  int Ls);

void QOP_D3_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_D3_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *M5[],
				QOP_D_Real *m[],
				int nmass[],
				QDP_D3_DiracFermion ***out[],
				QDP_D3_DiracFermion **in[],
				int nsrc,
				int Ls);

 /* fermion force routines */

void QOP_D3_dw_force(QOP_info_t *info,
		     QOP_D3_GaugeField *gauge,
		     QOP_D3_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     QOP_D_Real eps,
		     QOP_D3_DiracFermion *in_pt);

void QOP_D3_dw_force_multi(QOP_info_t *info,
			   QOP_D3_GaugeField *gauge,
			   QOP_D3_Force *force,
			   QOP_dw_coeffs_t *coef,
			   QOP_D_Real eps[],
			   QOP_D3_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == 'D'
#  if QOP_Colors == 3
#    include <qop_d3_generic.h>
#  endif
#  include <qop_d3_precision_generic.h>
#endif
#if QOP_Colors == 3
#  include <qop_d3_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_D3_H */
