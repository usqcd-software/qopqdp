// DO NOT EDIT
// generated from qop_pc.h
#ifndef _QOP_DN_H
#define _QOP_DN_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QOP_DN_ColorVector_struct   QOP_DN_ColorVector;
typedef struct QOP_DN_DiracFermion_struct  QOP_DN_DiracFermion;
typedef struct QOP_DN_GaugeField_struct    QOP_DN_GaugeField;
typedef struct QOP_DN_Force_struct         QOP_DN_Force;

typedef struct QOP_DN_FermionLinksAsqtad_struct  QOP_DN_FermionLinksAsqtad;
typedef struct QOP_DN_FermionLinksHisq_struct    QOP_DN_FermionLinksHisq;
typedef struct QOP_DN_FermionLinksWilson_struct  QOP_DN_FermionLinksWilson;
typedef struct QOP_DN_FermionLinksDW_struct      QOP_DN_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

#define QOP_DN_qla_type_V QLA_DN_ColorVector
#define QOP_DN_qla_type_D QLA_DN_DiracFermion
#define QOP_DN_qla_type_M QLA_DN_ColorMatrix
#define QOP_DN_raw_size(T) (QDP_sites_on_node*sizeof(QOP_DN_qla_type_##T))
#define QOP_DN_raw_size_V(evenodd) QOP_DN_raw_size(V)
#define QOP_DN_raw_size_D(evenodd) QOP_DN_raw_size(D)
#define QOP_DN_raw_size_G(evenodd) QOP_DN_raw_size(M)
#define QOP_DN_raw_size_F(evenodd) QOP_DN_raw_size(M)

#define QOP_DN_elem(T, raw, i, ...) QLA_DN_elem_##T(((QOP_DN_qla_type_##T *)raw)[i], __ARGV__)
#define QOP_DN_set(T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_DN_elem(T, raw, i, __ARGV__), re, im)
#define QOP_DN_raw_set_V(raw, evenodd, i, ic, re, im) QOP_DN_set(V, raw, i, re, im, ic)
#define QOP_DN_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_DN_set(D, raw, i, re, im, ic, is)
#define QOP_DN_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_DN_set(M, raw, i, re, im, ic, jc)
#define QOP_DN_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_DN_set(M, raw, i, re, im, ic, jc)
#define QOP_DN_get(T, re, im, raw, i, ...) {			\
    QLA_D_Complex _c = QOP_DN_elem(T, raw, i, __ARGV__);	\
    re = QLA_real(_c); im = QLA_imag(_c);			\
  }
#define QOP_DN_raw_get_V(re, im, raw, evenodd, i, ic)     QOP_DN_set(V, re, im, raw, i, ic)
#define QOP_DN_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_DN_set(D, re, im, raw, i, ic, is)
#define QOP_DN_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_DN_set(M, re, im, raw, i, ic, jc)
#define QOP_DN_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_DN_set(M, re, im, raw, i, ic, jc)

/* create a QOP field with a copy of the raw source field */
QOP_DN_ColorVector  *QOP_DN_create_V_from_raw(int nc, QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_DN_DiracFermion *QOP_DN_create_D_from_raw(int nc, QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_DN_GaugeField   *QOP_DN_create_G_from_raw(int nc, QOP_D_Real *links[], QOP_evenodd_t evenodd);
QOP_DN_Force        *QOP_DN_create_F_from_raw(int nc, QOP_D_Real *force[], QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_DN_extract_V_to_raw(QOP_D_Real *dest, QOP_DN_ColorVector *src, QOP_evenodd_t evenodd);
void QOP_DN_extract_D_to_raw(QOP_D_Real *dest, QOP_DN_DiracFermion *src, QOP_evenodd_t evenodd);
void QOP_DN_extract_G_to_raw(QOP_D_Real *dest[], QOP_DN_GaugeField *src, QOP_evenodd_t evenodd);
void QOP_DN_extract_F_to_raw(QOP_D_Real *dest[], QOP_DN_Force *src, QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_DN_destroy_V(QOP_DN_ColorVector *field);
void QOP_DN_destroy_D(QOP_DN_DiracFermion *field);
void QOP_DN_destroy_G(QOP_DN_GaugeField *field);
void QOP_DN_destroy_F(QOP_DN_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_DN_ColorVector  *QOP_DN_convert_V_from_raw(int nc, QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_DN_DiracFermion *QOP_DN_convert_D_from_raw(int nc, QOP_D_Real *src, QOP_evenodd_t evenodd);
QOP_DN_GaugeField   *QOP_DN_convert_G_from_raw(int nc, QOP_D_Real *links[], QOP_evenodd_t evenodd);
QOP_DN_Force        *QOP_DN_convert_F_from_raw(int nc, QOP_D_Real *force[], QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
QOP_D_Real  *QOP_DN_convert_V_to_raw(QOP_DN_ColorVector *src, QOP_evenodd_t evenodd);
QOP_D_Real  *QOP_DN_convert_D_to_raw(QOP_DN_DiracFermion *src, QOP_evenodd_t evenodd);
QOP_D_Real **QOP_DN_convert_G_to_raw(QOP_DN_GaugeField *src, QOP_evenodd_t evenodd);
QOP_D_Real **QOP_DN_convert_F_to_raw(QOP_DN_Force *src, QOP_evenodd_t evenodd);

QOP_DN_ColorVector  *QOP_DN_create_V_from_qdp(QDP_DN_ColorVector *src);
QOP_DN_DiracFermion *QOP_DN_create_D_from_qdp(QDP_DN_DiracFermion *src);
QOP_DN_GaugeField   *QOP_DN_create_G_from_qdp(QDP_DN_ColorMatrix *src[]);
QOP_DN_Force        *QOP_DN_create_F_from_qdp(QDP_DN_ColorMatrix *src[]);

void QOP_DN_extract_V_to_qdp(QDP_DN_ColorVector *d, QOP_DN_ColorVector *src);
void QOP_DN_extract_D_to_qdp(QDP_DN_DiracFermion *d, QOP_DN_DiracFermion *src);
void QOP_DN_extract_G_to_qdp(QDP_DN_ColorMatrix *d[], QOP_DN_GaugeField *src);
void QOP_DN_extract_F_to_qdp(QDP_DN_ColorMatrix *d[], QOP_DN_Force *src);

QOP_DN_ColorVector  *QOP_DN_convert_V_from_qdp(QDP_DN_ColorVector *src);
QOP_DN_DiracFermion *QOP_DN_convert_D_from_qdp(QDP_DN_DiracFermion *src);
QOP_DN_GaugeField   *QOP_DN_convert_G_from_qdp(QDP_DN_ColorMatrix *src[]);
QOP_DN_Force        *QOP_DN_convert_F_from_qdp(QDP_DN_ColorMatrix *src[]);

QDP_DN_ColorVector   *QOP_DN_convert_V_to_qdp(QOP_DN_ColorVector *src);
QDP_DN_DiracFermion  *QOP_DN_convert_D_to_qdp(QOP_DN_DiracFermion *src);
QDP_DN_ColorMatrix  **QOP_DN_convert_G_to_qdp(QOP_DN_GaugeField *src);
QDP_DN_ColorMatrix  **QOP_DN_convert_F_to_qdp(QOP_DN_Force *src);


  /********************/
  /*  Gauge routines  */
  /********************/

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */
void QOP_DN_rephase_G(QOP_DN_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_DN_rephase_G_qdp(QDP_DN_ColorMatrix *links[],
			  int *r0,
			  QOP_bc_t *bc,
			  QOP_staggered_sign_t *sign);

void QOP_DN_smear_fat7l_qdp(QOP_info_t *info, QDP_DN_ColorMatrix *sg[],
			    QDP_DN_ColorMatrix *g[],
			    QOP_asqtad_coeffs_t *coeffs);

void QOP_DN_gauge_deriv_multi_qdp(QOP_info_t *info,
				  QDP_DN_ColorMatrix *deriv[],
				  QOP_DN_GaugeField *g[],
				  QDP_DN_ColorMatrix **chain[],
				  int n, int doLastScale);

void QOP_DN_gauge_force_multi_qdp(QOP_info_t *info, QDP_DN_ColorMatrix *f[],
				  QOP_DN_GaugeField *g[],
				  QDP_DN_ColorMatrix **chain[], int n);

void QOP_DN_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_DN_GaugeField *gauge,
					QOP_D_Real *acts, QOP_D_Real *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_DN_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info,
					   QOP_DN_GaugeField *gauge,
					   QDP_DN_ColorMatrix *deriv[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_D_Real eps, int doLastScale);

void QOP_DN_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, 
					   QOP_DN_GaugeField *gauge, 
					   QDP_DN_ColorMatrix *force[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_D_Real eps);

void QOP_DN_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_DN_GaugeField *gauge, 
				       QOP_DN_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       QOP_D_Real eps);

void QOP_DN_projectU_qdp(QOP_info_t *info,
			 QDP_DN_ColorMatrix *pU,
			 QDP_DN_ColorMatrix *U);

void QOP_DN_projectU_deriv_qdp(QOP_info_t *info,
			       QDP_DN_ColorMatrix *f,
			       QDP_DN_ColorMatrix *pU,
			       QDP_DN_ColorMatrix *U,
			       QDP_DN_ColorMatrix *chain);

void QOP_DN_u3reunit(QOP_info_t *info, QDP_DN_ColorMatrix *U,
		     QDP_DN_ColorMatrix *V);

void QOP_DN_su3reunit(QOP_info_t *info, QDP_DN_ColorMatrix *U,
		      QDP_DN_ColorMatrix *Ur);

void QOP_DN_hisq_force_multi_reunit(QOP_info_t *info,
				    QDP_DN_ColorMatrix *gf[4],
				    QDP_DN_ColorMatrix *force_accum[4],
				    QDP_DN_ColorMatrix *force_accum_old[4]);

void QOP_DN_staples(QOP_info_t *info, int nout, int nin,
		    QDP_DN_ColorMatrix *out[], QDP_DN_ColorMatrix *in[],
		    int nstaples[], int *topdir[], int *sidedir[],
		    int *toplinknum[], int *sidelinknum[], QOP_D_Real *coef[]);

void QOP_DN_staples_deriv(QOP_info_t *info, int nout, int nin,
			  QDP_DN_ColorMatrix *deriv[],
			  QDP_DN_ColorMatrix *chain[],
			  QDP_DN_ColorMatrix *in[],
			  int nstaples[], int *topdir[], int *sidedir[],
			  int *toplinknum[], int *sidelinknum[],
			  QOP_D_Real *coef[]);

  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_DN_FermionLinksAsqtad *
  QOP_DN_asqtad_create_L_from_raw(int nc, QOP_D_Real *fatlinks[],
				  QOP_D_Real *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_DN_FermionLinksAsqtad *
  QOP_DN_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_DN_GaugeField *gauge);

void QOP_DN_asqtad_extract_L_to_raw(QOP_D_Real *fatlinks[],
				    QOP_D_Real *longlinks[],
				    QOP_DN_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_DN_asqtad_destroy_L(QOP_DN_FermionLinksAsqtad *field);

QOP_DN_FermionLinksAsqtad *
  QOP_DN_asqtad_convert_L_from_raw(int nc, QOP_D_Real *fatlinks[],
				   QOP_D_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_DN_asqtad_convert_L_to_raw(QOP_D_Real ***fatlinks,
				    QOP_D_Real ***longlinks,
				    QOP_DN_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_DN_asqtad_load_L_from_raw(QOP_DN_FermionLinksAsqtad *asqtad,
				   QOP_D_Real *fatlinks[],
				   QOP_D_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_DN_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_DN_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_DN_GaugeField *gauge);

void QOP_DN_asqtad_rephase_L(QOP_DN_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

QOP_DN_FermionLinksAsqtad *
  QOP_DN_asqtad_create_L_from_qdp(QDP_DN_ColorMatrix *fatlinks[],
				  QDP_DN_ColorMatrix *longlinks[]);

void QOP_DN_asqtad_extract_L_to_qdp(QDP_DN_ColorMatrix *fatlinks[],
				    QDP_DN_ColorMatrix *longlinks[],
				    QOP_DN_FermionLinksAsqtad *src);

QOP_DN_FermionLinksAsqtad *
  QOP_DN_asqtad_convert_L_from_qdp(QDP_DN_ColorMatrix *fatlinks[],
				   QDP_DN_ColorMatrix *longlinks[]);

void QOP_DN_asqtad_convert_L_to_qdp(QDP_DN_ColorMatrix ***fatlinks,
				    QDP_DN_ColorMatrix ***longlinks,
				    QOP_DN_FermionLinksAsqtad *src);

void QOP_DN_asqtad_load_L_from_qdp(QOP_DN_FermionLinksAsqtad *asqtad,
				   QDP_DN_ColorMatrix *fatlinks[],
				   QDP_DN_ColorMatrix *longlinks[]);

void QOP_DN_asqtad_rephase_field_L_qdp(QOP_DN_FermionLinksAsqtad *fla,
				       QDP_D_Complex *fatphase[],
				       QDP_D_Complex *longphase[]);


  /* inverter routines */

void QOP_DN_asqtad_dslash(QOP_info_t *info,
			  QOP_DN_FermionLinksAsqtad *asqtad,
			  QOP_D_Real mass,
			  QOP_DN_ColorVector *out,
			  QOP_DN_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_DN_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_DN_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_DN_ColorVector *out,
			      QOP_DN_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_DN_asqtad_diaginv(QOP_info_t *info,
			   QOP_DN_FermionLinksAsqtad *asqtad,
			   QOP_D_Real mass,
			   QOP_DN_ColorVector *out,
			   QOP_DN_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_DN_asqtad_invert(QOP_info_t *info,
			  QOP_DN_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real mass,
			  QOP_DN_ColorVector *out_pt,
			  QOP_DN_ColorVector *in_pt);

void QOP_DN_asqtad_invert_threaded(QOP_info_t *info,
				   QOP_DN_FermionLinksAsqtad *asqtad,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg,
				   QOP_D_Real mass,
				   QOP_DN_ColorVector *out_pt,
				   QOP_DN_ColorVector *in_pt,
				   int nthreads);

void QOP_DN_asqtad_invert_multi(QOP_info_t *info,
				QOP_DN_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *masses[],
				int nmass[],
				QOP_DN_ColorVector **out_pt[],
				QOP_DN_ColorVector *in_pt[],
				int nsrc);

void QOP_DN_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_DN_FermionLinksAsqtad *asqtad,
			      QOP_D_Real mass,
			      QDP_DN_ColorVector *out,
			      QDP_DN_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_DN_asqtad_dslash_dir_qdp(QOP_info_t *info,
				  QOP_DN_FermionLinksAsqtad *asqtad,
				  int dir, int fb,
				  double wtfat, double wtlong,
				  QDP_DN_ColorVector *out,
				  QDP_DN_ColorVector *in,
				  QOP_evenodd_t eo_out);

void QOP_DN_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_DN_FermionLinksAsqtad *asqtad,
			       QOP_D_Real mass,
			       QDP_DN_ColorVector *out,
			       QDP_DN_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_DN_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_DN_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_D_Real mass,
			      QDP_DN_ColorVector *out,
			      QDP_DN_ColorVector *in);

void QOP_DN_asqtad_invert_threaded_qdp(QOP_info_t *info,
				       QOP_DN_FermionLinksAsqtad *asqtad,
				       QOP_invert_arg_t *inv_arg,
				       QOP_resid_arg_t *res_arg,
				       QOP_D_Real mass,
				       QDP_DN_ColorVector *out,
				       QDP_DN_ColorVector *in,
				       int nthreads);

void QOP_DN_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_DN_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_D_Real *masses[],
				    int nmass[],
				    QDP_DN_ColorVector **out[],
				    QDP_DN_ColorVector *in[],
				    int nsrc);

void QOP_DN_asqtad_get_eigcg(QOP_DN_FermionLinksAsqtad *asqtad,
			     QDP_DN_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

  /* fermion force routines */

void QOP_DN_asqtad_deriv(QOP_info_t *info, QDP_DN_ColorMatrix *gauge[],
			 QDP_DN_ColorMatrix *force[],
			 QOP_asqtad_coeffs_t *coef,
			 QDP_DN_ColorMatrix *mid_fat[],
			 QDP_DN_ColorMatrix *mid_naik[]);

void QOP_DN_asqtad_force(QOP_info_t *info,
			 QOP_DN_GaugeField *gauge,
			 QOP_DN_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 QOP_D_Real eps,
			 QOP_DN_ColorVector *in_pt);

void QOP_DN_asqtad_force_multi(QOP_info_t *info,
			       QOP_DN_GaugeField *gauge,
			       QOP_DN_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       QOP_D_Real eps[],
			       QOP_DN_ColorVector *in_pt[],
			       int nsrc);

void QOP_DN_asqtad_force_multi_qdp(QOP_info_t *info,
				   QDP_DN_ColorMatrix *links[],
				   QDP_DN_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_D_Real eps[],
				   QDP_DN_ColorVector *in_pt[],
				   int nsrc);

void QOP_DN_asqtad_deriv_multi_qdp(QOP_info_t *info,
				   QDP_DN_ColorMatrix *links[],
				   QDP_DN_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_D_Real eps[],
				   QDP_DN_ColorVector *in_pt[],
				   int nsrc);

  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* single precision */

QOP_DN_FermionLinksHisq *
  QOP_DN_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_DN_GaugeField *gauge);

void QOP_DN_hisq_destroy_L(QOP_DN_FermionLinksHisq *field);

QOP_DN_FermionLinksAsqtad **
  QOP_DN_get_asqtad_links_from_hisq(QOP_DN_FermionLinksHisq *hl);
  
QOP_DN_FermionLinksAsqtad *
  QOP_DN_get_asqtad_deps_links_from_hisq(QOP_DN_FermionLinksHisq *hl);

  /* fermion force routines */

void QOP_DN_hisq_force_multi(QOP_info_t *info,
			     QOP_DN_FermionLinksHisq *flh,
			     QOP_DN_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     QOP_D_Real eps[],
			     QOP_DN_ColorVector *in_pt[],
			     int *n_orders_naik);

void QOP_DN_hisq_deriv_multi_qdp(QOP_info_t *info,
				 QOP_DN_FermionLinksHisq *flh,
				 QDP_DN_ColorMatrix *deriv[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_D_Real eps[],
				 QDP_DN_ColorVector *in_pt[],
				 int *n_orders_naik,
				 int doLastScale);

void QOP_DN_hisq_force_multi_qdp(QOP_info_t *info,
				 QOP_DN_FermionLinksHisq *flh,
				 QDP_DN_ColorMatrix *force[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_D_Real eps[],
				 QDP_DN_ColorVector *in_pt[],
				 int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_create_L_from_raw(int nc, QOP_D_Real *links[], QOP_D_Real *clov,
				  QOP_evenodd_t evenodd);

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_DN_GaugeField *gauge);

void QOP_DN_wilson_extract_L_to_raw(QOP_D_Real *links[], QOP_D_Real *clov,
				    QOP_DN_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_DN_wilson_destroy_L(QOP_DN_FermionLinksWilson *field);

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_convert_L_from_raw(int nc, QOP_D_Real *links[],
				   QOP_D_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_DN_wilson_convert_L_to_raw(QOP_D_Real ***links, QOP_D_Real **clov,
				    QOP_DN_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_DN_GaugeField *gauge);

QOP_DN_GaugeField *
  QOP_DN_wilson_convert_L_to_G(QOP_DN_FermionLinksWilson *links);

void QOP_DN_wilson_load_L_from_raw(QOP_DN_FermionLinksWilson *wilson,
				   QOP_D_Real *links[], QOP_D_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_DN_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_DN_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_DN_GaugeField *gauge);

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_create_L_from_qdp(QDP_DN_ColorMatrix *links[],
				  QDP_DN_DiracPropagator *clov);

void QOP_DN_wilson_extract_L_to_qdp(QDP_DN_ColorMatrix *links[],
				    QDP_DN_DiracPropagator *clov,
				    QOP_DN_FermionLinksWilson *src);

QOP_DN_FermionLinksWilson *
  QOP_DN_wilson_convert_L_from_qdp(QDP_DN_ColorMatrix *links[],
				   QDP_DN_DiracPropagator *clov);

void QOP_DN_wilson_convert_L_to_qdp(QDP_DN_ColorMatrix ***links,
				    QDP_DN_DiracPropagator **clov,
				    QOP_DN_FermionLinksWilson *src);

void QOP_DN_wilson_load_L_from_qdp(QOP_DN_FermionLinksWilson *wilson,
				   QDP_DN_ColorMatrix *links[],
				   QDP_DN_DiracPropagator *clov);


  /* inverter routines */

void QOP_DN_wilson_dslash(QOP_info_t *info,
			  QOP_DN_FermionLinksWilson *flw,
			  QOP_D_Real kappa,
			  int sign,
			  QOP_DN_DiracFermion *out,
			  QOP_DN_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_DN_wilson_diaginv(QOP_info_t *info,
			   QOP_DN_FermionLinksWilson *flw,
			   QOP_D_Real kappa,
			   QOP_DN_DiracFermion *out,
			   QOP_DN_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_DN_wilson_invert(QOP_info_t *info,
			  QOP_DN_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real kappa,
			  QOP_DN_DiracFermion *out_pt,
			  QOP_DN_DiracFermion *in_pt);

void QOP_DN_wilson_invert_multi(QOP_info_t *info,
				QOP_DN_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *kappas[],
				int nkappa[],
				QOP_DN_DiracFermion **out_pt[],
				QOP_DN_DiracFermion *in_pt[],
				int nsrc);

void QOP_DN_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_DN_FermionLinksWilson *flw,
			      QOP_D_Real kappa,
			      int sign,
			      QDP_DN_DiracFermion *out,
			      QDP_DN_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_DN_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_DN_FermionLinksWilson *flw,
			       QOP_D_Real kappa,
			       QDP_DN_DiracFermion *out,
			       QDP_DN_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_DN_wilson_invert_qdp(QOP_info_t *info,
			      QOP_DN_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_D_Real kappa,
			      QDP_DN_DiracFermion *out_pt,
			      QDP_DN_DiracFermion *in_pt);

void QOP_DN_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_DN_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_D_Real *kappas[],
				    int nkappa[],
				    QDP_DN_DiracFermion **out_pt[],
				    QDP_DN_DiracFermion *in_pt[],
				    int nsrc);

void QOP_DN_wilson_invert_ne_qdp(QOP_info_t *info,
				 QOP_DN_FermionLinksWilson *flw,
				 QOP_invert_arg_t *inv_arg,
				 QOP_resid_arg_t *res_arg,
				 QOP_D_Real kappa,
				 QDP_DN_DiracFermion *out,
				 QDP_DN_DiracFermion *in);

  /* fermion force routines */

void QOP_DN_wilson_force(QOP_info_t *info,
			 QOP_DN_GaugeField *gauge,
			 QOP_DN_Force *force,
			 QOP_wilson_coeffs_t *coeffs,
			 QOP_D_Real eps,
			 QOP_DN_DiracFermion *in_pt);

void QOP_DN_wilson_deriv_multi_qdp(QOP_info_t *info,
				   QOP_DN_FermionLinksWilson *flw,
				   QDP_DN_ColorMatrix *deriv[],
				   QOP_D_Real eps[],
				   QDP_DN_DiracFermion *x[],
				   QDP_DN_DiracFermion *y[],
				   int n);

void QOP_DN_wilson_force_multi(QOP_info_t *info,
			       QOP_DN_GaugeField *gauge,
			       QOP_DN_Force *force,
			       QOP_wilson_coeffs_t *coef,
			       QOP_D_Real eps[],
			       QOP_DN_DiracFermion *in_pt[],
			       int nsrc);

void QOP_DN_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
					QOP_DN_FermionLinksWilson *flw,
					QDP_DN_ColorMatrix *deriv[],
					QOP_D_Real kappa[],
					QOP_D_Real eps[],
					QDP_DN_DiracFermion *x[],
					QDP_DN_DiracFermion *y[],
					int n);

void QOP_DN_wilson_force_prec_multi_qdp(QOP_info_t *info,
					QOP_DN_FermionLinksWilson *flw,
					QDP_DN_ColorMatrix *force[],
					QOP_D_Real kappa[],
					QOP_D_Real eps[],
					QDP_DN_DiracFermion *x[],
					QDP_DN_DiracFermion *y[],
					int n);

  // new fermilab action IFLA -- added by bugra --------------- :

void QOP_DN_wilson_ifla_dslash(QOP_info_t *info,
			       QOP_DN_FermionLinksWilson *flw,
			       QOP_D_Real kappa,
			       int sign,
			       QOP_wilson_ifla_coeffs_t *coeffs,
			       QOP_DN_DiracFermion *out,
			       QOP_DN_DiracFermion *in,
			       QOP_evenodd_t eo_out,
			       QOP_evenodd_t eo_in);

void QOP_DN_wilson_ifla_dslash_qdp(QOP_info_t *info,
				   QOP_DN_FermionLinksWilson *flw,
				   QOP_D_Real kappa,
				   int sign,
				   QOP_wilson_ifla_coeffs_t *coeffs,
				   QDP_DN_DiracFermion *out,
				   QDP_DN_DiracFermion *in,
				   QOP_evenodd_t eo_out,
				   QOP_evenodd_t eo_in);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_DN_FermionLinksDW *
  QOP_DN_dw_create_L_from_raw(int nc, QOP_D_Real *links[],
			      QOP_evenodd_t evenodd);

QOP_DN_FermionLinksDW *
  QOP_DN_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_DN_GaugeField *gauge);

void QOP_DN_dw_extract_L_to_raw(QOP_D_Real *links[],
				QOP_DN_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_DN_dw_destroy_L(QOP_DN_FermionLinksDW *field);

QOP_DN_FermionLinksDW *
  QOP_DN_dw_convert_L_from_raw(int nc, QOP_D_Real *links[],
			       QOP_evenodd_t evenodd);

void QOP_DN_dw_convert_L_to_raw(QOP_D_Real ***links,
				QOP_DN_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_DN_FermionLinksDW *
  QOP_DN_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_DN_GaugeField *gauge);

QOP_DN_GaugeField *
  QOP_DN_dw_convert_L_to_G(QOP_DN_FermionLinksDW *links);

void QOP_DN_dw_load_L_from_raw(QOP_DN_FermionLinksDW *dw,
			       QOP_D_Real *links[], QOP_evenodd_t evenodd);

void QOP_DN_dw_load_L_from_G(QOP_info_t *info,
			     QOP_DN_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_DN_GaugeField *gauge);

QOP_DN_FermionLinksDW *
  QOP_DN_dw_create_L_from_qdp(QDP_DN_ColorMatrix *links[]);

void QOP_DN_dw_extract_L_to_qdp(QDP_DN_ColorMatrix *links[],
				QOP_DN_FermionLinksDW *src);

QOP_DN_FermionLinksDW *
  QOP_DN_dw_convert_L_from_qdp(QDP_DN_ColorMatrix *links[]);

void QOP_DN_dw_convert_L_to_qdp(QDP_DN_ColorMatrix ***links,
				QOP_DN_FermionLinksDW *src);

void QOP_DN_dw_load_L_from_qdp(QOP_DN_FermionLinksDW *dw,
			       QDP_DN_ColorMatrix *links[]);

  /* inverter routines */

void QOP_DN_dw_dslash(QOP_info_t *info,
		      QOP_DN_FermionLinksDW *links,
		      QOP_D_Real M5,
		      QOP_D_Real m,
		      int sign,
		      QOP_DN_DiracFermion *out_pt[],
		      QOP_DN_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_DN_dw_dslash2(QOP_info_t *info,
		       QOP_DN_FermionLinksDW *links,
		       QOP_D_Real M5,
		       QOP_D_Real m,
		       QOP_DN_DiracFermion *out_pt[],
		       QOP_DN_DiracFermion *in_pt[],
		       int Ls,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in);

void QOP_DN_dw_invert(QOP_info_t *info,
		      QOP_DN_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QOP_D_Real M5,
		      QOP_D_Real m,
		      QOP_DN_DiracFermion *out_pt[],
		      QOP_DN_DiracFermion *in_pt[],
		      int Ls);

void QOP_DN_dw_invert_multi(QOP_info_t *info,
			    QOP_DN_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_D_Real *M5[],
			    QOP_D_Real *m[],
			    int nmass[],
			    QOP_DN_DiracFermion ***out_pt[],
			    QOP_DN_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_DN_dw_dslash_qdp(QOP_info_t *info,
			  QOP_DN_FermionLinksDW *links,
			  QOP_D_Real M5, 
			  QOP_D_Real m,
			  int sign,
			  QDP_DN_DiracFermion *out_pt[],
			  QDP_DN_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_DN_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_DN_FermionLinksDW *links,
			   QOP_D_Real M5,
			   QOP_D_Real m,
			   QDP_DN_DiracFermion *out_pt[],
			   QDP_DN_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);

void QOP_DN_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_DN_FermionLinksDW *fldw,
			   QOP_D_Real M5,
			   QOP_D_Real m,
			   QDP_DN_DiracFermion **out,
			   QDP_DN_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_DN_dw_invert_qdp(QOP_info_t *info,
			  QOP_DN_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_D_Real M5,
			  QOP_D_Real m,
			  QDP_DN_DiracFermion *out[],
			  QDP_DN_DiracFermion *in[],
			  int Ls);

void QOP_DN_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_DN_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_D_Real *M5[],
				QOP_D_Real *m[],
				int nmass[],
				QDP_DN_DiracFermion ***out[],
				QDP_DN_DiracFermion **in[],
				int nsrc,
				int Ls);

 /* fermion force routines */

void QOP_DN_dw_force(QOP_info_t *info,
		     QOP_DN_GaugeField *gauge,
		     QOP_DN_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     QOP_D_Real eps,
		     QOP_DN_DiracFermion *in_pt);

void QOP_DN_dw_force_multi(QOP_info_t *info,
			   QOP_DN_GaugeField *gauge,
			   QOP_DN_Force *force,
			   QOP_dw_coeffs_t *coef,
			   QOP_D_Real eps[],
			   QOP_DN_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == 'D'
#  if QOP_Colors == 'N'
#    include <qop_dn_generic.h>
#  endif
#  include <qop_dn_precision_generic.h>
#endif
#if QOP_Colors == 'N'
#  include <qop_dn_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_DN_H */
