// DO NOT EDIT
// generated from qop_pc.h
#ifndef _QOP_F1_H
#define _QOP_F1_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QOP_F1_ColorVector_struct   QOP_F1_ColorVector;
typedef struct QOP_F1_DiracFermion_struct  QOP_F1_DiracFermion;
typedef struct QOP_F1_GaugeField_struct    QOP_F1_GaugeField;
typedef struct QOP_F1_Force_struct         QOP_F1_Force;

typedef struct QOP_F1_FermionLinksAsqtad_struct  QOP_F1_FermionLinksAsqtad;
typedef struct QOP_F1_FermionLinksHisq_struct    QOP_F1_FermionLinksHisq;
typedef struct QOP_F1_FermionLinksWilson_struct  QOP_F1_FermionLinksWilson;
typedef struct QOP_F1_FermionLinksDW_struct      QOP_F1_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

#define QOP_F1_qla_type_V QLA_F1_ColorVector
#define QOP_F1_qla_type_D QLA_F1_DiracFermion
#define QOP_F1_qla_type_M QLA_F1_ColorMatrix
#define QOP_F1_raw_size(T) (QDP_sites_on_node*sizeof(QOP_F1_qla_type_##T))
#define QOP_F1_raw_size_V(evenodd) QOP_F1_raw_size(V)
#define QOP_F1_raw_size_D(evenodd) QOP_F1_raw_size(D)
#define QOP_F1_raw_size_G(evenodd) QOP_F1_raw_size(M)
#define QOP_F1_raw_size_F(evenodd) QOP_F1_raw_size(M)

#define QOP_F1_elem(T, raw, i, ...) QLA_F1_elem_##T(((QOP_F1_qla_type_##T *)raw)[i], __ARGV__)
#define QOP_F1_set(T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_F1_elem(T, raw, i, __ARGV__), re, im)
#define QOP_F1_raw_set_V(raw, evenodd, i, ic, re, im) QOP_F1_set(V, raw, i, re, im, ic)
#define QOP_F1_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_F1_set(D, raw, i, re, im, ic, is)
#define QOP_F1_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_F1_set(M, raw, i, re, im, ic, jc)
#define QOP_F1_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_F1_set(M, raw, i, re, im, ic, jc)
#define QOP_F1_get(T, re, im, raw, i, ...) {			\
    QLA_F_Complex _c = QOP_F1_elem(T, raw, i, __ARGV__);	\
    re = QLA_real(_c); im = QLA_imag(_c);			\
  }
#define QOP_F1_raw_get_V(re, im, raw, evenodd, i, ic)     QOP_F1_set(V, re, im, raw, i, ic)
#define QOP_F1_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_F1_set(D, re, im, raw, i, ic, is)
#define QOP_F1_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_F1_set(M, re, im, raw, i, ic, jc)
#define QOP_F1_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_F1_set(M, re, im, raw, i, ic, jc)

/* create a QOP field with a copy of the raw source field */
QOP_F1_ColorVector  *QOP_F1_create_V_from_raw( QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_F1_DiracFermion *QOP_F1_create_D_from_raw( QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_F1_GaugeField   *QOP_F1_create_G_from_raw( QOP_F_Real *links[], QOP_evenodd_t evenodd);
QOP_F1_Force        *QOP_F1_create_F_from_raw( QOP_F_Real *force[], QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_F1_extract_V_to_raw(QOP_F_Real *dest, QOP_F1_ColorVector *src, QOP_evenodd_t evenodd);
void QOP_F1_extract_D_to_raw(QOP_F_Real *dest, QOP_F1_DiracFermion *src, QOP_evenodd_t evenodd);
void QOP_F1_extract_G_to_raw(QOP_F_Real *dest[], QOP_F1_GaugeField *src, QOP_evenodd_t evenodd);
void QOP_F1_extract_F_to_raw(QOP_F_Real *dest[], QOP_F1_Force *src, QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_F1_destroy_V(QOP_F1_ColorVector *field);
void QOP_F1_destroy_D(QOP_F1_DiracFermion *field);
void QOP_F1_destroy_G(QOP_F1_GaugeField *field);
void QOP_F1_destroy_F(QOP_F1_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_F1_ColorVector  *QOP_F1_convert_V_from_raw( QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_F1_DiracFermion *QOP_F1_convert_D_from_raw( QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_F1_GaugeField   *QOP_F1_convert_G_from_raw( QOP_F_Real *links[], QOP_evenodd_t evenodd);
QOP_F1_Force        *QOP_F1_convert_F_from_raw( QOP_F_Real *force[], QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
QOP_F_Real  *QOP_F1_convert_V_to_raw(QOP_F1_ColorVector *src, QOP_evenodd_t evenodd);
QOP_F_Real  *QOP_F1_convert_D_to_raw(QOP_F1_DiracFermion *src, QOP_evenodd_t evenodd);
QOP_F_Real **QOP_F1_convert_G_to_raw(QOP_F1_GaugeField *src, QOP_evenodd_t evenodd);
QOP_F_Real **QOP_F1_convert_F_to_raw(QOP_F1_Force *src, QOP_evenodd_t evenodd);

QOP_F1_ColorVector  *QOP_F1_create_V_from_qdp(QDP_F1_ColorVector *src);
QOP_F1_DiracFermion *QOP_F1_create_D_from_qdp(QDP_F1_DiracFermion *src);
QOP_F1_GaugeField   *QOP_F1_create_G_from_qdp(QDP_F1_ColorMatrix *src[]);
QOP_F1_Force        *QOP_F1_create_F_from_qdp(QDP_F1_ColorMatrix *src[]);

void QOP_F1_extract_V_to_qdp(QDP_F1_ColorVector *d, QOP_F1_ColorVector *src);
void QOP_F1_extract_D_to_qdp(QDP_F1_DiracFermion *d, QOP_F1_DiracFermion *src);
void QOP_F1_extract_G_to_qdp(QDP_F1_ColorMatrix *d[], QOP_F1_GaugeField *src);
void QOP_F1_extract_F_to_qdp(QDP_F1_ColorMatrix *d[], QOP_F1_Force *src);

QOP_F1_ColorVector  *QOP_F1_convert_V_from_qdp(QDP_F1_ColorVector *src);
QOP_F1_DiracFermion *QOP_F1_convert_D_from_qdp(QDP_F1_DiracFermion *src);
QOP_F1_GaugeField   *QOP_F1_convert_G_from_qdp(QDP_F1_ColorMatrix *src[]);
QOP_F1_Force        *QOP_F1_convert_F_from_qdp(QDP_F1_ColorMatrix *src[]);

QDP_F1_ColorVector   *QOP_F1_convert_V_to_qdp(QOP_F1_ColorVector *src);
QDP_F1_DiracFermion  *QOP_F1_convert_D_to_qdp(QOP_F1_DiracFermion *src);
QDP_F1_ColorMatrix  **QOP_F1_convert_G_to_qdp(QOP_F1_GaugeField *src);
QDP_F1_ColorMatrix  **QOP_F1_convert_F_to_qdp(QOP_F1_Force *src);


  /********************/
  /*  Gauge routines  */
  /********************/

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */
void QOP_F1_rephase_G(QOP_F1_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_F1_rephase_G_qdp(QDP_F1_ColorMatrix *links[],
			  int *r0,
			  QOP_bc_t *bc,
			  QOP_staggered_sign_t *sign);

void QOP_F1_smear_fat7l_qdp(QOP_info_t *info, QDP_F1_ColorMatrix *sg[],
			    QDP_F1_ColorMatrix *g[],
			    QOP_asqtad_coeffs_t *coeffs);

void QOP_F1_gauge_deriv_multi_qdp(QOP_info_t *info,
				  QDP_F1_ColorMatrix *deriv[],
				  QOP_F1_GaugeField *g[],
				  QDP_F1_ColorMatrix **chain[],
				  int n, int doLastScale);

void QOP_F1_gauge_force_multi_qdp(QOP_info_t *info, QDP_F1_ColorMatrix *f[],
				  QOP_F1_GaugeField *g[],
				  QDP_F1_ColorMatrix **chain[], int n);

void QOP_F1_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_F1_GaugeField *gauge,
					QOP_F_Real *acts, QOP_F_Real *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_F1_symanzik_1loop_gauge_action_qdp(QOP_info_t *info,
					    QDP_F1_ColorMatrix *links[],
					    QOP_F_Real *acts, QOP_F_Real *actt,
					    QOP_gauge_coeffs_t *coeffs);

void QOP_F1_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_F1_GaugeField *gauge, 
				       QOP_F1_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       QOP_F_Real eps);

void QOP_F1_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, 
					   QDP_F1_ColorMatrix *links[],
					   QDP_F1_ColorMatrix *force[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_F_Real eps);

void QOP_F1_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info,
					   QDP_F1_ColorMatrix *links[],
					   QDP_F1_ColorMatrix *deriv[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_F_Real eps, int doLastScale);

void QOP_F1_symanzik_1loop_gauge_heatbath_qdp(QOP_info_t *info,
					      QDP_F1_ColorMatrix *links[],
					      QLA_F_Real beta,
					      QOP_gauge_coeffs_t *coeffs,
					      QDP_RandomState *rs0,
					      int nup, int nhb, int nover);

void QOP_F1_symanzik_1loop_gauge_staple_qdp(QOP_info_t *info,
					    QDP_F1_ColorMatrix *links[],
					    QDP_F1_ColorMatrix *staple,
					    int mu,
					    QOP_gauge_coeffs_t *coeffs,
					    QDP_Subset subs[], int subi);

void QOP_F1_projectU_qdp(QOP_info_t *info,
			 QDP_F1_ColorMatrix *pU,
			 QDP_F1_ColorMatrix *U);

void QOP_F1_projectU_deriv_qdp(QOP_info_t *info,
			       QDP_F1_ColorMatrix *f,
			       QDP_F1_ColorMatrix *pU,
			       QDP_F1_ColorMatrix *U,
			       QDP_F1_ColorMatrix *chain);

void QOP_F1_u3reunit(QOP_info_t *info, QDP_F1_ColorMatrix *U,
		     QDP_F1_ColorMatrix *V);

void QOP_F1_su3reunit(QOP_info_t *info, QDP_F1_ColorMatrix *U,
		      QDP_F1_ColorMatrix *Ur);

void QOP_F1_hisq_force_multi_reunit(QOP_info_t *info,
				    QDP_F1_ColorMatrix *gf[4],
				    QDP_F1_ColorMatrix *force_accum[4],
				    QDP_F1_ColorMatrix *force_accum_old[4]);

void QOP_F1_staples(QOP_info_t *info, int nout, int nin,
		    QDP_F1_ColorMatrix *out[], QDP_F1_ColorMatrix *in[],
		    int nstaples[], int *topdir[], int *sidedir[],
		    int *toplinknum[], int *sidelinknum[], QOP_F_Real *coef[]);

void QOP_F1_staples_deriv(QOP_info_t *info, int nout, int nin,
			  QDP_F1_ColorMatrix *deriv[],
			  QDP_F1_ColorMatrix *chain[],
			  QDP_F1_ColorMatrix *in[],
			  int nstaples[], int *topdir[], int *sidedir[],
			  int *toplinknum[], int *sidelinknum[],
			  QOP_F_Real *coef[]);

  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_F1_FermionLinksAsqtad *
  QOP_F1_asqtad_create_L_from_raw( QOP_F_Real *fatlinks[],
				  QOP_F_Real *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_F1_FermionLinksAsqtad *
  QOP_F1_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_F1_GaugeField *gauge);

void QOP_F1_asqtad_extract_L_to_raw(QOP_F_Real *fatlinks[],
				    QOP_F_Real *longlinks[],
				    QOP_F1_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_F1_asqtad_destroy_L(QOP_F1_FermionLinksAsqtad *field);

QOP_F1_FermionLinksAsqtad *
  QOP_F1_asqtad_convert_L_from_raw( QOP_F_Real *fatlinks[],
				   QOP_F_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_F1_asqtad_convert_L_to_raw(QOP_F_Real ***fatlinks,
				    QOP_F_Real ***longlinks,
				    QOP_F1_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_F1_asqtad_load_L_from_raw(QOP_F1_FermionLinksAsqtad *asqtad,
				   QOP_F_Real *fatlinks[],
				   QOP_F_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_F1_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_F1_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_F1_GaugeField *gauge);

void QOP_F1_asqtad_rephase_L(QOP_F1_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

QOP_F1_FermionLinksAsqtad *
  QOP_F1_asqtad_create_L_from_qdp(QDP_F1_ColorMatrix *fatlinks[],
				  QDP_F1_ColorMatrix *longlinks[]);

void QOP_F1_asqtad_extract_L_to_qdp(QDP_F1_ColorMatrix *fatlinks[],
				    QDP_F1_ColorMatrix *longlinks[],
				    QOP_F1_FermionLinksAsqtad *src);

QOP_F1_FermionLinksAsqtad *
  QOP_F1_asqtad_convert_L_from_qdp(QDP_F1_ColorMatrix *fatlinks[],
				   QDP_F1_ColorMatrix *longlinks[]);

void QOP_F1_asqtad_convert_L_to_qdp(QDP_F1_ColorMatrix ***fatlinks,
				    QDP_F1_ColorMatrix ***longlinks,
				    QOP_F1_FermionLinksAsqtad *src);

void QOP_F1_asqtad_load_L_from_qdp(QOP_F1_FermionLinksAsqtad *asqtad,
				   QDP_F1_ColorMatrix *fatlinks[],
				   QDP_F1_ColorMatrix *longlinks[]);

void QOP_F1_asqtad_rephase_field_L_qdp(QOP_F1_FermionLinksAsqtad *fla,
				       QDP_F_Complex *fatphase[],
				       QDP_F_Complex *longphase[]);


  /* inverter routines */

void QOP_F1_asqtad_dslash(QOP_info_t *info,
			  QOP_F1_FermionLinksAsqtad *asqtad,
			  QOP_F_Real mass,
			  QOP_F1_ColorVector *out,
			  QOP_F1_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_F1_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_F1_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_F1_ColorVector *out,
			      QOP_F1_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_F1_asqtad_diaginv(QOP_info_t *info,
			   QOP_F1_FermionLinksAsqtad *asqtad,
			   QOP_F_Real mass,
			   QOP_F1_ColorVector *out,
			   QOP_F1_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_F1_asqtad_invert(QOP_info_t *info,
			  QOP_F1_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real mass,
			  QOP_F1_ColorVector *out_pt,
			  QOP_F1_ColorVector *in_pt);

void QOP_F1_asqtad_invert_threaded(QOP_info_t *info,
				   QOP_F1_FermionLinksAsqtad *asqtad,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg,
				   QOP_F_Real mass,
				   QOP_F1_ColorVector *out_pt,
				   QOP_F1_ColorVector *in_pt,
				   int nthreads);

void QOP_F1_asqtad_invert_multi(QOP_info_t *info,
				QOP_F1_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *masses[],
				int nmass[],
				QOP_F1_ColorVector **out_pt[],
				QOP_F1_ColorVector *in_pt[],
				int nsrc);

void QOP_F1_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_F1_FermionLinksAsqtad *asqtad,
			      QOP_F_Real mass,
			      QDP_F1_ColorVector *out,
			      QDP_F1_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_F1_asqtad_dslash_dir_qdp(QOP_info_t *info,
				  QOP_F1_FermionLinksAsqtad *asqtad,
				  int dir, int fb,
				  double wtfat, double wtlong,
				  QDP_F1_ColorVector *out,
				  QDP_F1_ColorVector *in,
				  QOP_evenodd_t eo_out);

void QOP_F1_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_F1_FermionLinksAsqtad *asqtad,
			       QOP_F_Real mass,
			       QDP_F1_ColorVector *out,
			       QDP_F1_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_F1_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_F1_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_F_Real mass,
			      QDP_F1_ColorVector *out,
			      QDP_F1_ColorVector *in);

void QOP_F1_asqtad_invert_threaded_qdp(QOP_info_t *info,
				       QOP_F1_FermionLinksAsqtad *asqtad,
				       QOP_invert_arg_t *inv_arg,
				       QOP_resid_arg_t *res_arg,
				       QOP_F_Real mass,
				       QDP_F1_ColorVector *out,
				       QDP_F1_ColorVector *in,
				       int nthreads);

void QOP_F1_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_F1_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_F_Real *masses[],
				    int nmass[],
				    QDP_F1_ColorVector **out[],
				    QDP_F1_ColorVector *in[],
				    int nsrc);

void QOP_F1_asqtad_get_eigcg(QOP_F1_FermionLinksAsqtad *asqtad,
			     QDP_F1_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

  /* fermion force routines */

void QOP_F1_asqtad_deriv(QOP_info_t *info, QDP_F1_ColorMatrix *gauge[],
			 QDP_F1_ColorMatrix *force[],
			 QOP_asqtad_coeffs_t *coef,
			 QDP_F1_ColorMatrix *mid_fat[],
			 QDP_F1_ColorMatrix *mid_naik[]);

void QOP_F1_asqtad_force(QOP_info_t *info,
			 QOP_F1_GaugeField *gauge,
			 QOP_F1_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 QOP_F_Real eps,
			 QOP_F1_ColorVector *in_pt);

void QOP_F1_asqtad_force_multi(QOP_info_t *info,
			       QOP_F1_GaugeField *gauge,
			       QOP_F1_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       QOP_F_Real eps[],
			       QOP_F1_ColorVector *in_pt[],
			       int nsrc);

void QOP_F1_asqtad_force_multi_qdp(QOP_info_t *info,
				   QDP_F1_ColorMatrix *links[],
				   QDP_F1_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_F_Real eps[],
				   QDP_F1_ColorVector *in_pt[],
				   int nsrc);

void QOP_F1_asqtad_deriv_multi_qdp(QOP_info_t *info,
				   QDP_F1_ColorMatrix *links[],
				   QDP_F1_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_F_Real eps[],
				   QDP_F1_ColorVector *in_pt[],
				   int nsrc);

  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* single precision */

QOP_F1_FermionLinksHisq *
  QOP_F1_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_F1_GaugeField *gauge);

void QOP_F1_hisq_destroy_L(QOP_F1_FermionLinksHisq *field);

QOP_F1_FermionLinksAsqtad **
  QOP_F1_get_asqtad_links_from_hisq(QOP_F1_FermionLinksHisq *hl);
  
QOP_F1_FermionLinksAsqtad *
  QOP_F1_get_asqtad_deps_links_from_hisq(QOP_F1_FermionLinksHisq *hl);

  /* fermion force routines */

void QOP_F1_hisq_force_multi(QOP_info_t *info,
			     QOP_F1_FermionLinksHisq *flh,
			     QOP_F1_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     QOP_F_Real eps[],
			     QOP_F1_ColorVector *in_pt[],
			     int *n_orders_naik);

void QOP_F1_hisq_deriv_multi_qdp(QOP_info_t *info,
				 QOP_F1_FermionLinksHisq *flh,
				 QDP_F1_ColorMatrix *deriv[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_F_Real eps[],
				 QDP_F1_ColorVector *in_pt[],
				 int *n_orders_naik,
				 int doLastScale);

void QOP_F1_hisq_force_multi_qdp(QOP_info_t *info,
				 QOP_F1_FermionLinksHisq *flh,
				 QDP_F1_ColorMatrix *force[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_F_Real eps[],
				 QDP_F1_ColorVector *in_pt[],
				 int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_create_L_from_raw( QOP_F_Real *links[], QOP_F_Real *clov,
				  QOP_evenodd_t evenodd);

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_F1_GaugeField *gauge);

void QOP_F1_wilson_extract_L_to_raw(QOP_F_Real *links[], QOP_F_Real *clov,
				    QOP_F1_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_F1_wilson_destroy_L(QOP_F1_FermionLinksWilson *field);

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_convert_L_from_raw( QOP_F_Real *links[],
				   QOP_F_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_F1_wilson_convert_L_to_raw(QOP_F_Real ***links, QOP_F_Real **clov,
				    QOP_F1_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_F1_GaugeField *gauge);

QOP_F1_GaugeField *
  QOP_F1_wilson_convert_L_to_G(QOP_F1_FermionLinksWilson *links);

void QOP_F1_wilson_load_L_from_raw(QOP_F1_FermionLinksWilson *wilson,
				   QOP_F_Real *links[], QOP_F_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_F1_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_F1_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_F1_GaugeField *gauge);

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_create_L_from_qdp(QDP_F1_ColorMatrix *links[],
				  QDP_F1_DiracPropagator *clov);

void QOP_F1_wilson_extract_L_to_qdp(QDP_F1_ColorMatrix *links[],
				    QDP_F1_DiracPropagator *clov,
				    QOP_F1_FermionLinksWilson *src);

QOP_F1_FermionLinksWilson *
  QOP_F1_wilson_convert_L_from_qdp(QDP_F1_ColorMatrix *links[],
				   QDP_F1_DiracPropagator *clov);

void QOP_F1_wilson_convert_L_to_qdp(QDP_F1_ColorMatrix ***links,
				    QDP_F1_DiracPropagator **clov,
				    QOP_F1_FermionLinksWilson *src);

void QOP_F1_wilson_load_L_from_qdp(QOP_F1_FermionLinksWilson *wilson,
				   QDP_F1_ColorMatrix *links[],
				   QDP_F1_DiracPropagator *clov);


  /* inverter routines */

void QOP_F1_wilson_dslash(QOP_info_t *info,
			  QOP_F1_FermionLinksWilson *flw,
			  QOP_F_Real kappa,
			  int sign,
			  QOP_F1_DiracFermion *out,
			  QOP_F1_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_F1_wilson_diaginv(QOP_info_t *info,
			   QOP_F1_FermionLinksWilson *flw,
			   QOP_F_Real kappa,
			   QOP_F1_DiracFermion *out,
			   QOP_F1_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_F1_wilson_invert(QOP_info_t *info,
			  QOP_F1_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real kappa,
			  QOP_F1_DiracFermion *out_pt,
			  QOP_F1_DiracFermion *in_pt);

void QOP_F1_wilson_invert_multi(QOP_info_t *info,
				QOP_F1_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *kappas[],
				int nkappa[],
				QOP_F1_DiracFermion **out_pt[],
				QOP_F1_DiracFermion *in_pt[],
				int nsrc);

void QOP_F1_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_F1_FermionLinksWilson *flw,
			      QOP_F_Real kappa,
			      int sign,
			      QDP_F1_DiracFermion *out,
			      QDP_F1_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_F1_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_F1_FermionLinksWilson *flw,
			       QOP_F_Real kappa,
			       QDP_F1_DiracFermion *out,
			       QDP_F1_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_F1_wilson_invert_qdp(QOP_info_t *info,
			      QOP_F1_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_F_Real kappa,
			      QDP_F1_DiracFermion *out_pt,
			      QDP_F1_DiracFermion *in_pt);

void QOP_F1_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_F1_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_F_Real *kappas[],
				    int nkappa[],
				    QDP_F1_DiracFermion **out_pt[],
				    QDP_F1_DiracFermion *in_pt[],
				    int nsrc);

void QOP_F1_wilson_invert_ne_qdp(QOP_info_t *info,
				 QOP_F1_FermionLinksWilson *flw,
				 QOP_invert_arg_t *inv_arg,
				 QOP_resid_arg_t *res_arg,
				 QOP_F_Real kappa,
				 QDP_F1_DiracFermion *out,
				 QDP_F1_DiracFermion *in);

  /* fermion force routines */

void QOP_F1_wilson_deriv_multi_qdp(QOP_info_t *info,
				   QOP_F1_FermionLinksWilson *flw,
				   QDP_F1_ColorMatrix *deriv[],
				   QOP_F_Real eps[],
				   QDP_F1_DiracFermion *x[],
				   QDP_F1_DiracFermion *y[],
				   int n);

void QOP_F1_wilson_force_multi_qdp(QOP_info_t *info,
				   QOP_F1_FermionLinksWilson *flw,
				   QDP_F1_ColorMatrix *force[],
				   QOP_F_Real eps[],
				   QDP_F1_DiracFermion *x[],
				   QDP_F1_DiracFermion *y[],
				   int n);

void QOP_F1_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
					QOP_F1_FermionLinksWilson *flw,
					QDP_F1_ColorMatrix *deriv[],
					QOP_F_Real kappa[],
					QOP_F_Real eps[],
					QDP_F1_DiracFermion *x[],
					QDP_F1_DiracFermion *y[],
					int n);

void QOP_F1_wilson_force_prec_multi_qdp(QOP_info_t *info,
					QOP_F1_FermionLinksWilson *flw,
					QDP_F1_ColorMatrix *force[],
					QOP_F_Real kappa[],
					QOP_F_Real eps[],
					QDP_F1_DiracFermion *x[],
					QDP_F1_DiracFermion *y[],
					int n);

  // new fermilab action IFLA -- added by bugra --------------- :

void QOP_F1_wilson_ifla_dslash(QOP_info_t *info,
			       QOP_F1_FermionLinksWilson *flw,
			       QOP_F_Real kappa,
			       int sign,
			       QOP_wilson_ifla_coeffs_t *coeffs,
			       QOP_F1_DiracFermion *out,
			       QOP_F1_DiracFermion *in,
			       QOP_evenodd_t eo_out,
			       QOP_evenodd_t eo_in);

void QOP_F1_wilson_ifla_dslash_qdp(QOP_info_t *info,
				   QOP_F1_FermionLinksWilson *flw,
				   QOP_F_Real kappa,
				   int sign,
				   QOP_wilson_ifla_coeffs_t *coeffs,
				   QDP_F1_DiracFermion *out,
				   QDP_F1_DiracFermion *in,
				   QOP_evenodd_t eo_out,
				   QOP_evenodd_t eo_in);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_F1_FermionLinksDW *
  QOP_F1_dw_create_L_from_raw( QOP_F_Real *links[],
			      QOP_evenodd_t evenodd);

QOP_F1_FermionLinksDW *
  QOP_F1_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_F1_GaugeField *gauge);

void QOP_F1_dw_extract_L_to_raw(QOP_F_Real *links[],
				QOP_F1_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_F1_dw_destroy_L(QOP_F1_FermionLinksDW *field);

QOP_F1_FermionLinksDW *
  QOP_F1_dw_convert_L_from_raw( QOP_F_Real *links[],
			       QOP_evenodd_t evenodd);

void QOP_F1_dw_convert_L_to_raw(QOP_F_Real ***links,
				QOP_F1_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_F1_FermionLinksDW *
  QOP_F1_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_F1_GaugeField *gauge);

QOP_F1_GaugeField *
  QOP_F1_dw_convert_L_to_G(QOP_F1_FermionLinksDW *links);

void QOP_F1_dw_load_L_from_raw(QOP_F1_FermionLinksDW *dw,
			       QOP_F_Real *links[], QOP_evenodd_t evenodd);

void QOP_F1_dw_load_L_from_G(QOP_info_t *info,
			     QOP_F1_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_F1_GaugeField *gauge);

QOP_F1_FermionLinksDW *
  QOP_F1_dw_create_L_from_qdp(QDP_F1_ColorMatrix *links[]);

void QOP_F1_dw_extract_L_to_qdp(QDP_F1_ColorMatrix *links[],
				QOP_F1_FermionLinksDW *src);

QOP_F1_FermionLinksDW *
  QOP_F1_dw_convert_L_from_qdp(QDP_F1_ColorMatrix *links[]);

void QOP_F1_dw_convert_L_to_qdp(QDP_F1_ColorMatrix ***links,
				QOP_F1_FermionLinksDW *src);

void QOP_F1_dw_load_L_from_qdp(QOP_F1_FermionLinksDW *dw,
			       QDP_F1_ColorMatrix *links[]);

  /* inverter routines */

void QOP_F1_dw_dslash(QOP_info_t *info,
		      QOP_F1_FermionLinksDW *links,
		      QOP_F_Real M5,
		      QOP_F_Real m,
		      int sign,
		      QOP_F1_DiracFermion *out_pt[],
		      QOP_F1_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_F1_dw_dslash2(QOP_info_t *info,
		       QOP_F1_FermionLinksDW *links,
		       QOP_F_Real M5,
		       QOP_F_Real m,
		       QOP_F1_DiracFermion *out_pt[],
		       QOP_F1_DiracFermion *in_pt[],
		       int Ls,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in);

void QOP_F1_dw_invert(QOP_info_t *info,
		      QOP_F1_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QOP_F_Real M5,
		      QOP_F_Real m,
		      QOP_F1_DiracFermion *out_pt[],
		      QOP_F1_DiracFermion *in_pt[],
		      int Ls);

void QOP_F1_dw_invert_multi(QOP_info_t *info,
			    QOP_F1_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_F_Real *M5[],
			    QOP_F_Real *m[],
			    int nmass[],
			    QOP_F1_DiracFermion ***out_pt[],
			    QOP_F1_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_F1_dw_dslash_qdp(QOP_info_t *info,
			  QOP_F1_FermionLinksDW *links,
			  QOP_F_Real M5, 
			  QOP_F_Real m,
			  int sign,
			  QDP_F1_DiracFermion *out_pt[],
			  QDP_F1_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_F1_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_F1_FermionLinksDW *links,
			   QOP_F_Real M5,
			   QOP_F_Real m,
			   QDP_F1_DiracFermion *out_pt[],
			   QDP_F1_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);

void QOP_F1_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_F1_FermionLinksDW *fldw,
			   QOP_F_Real M5,
			   QOP_F_Real m,
			   QDP_F1_DiracFermion **out,
			   QDP_F1_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_F1_dw_invert_qdp(QOP_info_t *info,
			  QOP_F1_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real M5,
			  QOP_F_Real m,
			  QDP_F1_DiracFermion *out[],
			  QDP_F1_DiracFermion *in[],
			  int Ls);

void QOP_F1_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_F1_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *M5[],
				QOP_F_Real *m[],
				int nmass[],
				QDP_F1_DiracFermion ***out[],
				QDP_F1_DiracFermion **in[],
				int nsrc,
				int Ls);

 /* fermion force routines */

void QOP_F1_dw_force(QOP_info_t *info,
		     QOP_F1_GaugeField *gauge,
		     QOP_F1_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     QOP_F_Real eps,
		     QOP_F1_DiracFermion *in_pt);

void QOP_F1_dw_force_multi(QOP_info_t *info,
			   QOP_F1_GaugeField *gauge,
			   QOP_F1_Force *force,
			   QOP_dw_coeffs_t *coef,
			   QOP_F_Real eps[],
			   QOP_F1_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == 'F'
#  if QOP_Colors == 1
#    include <qop_f1_generic.h>
#  endif
#  include <qop_f1_precision_generic.h>
#endif
#if QOP_Colors == 1
#  include <qop_f1_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_F1_H */
