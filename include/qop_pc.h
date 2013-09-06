#ifndef _QOP_PC_H
#define _QOP_PC_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QOP_IPC_ColorVector_struct   QOP_PC_ColorVector;
typedef struct QOP_IPC_DiracFermion_struct  QOP_PC_DiracFermion;
typedef struct QOP_IPC_GaugeField_struct    QOP_PC_GaugeField;
typedef struct QOP_IPC_Force_struct         QOP_PC_Force;

typedef struct QOP_IPC_FermionLinksAsqtad_struct  QOP_PC_FermionLinksAsqtad;
typedef struct QOP_IPC_FermionLinksHisq_struct    QOP_PC_FermionLinksHisq;
typedef struct QOP_IPC_FermionLinksWilson_struct  QOP_PC_FermionLinksWilson;
typedef struct QOP_IPC_FermionLinksDW_struct      QOP_PC_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

#define QOP_PC_qla_type_V QLA_PC_ColorVector
#define QOP_PC_qla_type_D QLA_PC_DiracFermion
#define QOP_PC_qla_type_M QLA_PC_ColorMatrix
#define QOP_PC_raw_size(T) (QDP_sites_on_node*sizeof(QOP_PC_qla_type_##T))
#define QOP_PC_raw_size_V(evenodd) QOP_PC_raw_size(V)
#define QOP_PC_raw_size_D(evenodd) QOP_PC_raw_size(D)
#define QOP_PC_raw_size_G(evenodd) QOP_PC_raw_size(M)
#define QOP_PC_raw_size_F(evenodd) QOP_PC_raw_size(M)

#define QOP_PC_elem(T, raw, i, ...) QLA_PC_elem_##T(((QOP_PC_qla_type_##T *)raw)[i], __ARGV__)
#define QOP_PC_set(T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_PC_elem(T, raw, i, __ARGV__), re, im)
#define QOP_PC_raw_set_V(raw, evenodd, i, ic, re, im) QOP_PC_set(V, raw, i, re, im, ic)
#define QOP_PC_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_PC_set(D, raw, i, re, im, ic, is)
#define QOP_PC_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_PC_set(M, raw, i, re, im, ic, jc)
#define QOP_PC_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_PC_set(M, raw, i, re, im, ic, jc)
#define QOP_PC_get(T, re, im, raw, i, ...) {			\
    QLA_P_Complex _c = QOP_PC_elem(T, raw, i, __ARGV__);	\
    re = QLA_real(_c); im = QLA_imag(_c);			\
  }
#define QOP_PC_raw_get_V(re, im, raw, evenodd, i, ic)     QOP_PC_set(V, re, im, raw, i, ic)
#define QOP_PC_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_PC_set(D, re, im, raw, i, ic, is)
#define QOP_PC_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_PC_set(M, re, im, raw, i, ic, jc)
#define QOP_PC_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_PC_set(M, re, im, raw, i, ic, jc)

/* create a QOP field with a copy of the raw source field */
QOP_PC_ColorVector  *QOP_PC_create_V_from_raw(NCPROT QOP_P_Real *src, QOP_evenodd_t evenodd);
QOP_PC_DiracFermion *QOP_PC_create_D_from_raw(NCPROT QOP_P_Real *src, QOP_evenodd_t evenodd);
QOP_PC_GaugeField   *QOP_PC_create_G_from_raw(NCPROT QOP_P_Real *links[], QOP_evenodd_t evenodd);
QOP_PC_Force        *QOP_PC_create_F_from_raw(NCPROT QOP_P_Real *force[], QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_PC_extract_V_to_raw(QOP_P_Real *dest, QOP_PC_ColorVector *src, QOP_evenodd_t evenodd);
void QOP_PC_extract_D_to_raw(QOP_P_Real *dest, QOP_PC_DiracFermion *src, QOP_evenodd_t evenodd);
void QOP_PC_extract_G_to_raw(QOP_P_Real *dest[], QOP_PC_GaugeField *src, QOP_evenodd_t evenodd);
void QOP_PC_extract_F_to_raw(QOP_P_Real *dest[], QOP_PC_Force *src, QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_PC_destroy_V(QOP_PC_ColorVector *field);
void QOP_PC_destroy_D(QOP_PC_DiracFermion *field);
void QOP_PC_destroy_G(QOP_PC_GaugeField *field);
void QOP_PC_destroy_F(QOP_PC_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_PC_ColorVector  *QOP_PC_convert_V_from_raw(NCPROT QOP_P_Real *src, QOP_evenodd_t evenodd);
QOP_PC_DiracFermion *QOP_PC_convert_D_from_raw(NCPROT QOP_P_Real *src, QOP_evenodd_t evenodd);
QOP_PC_GaugeField   *QOP_PC_convert_G_from_raw(NCPROT QOP_P_Real *links[], QOP_evenodd_t evenodd);
QOP_PC_Force        *QOP_PC_convert_F_from_raw(NCPROT QOP_P_Real *force[], QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
QOP_P_Real  *QOP_PC_convert_V_to_raw(QOP_PC_ColorVector *src, QOP_evenodd_t evenodd);
QOP_P_Real  *QOP_PC_convert_D_to_raw(QOP_PC_DiracFermion *src, QOP_evenodd_t evenodd);
QOP_P_Real **QOP_PC_convert_G_to_raw(QOP_PC_GaugeField *src, QOP_evenodd_t evenodd);
QOP_P_Real **QOP_PC_convert_F_to_raw(QOP_PC_Force *src, QOP_evenodd_t evenodd);

QOP_PC_ColorVector  *QOP_PC_create_V_from_qdp(QDP_PC_ColorVector *src);
QOP_PC_DiracFermion *QOP_PC_create_D_from_qdp(QDP_PC_DiracFermion *src);
QOP_PC_GaugeField   *QOP_PC_create_G_from_qdp(QDP_PC_ColorMatrix *src[]);
QOP_PC_Force        *QOP_PC_create_F_from_qdp(QDP_PC_ColorMatrix *src[]);

void QOP_PC_extract_V_to_qdp(QDP_PC_ColorVector *d, QOP_PC_ColorVector *src);
void QOP_PC_extract_D_to_qdp(QDP_PC_DiracFermion *d, QOP_PC_DiracFermion *src);
void QOP_PC_extract_G_to_qdp(QDP_PC_ColorMatrix *d[], QOP_PC_GaugeField *src);
void QOP_PC_extract_F_to_qdp(QDP_PC_ColorMatrix *d[], QOP_PC_Force *src);

QOP_PC_ColorVector  *QOP_PC_convert_V_from_qdp(QDP_PC_ColorVector *src);
QOP_PC_DiracFermion *QOP_PC_convert_D_from_qdp(QDP_PC_DiracFermion *src);
QOP_PC_GaugeField   *QOP_PC_convert_G_from_qdp(QDP_PC_ColorMatrix *src[]);
QOP_PC_Force        *QOP_PC_convert_F_from_qdp(QDP_PC_ColorMatrix *src[]);

QDP_PC_ColorVector   *QOP_PC_convert_V_to_qdp(QOP_PC_ColorVector *src);
QDP_PC_DiracFermion  *QOP_PC_convert_D_to_qdp(QOP_PC_DiracFermion *src);
QDP_PC_ColorMatrix  **QOP_PC_convert_G_to_qdp(QOP_PC_GaugeField *src);
QDP_PC_ColorMatrix  **QOP_PC_convert_F_to_qdp(QOP_PC_Force *src);


  /********************/
  /*  Gauge routines  */
  /********************/

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */
void QOP_PC_rephase_G(QOP_PC_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_PC_rephase_G_qdp(QDP_PC_ColorMatrix *links[],
			  int *r0,
			  QOP_bc_t *bc,
			  QOP_staggered_sign_t *sign);

void QOP_PC_smear_fat7l_qdp(QOP_info_t *info, QDP_PC_ColorMatrix *sg[],
			    QDP_PC_ColorMatrix *g[],
			    QOP_asqtad_coeffs_t *coeffs);

void QOP_PC_gauge_deriv_multi_qdp(QOP_info_t *info,
				  QDP_PC_ColorMatrix *deriv[],
				  QOP_PC_GaugeField *g[],
				  QDP_PC_ColorMatrix **chain[],
				  int n, int doLastScale);

void QOP_PC_gauge_force_multi_qdp(QOP_info_t *info, QDP_PC_ColorMatrix *f[],
				  QOP_PC_GaugeField *g[],
				  QDP_PC_ColorMatrix **chain[], int n);

void QOP_PC_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_PC_GaugeField *gauge,
					QOP_P_Real *acts, QOP_P_Real *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_PC_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info,
					   QOP_PC_GaugeField *gauge,
					   QDP_PC_ColorMatrix *deriv[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_P_Real eps, int doLastScale);

void QOP_PC_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, 
					   QOP_PC_GaugeField *gauge, 
					   QDP_PC_ColorMatrix *force[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_P_Real eps);

void QOP_PC_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_PC_GaugeField *gauge, 
				       QOP_PC_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       QOP_P_Real eps);

void QOP_PC_projectU_qdp(QOP_info_t *info,
			 QDP_PC_ColorMatrix *pU,
			 QDP_PC_ColorMatrix *U);

void QOP_PC_projectU_deriv_qdp(QOP_info_t *info,
			       QDP_PC_ColorMatrix *f,
			       QDP_PC_ColorMatrix *pU,
			       QDP_PC_ColorMatrix *U,
			       QDP_PC_ColorMatrix *chain);

void QOP_PC_u3reunit(QOP_info_t *info, QDP_PC_ColorMatrix *U,
		     QDP_PC_ColorMatrix *V);

void QOP_PC_su3reunit(QOP_info_t *info, QDP_PC_ColorMatrix *U,
		      QDP_PC_ColorMatrix *Ur);

void QOP_PC_hisq_force_multi_reunit(QOP_info_t *info,
				    QDP_PC_ColorMatrix *gf[4],
				    QDP_PC_ColorMatrix *force_accum[4],
				    QDP_PC_ColorMatrix *force_accum_old[4]);

void QOP_PC_staples(QOP_info_t *info, int nout, int nin,
		    QDP_PC_ColorMatrix *out[], QDP_PC_ColorMatrix *in[],
		    int nstaples[], int *topdir[], int *sidedir[],
		    int *toplinknum[], int *sidelinknum[], QOP_P_Real *coef[]);

void QOP_PC_staples_deriv(QOP_info_t *info, int nout, int nin,
			  QDP_PC_ColorMatrix *deriv[],
			  QDP_PC_ColorMatrix *chain[],
			  QDP_PC_ColorMatrix *in[],
			  int nstaples[], int *topdir[], int *sidedir[],
			  int *toplinknum[], int *sidelinknum[],
			  QOP_P_Real *coef[]);

  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_PC_FermionLinksAsqtad *
  QOP_PC_asqtad_create_L_from_raw(NCPROT QOP_P_Real *fatlinks[],
				  QOP_P_Real *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_PC_FermionLinksAsqtad *
  QOP_PC_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_PC_GaugeField *gauge);

void QOP_PC_asqtad_extract_L_to_raw(QOP_P_Real *fatlinks[],
				    QOP_P_Real *longlinks[],
				    QOP_PC_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_PC_asqtad_destroy_L(QOP_PC_FermionLinksAsqtad *field);

QOP_PC_FermionLinksAsqtad *
  QOP_PC_asqtad_convert_L_from_raw(NCPROT QOP_P_Real *fatlinks[],
				   QOP_P_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_PC_asqtad_convert_L_to_raw(QOP_P_Real ***fatlinks,
				    QOP_P_Real ***longlinks,
				    QOP_PC_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_PC_asqtad_load_L_from_raw(QOP_PC_FermionLinksAsqtad *asqtad,
				   QOP_P_Real *fatlinks[],
				   QOP_P_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_PC_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_PC_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_PC_GaugeField *gauge);

void QOP_PC_asqtad_rephase_L(QOP_PC_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

QOP_PC_FermionLinksAsqtad *
  QOP_PC_asqtad_create_L_from_qdp(QDP_PC_ColorMatrix *fatlinks[],
				  QDP_PC_ColorMatrix *longlinks[]);

void QOP_PC_asqtad_extract_L_to_qdp(QDP_PC_ColorMatrix *fatlinks[],
				    QDP_PC_ColorMatrix *longlinks[],
				    QOP_PC_FermionLinksAsqtad *src);

QOP_PC_FermionLinksAsqtad *
  QOP_PC_asqtad_convert_L_from_qdp(QDP_PC_ColorMatrix *fatlinks[],
				   QDP_PC_ColorMatrix *longlinks[]);

void QOP_PC_asqtad_convert_L_to_qdp(QDP_PC_ColorMatrix ***fatlinks,
				    QDP_PC_ColorMatrix ***longlinks,
				    QOP_PC_FermionLinksAsqtad *src);

void QOP_PC_asqtad_load_L_from_qdp(QOP_PC_FermionLinksAsqtad *asqtad,
				   QDP_PC_ColorMatrix *fatlinks[],
				   QDP_PC_ColorMatrix *longlinks[]);

void QOP_PC_asqtad_rephase_field_L_qdp(QOP_PC_FermionLinksAsqtad *fla,
				       QDP_P_Complex *fatphase[],
				       QDP_P_Complex *longphase[]);


  /* inverter routines */

void QOP_PC_asqtad_dslash(QOP_info_t *info,
			  QOP_PC_FermionLinksAsqtad *asqtad,
			  QOP_P_Real mass,
			  QOP_PC_ColorVector *out,
			  QOP_PC_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_PC_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_PC_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_PC_ColorVector *out,
			      QOP_PC_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_PC_asqtad_diaginv(QOP_info_t *info,
			   QOP_PC_FermionLinksAsqtad *asqtad,
			   QOP_P_Real mass,
			   QOP_PC_ColorVector *out,
			   QOP_PC_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_PC_asqtad_invert(QOP_info_t *info,
			  QOP_PC_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_P_Real mass,
			  QOP_PC_ColorVector *out_pt,
			  QOP_PC_ColorVector *in_pt);

void QOP_PC_asqtad_invert_threaded(QOP_info_t *info,
				   QOP_PC_FermionLinksAsqtad *asqtad,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg,
				   QOP_P_Real mass,
				   QOP_PC_ColorVector *out_pt,
				   QOP_PC_ColorVector *in_pt,
				   int nthreads);

void QOP_PC_asqtad_invert_multi(QOP_info_t *info,
				QOP_PC_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_P_Real *masses[],
				int nmass[],
				QOP_PC_ColorVector **out_pt[],
				QOP_PC_ColorVector *in_pt[],
				int nsrc);

void QOP_PC_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_PC_FermionLinksAsqtad *asqtad,
			      QOP_P_Real mass,
			      QDP_PC_ColorVector *out,
			      QDP_PC_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_PC_asqtad_dslash_dir_qdp(QOP_info_t *info,
				  QOP_PC_FermionLinksAsqtad *asqtad,
				  int dir, int fb,
				  double wtfat, double wtlong,
				  QDP_PC_ColorVector *out,
				  QDP_PC_ColorVector *in,
				  QOP_evenodd_t eo_out);

void QOP_PC_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_PC_FermionLinksAsqtad *asqtad,
			       QOP_P_Real mass,
			       QDP_PC_ColorVector *out,
			       QDP_PC_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_PC_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_PC_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_P_Real mass,
			      QDP_PC_ColorVector *out,
			      QDP_PC_ColorVector *in);

void QOP_PC_asqtad_invert_threaded_qdp(QOP_info_t *info,
				       QOP_PC_FermionLinksAsqtad *asqtad,
				       QOP_invert_arg_t *inv_arg,
				       QOP_resid_arg_t *res_arg,
				       QOP_P_Real mass,
				       QDP_PC_ColorVector *out,
				       QDP_PC_ColorVector *in,
				       int nthreads);

void QOP_PC_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_PC_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_P_Real *masses[],
				    int nmass[],
				    QDP_PC_ColorVector **out[],
				    QDP_PC_ColorVector *in[],
				    int nsrc);

void QOP_PC_asqtad_get_eigcg(QOP_PC_FermionLinksAsqtad *asqtad,
			     QDP_PC_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

  /* fermion force routines */

void QOP_PC_asqtad_deriv(QOP_info_t *info, QDP_PC_ColorMatrix *gauge[],
			 QDP_PC_ColorMatrix *force[],
			 QOP_asqtad_coeffs_t *coef,
			 QDP_PC_ColorMatrix *mid_fat[],
			 QDP_PC_ColorMatrix *mid_naik[]);

void QOP_PC_asqtad_force(QOP_info_t *info,
			 QOP_PC_GaugeField *gauge,
			 QOP_PC_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 QOP_P_Real eps,
			 QOP_PC_ColorVector *in_pt);

void QOP_PC_asqtad_force_multi(QOP_info_t *info,
			       QOP_PC_GaugeField *gauge,
			       QOP_PC_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       QOP_P_Real eps[],
			       QOP_PC_ColorVector *in_pt[],
			       int nsrc);

void QOP_PC_asqtad_force_multi_qdp(QOP_info_t *info,
				   QDP_PC_ColorMatrix *links[],
				   QDP_PC_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_P_Real eps[],
				   QDP_PC_ColorVector *in_pt[],
				   int nsrc);

void QOP_PC_asqtad_deriv_multi_qdp(QOP_info_t *info,
				   QDP_PC_ColorMatrix *links[],
				   QDP_PC_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_P_Real eps[],
				   QDP_PC_ColorVector *in_pt[],
				   int nsrc);

  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* single precision */

QOP_PC_FermionLinksHisq *
  QOP_PC_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_PC_GaugeField *gauge);

void QOP_PC_hisq_destroy_L(QOP_PC_FermionLinksHisq *field);

QOP_PC_FermionLinksAsqtad **
  QOP_PC_get_asqtad_links_from_hisq(QOP_PC_FermionLinksHisq *hl);
  
QOP_PC_FermionLinksAsqtad *
  QOP_PC_get_asqtad_deps_links_from_hisq(QOP_PC_FermionLinksHisq *hl);

  /* fermion force routines */

void QOP_PC_hisq_force_multi(QOP_info_t *info,
			     QOP_PC_FermionLinksHisq *flh,
			     QOP_PC_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     QOP_P_Real eps[],
			     QOP_PC_ColorVector *in_pt[],
			     int *n_orders_naik);

void QOP_PC_hisq_deriv_multi_qdp(QOP_info_t *info,
				 QOP_PC_FermionLinksHisq *flh,
				 QDP_PC_ColorMatrix *deriv[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_P_Real eps[],
				 QDP_PC_ColorVector *in_pt[],
				 int *n_orders_naik,
				 int doLastScale);

void QOP_PC_hisq_force_multi_qdp(QOP_info_t *info,
				 QOP_PC_FermionLinksHisq *flh,
				 QDP_PC_ColorMatrix *force[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_P_Real eps[],
				 QDP_PC_ColorVector *in_pt[],
				 int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_create_L_from_raw(NCPROT QOP_P_Real *links[], QOP_P_Real *clov,
				  QOP_evenodd_t evenodd);

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_PC_GaugeField *gauge);

void QOP_PC_wilson_extract_L_to_raw(QOP_P_Real *links[], QOP_P_Real *clov,
				    QOP_PC_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_PC_wilson_destroy_L(QOP_PC_FermionLinksWilson *field);

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_convert_L_from_raw(NCPROT QOP_P_Real *links[],
				   QOP_P_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_PC_wilson_convert_L_to_raw(QOP_P_Real ***links, QOP_P_Real **clov,
				    QOP_PC_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_PC_GaugeField *gauge);

QOP_PC_GaugeField *
  QOP_PC_wilson_convert_L_to_G(QOP_PC_FermionLinksWilson *links);

void QOP_PC_wilson_load_L_from_raw(QOP_PC_FermionLinksWilson *wilson,
				   QOP_P_Real *links[], QOP_P_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_PC_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_PC_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_PC_GaugeField *gauge);

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_create_L_from_qdp(QDP_PC_ColorMatrix *links[],
				  QDP_PC_DiracPropagator *clov);

void QOP_PC_wilson_extract_L_to_qdp(QDP_PC_ColorMatrix *links[],
				    QDP_PC_DiracPropagator *clov,
				    QOP_PC_FermionLinksWilson *src);

QOP_PC_FermionLinksWilson *
  QOP_PC_wilson_convert_L_from_qdp(QDP_PC_ColorMatrix *links[],
				   QDP_PC_DiracPropagator *clov);

void QOP_PC_wilson_convert_L_to_qdp(QDP_PC_ColorMatrix ***links,
				    QDP_PC_DiracPropagator **clov,
				    QOP_PC_FermionLinksWilson *src);

void QOP_PC_wilson_load_L_from_qdp(QOP_PC_FermionLinksWilson *wilson,
				   QDP_PC_ColorMatrix *links[],
				   QDP_PC_DiracPropagator *clov);


  /* inverter routines */

void QOP_PC_wilson_dslash(QOP_info_t *info,
			  QOP_PC_FermionLinksWilson *flw,
			  QOP_P_Real kappa,
			  int sign,
			  QOP_PC_DiracFermion *out,
			  QOP_PC_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_PC_wilson_diaginv(QOP_info_t *info,
			   QOP_PC_FermionLinksWilson *flw,
			   QOP_P_Real kappa,
			   QOP_PC_DiracFermion *out,
			   QOP_PC_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_PC_wilson_invert(QOP_info_t *info,
			  QOP_PC_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_P_Real kappa,
			  QOP_PC_DiracFermion *out_pt,
			  QOP_PC_DiracFermion *in_pt);

void QOP_PC_wilson_invert_multi(QOP_info_t *info,
				QOP_PC_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_P_Real *kappas[],
				int nkappa[],
				QOP_PC_DiracFermion **out_pt[],
				QOP_PC_DiracFermion *in_pt[],
				int nsrc);

void QOP_PC_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_PC_FermionLinksWilson *flw,
			      QOP_P_Real kappa,
			      int sign,
			      QDP_PC_DiracFermion *out,
			      QDP_PC_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_PC_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_PC_FermionLinksWilson *flw,
			       QOP_P_Real kappa,
			       QDP_PC_DiracFermion *out,
			       QDP_PC_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_PC_wilson_invert_qdp(QOP_info_t *info,
			      QOP_PC_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_P_Real kappa,
			      QDP_PC_DiracFermion *out_pt,
			      QDP_PC_DiracFermion *in_pt);

void QOP_PC_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_PC_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_P_Real *kappas[],
				    int nkappa[],
				    QDP_PC_DiracFermion **out_pt[],
				    QDP_PC_DiracFermion *in_pt[],
				    int nsrc);

void QOP_PC_wilson_invert_ne_qdp(QOP_info_t *info,
				 QOP_PC_FermionLinksWilson *flw,
				 QOP_invert_arg_t *inv_arg,
				 QOP_resid_arg_t *res_arg,
				 QOP_P_Real kappa,
				 QDP_PC_DiracFermion *out,
				 QDP_PC_DiracFermion *in);

  /* fermion force routines */

void QOP_PC_wilson_deriv_multi_qdp(QOP_info_t *info,
				   QOP_PC_FermionLinksWilson *flw,
				   QDP_PC_ColorMatrix *deriv[],
				   QOP_P_Real eps[],
				   QDP_PC_DiracFermion *x[],
				   QDP_PC_DiracFermion *y[],
				   int n);

void QOP_PC_wilson_force_multi_qdp(QOP_info_t *info,
				   QOP_PC_FermionLinksWilson *flw,
				   QDP_PC_ColorMatrix *force[],
				   QOP_P_Real eps[],
				   QDP_PC_DiracFermion *x[],
				   QDP_PC_DiracFermion *y[],
				   int n);

void QOP_PC_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
					QOP_PC_FermionLinksWilson *flw,
					QDP_PC_ColorMatrix *deriv[],
					QOP_P_Real kappa[],
					QOP_P_Real eps[],
					QDP_PC_DiracFermion *x[],
					QDP_PC_DiracFermion *y[],
					int n);

void QOP_PC_wilson_force_prec_multi_qdp(QOP_info_t *info,
					QOP_PC_FermionLinksWilson *flw,
					QDP_PC_ColorMatrix *force[],
					QOP_P_Real kappa[],
					QOP_P_Real eps[],
					QDP_PC_DiracFermion *x[],
					QDP_PC_DiracFermion *y[],
					int n);

  // new fermilab action IFLA -- added by bugra --------------- :

void QOP_PC_wilson_ifla_dslash(QOP_info_t *info,
			       QOP_PC_FermionLinksWilson *flw,
			       QOP_P_Real kappa,
			       int sign,
			       QOP_wilson_ifla_coeffs_t *coeffs,
			       QOP_PC_DiracFermion *out,
			       QOP_PC_DiracFermion *in,
			       QOP_evenodd_t eo_out,
			       QOP_evenodd_t eo_in);

void QOP_PC_wilson_ifla_dslash_qdp(QOP_info_t *info,
				   QOP_PC_FermionLinksWilson *flw,
				   QOP_P_Real kappa,
				   int sign,
				   QOP_wilson_ifla_coeffs_t *coeffs,
				   QDP_PC_DiracFermion *out,
				   QDP_PC_DiracFermion *in,
				   QOP_evenodd_t eo_out,
				   QOP_evenodd_t eo_in);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_PC_FermionLinksDW *
  QOP_PC_dw_create_L_from_raw(NCPROT QOP_P_Real *links[],
			      QOP_evenodd_t evenodd);

QOP_PC_FermionLinksDW *
  QOP_PC_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_PC_GaugeField *gauge);

void QOP_PC_dw_extract_L_to_raw(QOP_P_Real *links[],
				QOP_PC_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_PC_dw_destroy_L(QOP_PC_FermionLinksDW *field);

QOP_PC_FermionLinksDW *
  QOP_PC_dw_convert_L_from_raw(NCPROT QOP_P_Real *links[],
			       QOP_evenodd_t evenodd);

void QOP_PC_dw_convert_L_to_raw(QOP_P_Real ***links,
				QOP_PC_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_PC_FermionLinksDW *
  QOP_PC_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_PC_GaugeField *gauge);

QOP_PC_GaugeField *
  QOP_PC_dw_convert_L_to_G(QOP_PC_FermionLinksDW *links);

void QOP_PC_dw_load_L_from_raw(QOP_PC_FermionLinksDW *dw,
			       QOP_P_Real *links[], QOP_evenodd_t evenodd);

void QOP_PC_dw_load_L_from_G(QOP_info_t *info,
			     QOP_PC_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_PC_GaugeField *gauge);

QOP_PC_FermionLinksDW *
  QOP_PC_dw_create_L_from_qdp(QDP_PC_ColorMatrix *links[]);

void QOP_PC_dw_extract_L_to_qdp(QDP_PC_ColorMatrix *links[],
				QOP_PC_FermionLinksDW *src);

QOP_PC_FermionLinksDW *
  QOP_PC_dw_convert_L_from_qdp(QDP_PC_ColorMatrix *links[]);

void QOP_PC_dw_convert_L_to_qdp(QDP_PC_ColorMatrix ***links,
				QOP_PC_FermionLinksDW *src);

void QOP_PC_dw_load_L_from_qdp(QOP_PC_FermionLinksDW *dw,
			       QDP_PC_ColorMatrix *links[]);

  /* inverter routines */

void QOP_PC_dw_dslash(QOP_info_t *info,
		      QOP_PC_FermionLinksDW *links,
		      QOP_P_Real M5,
		      QOP_P_Real m,
		      int sign,
		      QOP_PC_DiracFermion *out_pt[],
		      QOP_PC_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_PC_dw_dslash2(QOP_info_t *info,
		       QOP_PC_FermionLinksDW *links,
		       QOP_P_Real M5,
		       QOP_P_Real m,
		       QOP_PC_DiracFermion *out_pt[],
		       QOP_PC_DiracFermion *in_pt[],
		       int Ls,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in);

void QOP_PC_dw_invert(QOP_info_t *info,
		      QOP_PC_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QOP_P_Real M5,
		      QOP_P_Real m,
		      QOP_PC_DiracFermion *out_pt[],
		      QOP_PC_DiracFermion *in_pt[],
		      int Ls);

void QOP_PC_dw_invert_multi(QOP_info_t *info,
			    QOP_PC_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_P_Real *M5[],
			    QOP_P_Real *m[],
			    int nmass[],
			    QOP_PC_DiracFermion ***out_pt[],
			    QOP_PC_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_PC_dw_dslash_qdp(QOP_info_t *info,
			  QOP_PC_FermionLinksDW *links,
			  QOP_P_Real M5, 
			  QOP_P_Real m,
			  int sign,
			  QDP_PC_DiracFermion *out_pt[],
			  QDP_PC_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_PC_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_PC_FermionLinksDW *links,
			   QOP_P_Real M5,
			   QOP_P_Real m,
			   QDP_PC_DiracFermion *out_pt[],
			   QDP_PC_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);

void QOP_PC_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_PC_FermionLinksDW *fldw,
			   QOP_P_Real M5,
			   QOP_P_Real m,
			   QDP_PC_DiracFermion **out,
			   QDP_PC_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_PC_dw_invert_qdp(QOP_info_t *info,
			  QOP_PC_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_P_Real M5,
			  QOP_P_Real m,
			  QDP_PC_DiracFermion *out[],
			  QDP_PC_DiracFermion *in[],
			  int Ls);

void QOP_PC_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_PC_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_P_Real *M5[],
				QOP_P_Real *m[],
				int nmass[],
				QDP_PC_DiracFermion ***out[],
				QDP_PC_DiracFermion **in[],
				int nsrc,
				int Ls);

 /* fermion force routines */

void QOP_PC_dw_force(QOP_info_t *info,
		     QOP_PC_GaugeField *gauge,
		     QOP_PC_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     QOP_P_Real eps,
		     QOP_PC_DiracFermion *in_pt);

void QOP_PC_dw_force_multi(QOP_info_t *info,
			   QOP_PC_GaugeField *gauge,
			   QOP_PC_Force *force,
			   QOP_dw_coeffs_t *coef,
			   QOP_P_Real eps[],
			   QOP_PC_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == _QOP_Precision
#  if QOP_Colors == _QOP_Colors
#    include <qop_pc_generic.h>
#  endif
#  include <qop_pc_precision_generic.h>
#endif
#if QOP_Colors == _QOP_Colors
#  include <qop_pc_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_PC_H */
