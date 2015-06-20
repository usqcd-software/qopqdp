// DO NOT EDIT
// generated from qop_pc.h
#ifndef _QOP_FN_H
#define _QOP_FN_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QOP_FN_ColorVector_struct   QOP_FN_ColorVector;
typedef struct QOP_FN_DiracFermion_struct  QOP_FN_DiracFermion;
typedef struct QOP_FN_GaugeField_struct    QOP_FN_GaugeField;
typedef struct QOP_FN_Force_struct         QOP_FN_Force;

typedef struct QOP_FN_FermionLinksAsqtad_struct  QOP_FN_FermionLinksAsqtad;
typedef struct QOP_FN_FermionLinksHisq_struct    QOP_FN_FermionLinksHisq;
typedef struct QOP_FN_FermionLinksWilson_struct  QOP_FN_FermionLinksWilson;
typedef struct QOP_FN_FermionLinksDW_struct      QOP_FN_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

#define QOP_FN_qla_type_V QLA_FN_ColorVector
#define QOP_FN_qla_type_D QLA_FN_DiracFermion
#define QOP_FN_qla_type_M QLA_FN_ColorMatrix
#define QOP_FN_raw_size(T) (QDP_sites_on_node*sizeof(QOP_FN_qla_type_##T))
#define QOP_FN_raw_size_V(evenodd) QOP_FN_raw_size(V)
#define QOP_FN_raw_size_D(evenodd) QOP_FN_raw_size(D)
#define QOP_FN_raw_size_G(evenodd) QOP_FN_raw_size(M)
#define QOP_FN_raw_size_F(evenodd) QOP_FN_raw_size(M)

#define QOP_FN_elem(T, raw, i, ...) QLA_FN_elem_##T(((QOP_FN_qla_type_##T *)raw)[i], __ARGV__)
#define QOP_FN_set(T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_FN_elem(T, raw, i, __ARGV__), re, im)
#define QOP_FN_raw_set_V(raw, evenodd, i, ic, re, im) QOP_FN_set(V, raw, i, re, im, ic)
#define QOP_FN_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_FN_set(D, raw, i, re, im, ic, is)
#define QOP_FN_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_FN_set(M, raw, i, re, im, ic, jc)
#define QOP_FN_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_FN_set(M, raw, i, re, im, ic, jc)
#define QOP_FN_get(T, re, im, raw, i, ...) {			\
    QLA_F_Complex _c = QOP_FN_elem(T, raw, i, __ARGV__);	\
    re = QLA_real(_c); im = QLA_imag(_c);			\
  }
#define QOP_FN_raw_get_V(re, im, raw, evenodd, i, ic)     QOP_FN_set(V, re, im, raw, i, ic)
#define QOP_FN_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_FN_set(D, re, im, raw, i, ic, is)
#define QOP_FN_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_FN_set(M, re, im, raw, i, ic, jc)
#define QOP_FN_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_FN_set(M, re, im, raw, i, ic, jc)

/* create a QOP field with a copy of the raw source field */
QOP_FN_ColorVector  *QOP_FN_create_V_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_FN_DiracFermion *QOP_FN_create_D_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_FN_GaugeField   *QOP_FN_create_G_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *links[], QOP_evenodd_t evenodd);
QOP_FN_Force        *QOP_FN_create_F_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *force[], QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_FN_extract_V_to_raw(QOP_F_Real *dest, QOP_FN_ColorVector *src, QOP_evenodd_t evenodd);
void QOP_FN_extract_D_to_raw(QOP_F_Real *dest, QOP_FN_DiracFermion *src, QOP_evenodd_t evenodd);
void QOP_FN_extract_G_to_raw(QOP_F_Real *dest[], QOP_FN_GaugeField *src, QOP_evenodd_t evenodd);
void QOP_FN_extract_F_to_raw(QOP_F_Real *dest[], QOP_FN_Force *src, QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_FN_destroy_V(QOP_FN_ColorVector *field);
void QOP_FN_destroy_D(QOP_FN_DiracFermion *field);
void QOP_FN_destroy_G(QOP_FN_GaugeField *field);
void QOP_FN_destroy_F(QOP_FN_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_FN_ColorVector  *QOP_FN_convert_V_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_FN_DiracFermion *QOP_FN_convert_D_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *src, QOP_evenodd_t evenodd);
QOP_FN_GaugeField   *QOP_FN_convert_G_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *links[], QOP_evenodd_t evenodd);
QOP_FN_Force        *QOP_FN_convert_F_from_raw(int nc, QDP_Lattice *lat, QOP_F_Real *force[], QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
QOP_F_Real  *QOP_FN_convert_V_to_raw(QDP_Lattice *lat, QOP_FN_ColorVector *src, QOP_evenodd_t evenodd);
QOP_F_Real  *QOP_FN_convert_D_to_raw(QDP_Lattice *lat, QOP_FN_DiracFermion *src, QOP_evenodd_t evenodd);
QOP_F_Real **QOP_FN_convert_G_to_raw(QDP_Lattice *lat, QOP_FN_GaugeField *src, QOP_evenodd_t evenodd);
QOP_F_Real **QOP_FN_convert_F_to_raw(QDP_Lattice *lat, QOP_FN_Force *src, QOP_evenodd_t evenodd);

QOP_FN_ColorVector  *QOP_FN_create_V_from_qdp(QDP_FN_ColorVector *src);
QOP_FN_DiracFermion *QOP_FN_create_D_from_qdp(QDP_FN_DiracFermion *src);
QOP_FN_GaugeField   *QOP_FN_create_G_from_qdp(QDP_FN_ColorMatrix *src[]);
QOP_FN_Force        *QOP_FN_create_F_from_qdp(QDP_FN_ColorMatrix *src[]);

void QOP_FN_extract_V_to_qdp(QDP_FN_ColorVector *d, QOP_FN_ColorVector *src);
void QOP_FN_extract_D_to_qdp(QDP_FN_DiracFermion *d, QOP_FN_DiracFermion *src);
void QOP_FN_extract_G_to_qdp(QDP_FN_ColorMatrix *d[], QOP_FN_GaugeField *src);
void QOP_FN_extract_F_to_qdp(QDP_FN_ColorMatrix *d[], QOP_FN_Force *src);

QOP_FN_ColorVector  *QOP_FN_convert_V_from_qdp(QDP_FN_ColorVector *src);
QOP_FN_DiracFermion *QOP_FN_convert_D_from_qdp(QDP_FN_DiracFermion *src);
QOP_FN_GaugeField   *QOP_FN_convert_G_from_qdp(QDP_FN_ColorMatrix *src[]);
QOP_FN_Force        *QOP_FN_convert_F_from_qdp(QDP_FN_ColorMatrix *src[]);

QDP_FN_ColorVector   *QOP_FN_convert_V_to_qdp(QOP_FN_ColorVector *src);
QDP_FN_DiracFermion  *QOP_FN_convert_D_to_qdp(QOP_FN_DiracFermion *src);
QDP_FN_ColorMatrix  **QOP_FN_convert_G_to_qdp(QOP_FN_GaugeField *src);
QDP_FN_ColorMatrix  **QOP_FN_convert_F_to_qdp(QOP_FN_Force *src);


  /********************/
  /*  Gauge routines  */
  /********************/

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */
void QOP_FN_rephase_G(QOP_FN_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_FN_rephase_G_qdp(QDP_FN_ColorMatrix *links[],
			  int *r0,
			  QOP_bc_t *bc,
			  QOP_staggered_sign_t *sign);

void QOP_FN_smear_fat7l_qdp(QOP_info_t *info, QDP_FN_ColorMatrix *sg[],
			    QDP_FN_ColorMatrix *g[],
			    QOP_asqtad_coeffs_t *coeffs);

void QOP_FN_gauge_deriv_multi_qdp(QOP_info_t *info,
				  QDP_FN_ColorMatrix *deriv[],
				  QOP_FN_GaugeField *g[],
				  QDP_FN_ColorMatrix **chain[],
				  int n, int doLastScale);

void QOP_FN_gauge_force_multi_qdp(QOP_info_t *info, QDP_FN_ColorMatrix *f[],
				  QOP_FN_GaugeField *g[],
				  QDP_FN_ColorMatrix **chain[], int n);

void QOP_FN_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_FN_GaugeField *gauge,
					QOP_F_Real *acts, QOP_F_Real *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_FN_symanzik_1loop_gauge_action_qdp(QOP_info_t *info,
					    QDP_FN_ColorMatrix *links[],
					    QOP_F_Real *acts, QOP_F_Real *actt,
					    QOP_gauge_coeffs_t *coeffs);

void QOP_FN_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_FN_GaugeField *gauge, 
				       QOP_FN_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       QOP_F_Real eps);

void QOP_FN_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, 
					   QDP_FN_ColorMatrix *links[],
					   QDP_FN_ColorMatrix *force[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_F_Real eps);

void QOP_FN_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info,
					   QDP_FN_ColorMatrix *links[],
					   QDP_FN_ColorMatrix *deriv[],
					   QOP_gauge_coeffs_t *coeffs,
					   QOP_F_Real eps, int doLastScale);

void QOP_FN_symanzik_1loop_gauge_heatbath_qdp(QOP_info_t *info,
					      QDP_FN_ColorMatrix *links[],
					      QLA_F_Real beta,
					      QOP_gauge_coeffs_t *coeffs,
					      QDP_RandomState *rs0,
					      int nup, int nhb, int nover);

void QOP_FN_symanzik_1loop_gauge_staple_qdp(QOP_info_t *info,
					    QDP_FN_ColorMatrix *links[],
					    QDP_FN_ColorMatrix *staple,
					    int mu,
					    QOP_gauge_coeffs_t *coeffs,
					    QDP_Subset subs[], int subi);

void QOP_FN_projectU_qdp(QOP_info_t *info,
			 QDP_FN_ColorMatrix *pU,
			 QDP_FN_ColorMatrix *U);

void QOP_FN_projectU_deriv_qdp(QOP_info_t *info,
			       QDP_FN_ColorMatrix *f,
			       QDP_FN_ColorMatrix *pU,
			       QDP_FN_ColorMatrix *U,
			       QDP_FN_ColorMatrix *chain);

void QOP_FN_u3reunit(QOP_info_t *info, QDP_FN_ColorMatrix *U,
		     QDP_FN_ColorMatrix *V);

void QOP_FN_su3reunit(QOP_info_t *info, QDP_FN_ColorMatrix *U,
		      QDP_FN_ColorMatrix *Ur);

void QOP_FN_hisq_force_multi_reunit(QOP_info_t *info,
				    QDP_FN_ColorMatrix *gf[4],
				    QDP_FN_ColorMatrix *force_accum[4],
				    QDP_FN_ColorMatrix *force_accum_old[4]);

void QOP_FN_staples(QOP_info_t *info, int nout, int nin,
		    QDP_FN_ColorMatrix *out[], QDP_FN_ColorMatrix *in[],
		    int nstaples[], int *topdir[], int *sidedir[],
		    int *toplinknum[], int *sidelinknum[], QOP_F_Real *coef[]);

void QOP_FN_staples_deriv(QOP_info_t *info, int nout, int nin,
			  QDP_FN_ColorMatrix *deriv[],
			  QDP_FN_ColorMatrix *chain[],
			  QDP_FN_ColorMatrix *in[],
			  int nstaples[], int *topdir[], int *sidedir[],
			  int *toplinknum[], int *sidelinknum[],
			  QOP_F_Real *coef[]);

  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_create_L_from_raw(int nc, QDP_Lattice *lat, 
                                  QOP_F_Real *fatlinks[],
				  QOP_F_Real *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_FN_GaugeField *gauge);

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_create_L_from_G2(QOP_info_t *info,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_FN_GaugeField *gFat,
				 QOP_FN_GaugeField *gLong);

void QOP_FN_asqtad_extract_L_to_raw(QOP_F_Real *fatlinks[],
				    QOP_F_Real *longlinks[],
				    QOP_FN_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_FN_asqtad_destroy_L(QOP_FN_FermionLinksAsqtad *field);

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_convert_L_from_raw(int nc, QOP_F_Real *fatlinks[],
				   QOP_F_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_FN_asqtad_convert_L_to_raw(QOP_F_Real ***fatlinks,
				    QOP_F_Real ***longlinks,
				    QOP_FN_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_FN_asqtad_load_L_from_raw(QOP_FN_FermionLinksAsqtad *asqtad,
				   QOP_F_Real *fatlinks[],
				   QOP_F_Real *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_FN_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_FN_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_FN_GaugeField *gauge);

void QOP_FN_asqtad_load_L_from_G2(QOP_info_t *info,
				  QOP_FN_FermionLinksAsqtad *asqtad,
				  QOP_asqtad_coeffs_t *coeffs,
				  QOP_FN_GaugeField *gFat,
				  QOP_FN_GaugeField *gLong);

void QOP_FN_asqtad_rephase_L(QOP_FN_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_create_L_from_qdp(QDP_FN_ColorMatrix *fatlinks[],
				  QDP_FN_ColorMatrix *longlinks[]);

void QOP_FN_asqtad_extract_L_to_qdp(QDP_FN_ColorMatrix *fatlinks[],
				    QDP_FN_ColorMatrix *longlinks[],
				    QOP_FN_FermionLinksAsqtad *src);

QOP_FN_FermionLinksAsqtad *
  QOP_FN_asqtad_convert_L_from_qdp(QDP_FN_ColorMatrix *fatlinks[],
				   QDP_FN_ColorMatrix *longlinks[]);

void QOP_FN_asqtad_convert_L_to_qdp(QDP_FN_ColorMatrix ***fatlinks,
				    QDP_FN_ColorMatrix ***longlinks,
				    QOP_FN_FermionLinksAsqtad *src);

void QOP_FN_asqtad_load_L_from_qdp(QOP_FN_FermionLinksAsqtad *asqtad,
				   QDP_FN_ColorMatrix *fatlinks[],
				   QDP_FN_ColorMatrix *longlinks[]);

void QOP_FN_asqtad_rephase_field_L_qdp(QOP_FN_FermionLinksAsqtad *fla,
				       QDP_F_Complex *fatphase[],
				       QDP_F_Complex *longphase[]);


  /* inverter routines */

void QOP_FN_asqtad_dslash(QOP_info_t *info,
			  QOP_FN_FermionLinksAsqtad *asqtad,
			  QOP_F_Real mass,
			  QOP_FN_ColorVector *out,
			  QOP_FN_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_FN_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_FN_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_FN_ColorVector *out,
			      QOP_FN_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_FN_asqtad_diaginv(QOP_info_t *info,
			   QOP_FN_FermionLinksAsqtad *asqtad,
			   QOP_F_Real mass,
			   QOP_FN_ColorVector *out,
			   QOP_FN_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_FN_asqtad_ddagd(QOP_info_t *info,
			 QOP_FN_FermionLinksAsqtad *asqtad,
			 QOP_F_Real mass,
			 QDP_FN_ColorVector *out,
			 QDP_FN_ColorVector *in,
			 QOP_evenodd_t eo);

QOP_F_Real QOP_FN_asqtad_ddagd_norm2(QOP_info_t *info,
				     QOP_FN_FermionLinksAsqtad *asqtad,
				     QOP_F_Real mass,
				     QDP_FN_ColorVector *out,
				     QDP_FN_ColorVector *in,
				     QOP_evenodd_t eo);

void QOP_FN_asqtad_solve_multi_qdp(QOP_info_t *info,
				   QOP_FN_FermionLinksAsqtad *fla,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg[],
				   QOP_F_Real masses[],
				   QDP_FN_ColorVector *out[],
				   QDP_FN_ColorVector *in[],
				   int nsolve);

void QOP_FN_asqtad_invert(QOP_info_t *info,
			  QOP_FN_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real mass,
			  QOP_FN_ColorVector *out_pt,
			  QOP_FN_ColorVector *in_pt);

void QOP_FN_asqtad_invert_threaded(QOP_info_t *info,
				   QOP_FN_FermionLinksAsqtad *asqtad,
				   QOP_invert_arg_t *inv_arg,
				   QOP_resid_arg_t *res_arg,
				   QOP_F_Real mass,
				   QOP_FN_ColorVector *out_pt,
				   QOP_FN_ColorVector *in_pt,
				   int nthreads);

void QOP_FN_asqtad_invert_multi(QOP_info_t *info,
				QOP_FN_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *masses[],
				int nmass[],
				QOP_FN_ColorVector **out_pt[],
				QOP_FN_ColorVector *in_pt[],
				int nsrc);

void QOP_FN_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_FN_FermionLinksAsqtad *asqtad,
			      QOP_F_Real mass,
			      QDP_FN_ColorVector *out,
			      QDP_FN_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_FN_asqtad_dslash_dir_qdp(QOP_info_t *info,
				  QOP_FN_FermionLinksAsqtad *asqtad,
				  int dir, int fb,
				  double wtfat, double wtlong,
				  QDP_FN_ColorVector *out,
				  QDP_FN_ColorVector *in,
				  QOP_evenodd_t eo_out);

void QOP_FN_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_FN_FermionLinksAsqtad *asqtad,
			       QOP_F_Real mass,
			       QDP_FN_ColorVector *out,
			       QDP_FN_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_FN_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_FN_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_F_Real mass,
			      QDP_FN_ColorVector *out,
			      QDP_FN_ColorVector *in);

void QOP_FN_asqtad_invert_threaded_qdp(QOP_info_t *info,
				       QOP_FN_FermionLinksAsqtad *asqtad,
				       QOP_invert_arg_t *inv_arg,
				       QOP_resid_arg_t *res_arg,
				       QOP_F_Real mass,
				       QDP_FN_ColorVector *out,
				       QDP_FN_ColorVector *in,
				       int nthreads);

void QOP_FN_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_FN_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_F_Real *masses[],
				    int nmass[],
				    QDP_FN_ColorVector **out[],
				    QDP_FN_ColorVector *in[],
				    int nsrc);

void QOP_FN_asqtad_get_eigcg(QOP_FN_FermionLinksAsqtad *asqtad,
			     QDP_FN_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

  /* fermion force routines */

void QOP_FN_asqtad_deriv(QOP_info_t *info, QDP_FN_ColorMatrix *gauge[],
			 QDP_FN_ColorMatrix *force[],
			 QOP_asqtad_coeffs_t *coef,
			 QDP_FN_ColorMatrix *mid_fat[],
			 QDP_FN_ColorMatrix *mid_naik[]);

void QOP_FN_asqtad_force(QOP_info_t *info,
			 QOP_FN_GaugeField *gauge,
			 QOP_FN_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 QOP_F_Real eps,
			 QOP_FN_ColorVector *in_pt);

void QOP_FN_asqtad_force_multi(QOP_info_t *info,
			       QOP_FN_GaugeField *gauge,
			       QOP_FN_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       QOP_F_Real eps[],
			       QOP_FN_ColorVector *in_pt[],
			       int nsrc);

void QOP_FN_asqtad_force_multi_qdp(QOP_info_t *info,
				   QDP_FN_ColorMatrix *links[],
				   QDP_FN_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_F_Real eps[],
				   QDP_FN_ColorVector *in_pt[],
				   int nsrc);

void QOP_FN_asqtad_deriv_multi_qdp(QOP_info_t *info,
				   QDP_FN_ColorMatrix *links[],
				   QDP_FN_ColorMatrix *force[],
				   QOP_asqtad_coeffs_t *coef,
				   QOP_F_Real eps[],
				   QDP_FN_ColorVector *in_pt[],
				   int nsrc);

  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* single precision */

QOP_FN_FermionLinksHisq *
  QOP_FN_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_FN_GaugeField *gauge);

void QOP_FN_hisq_destroy_L(QOP_FN_FermionLinksHisq *field);

QOP_FN_FermionLinksAsqtad **
  QOP_FN_get_asqtad_links_from_hisq(QOP_FN_FermionLinksHisq *hl);
  
QOP_FN_FermionLinksAsqtad *
  QOP_FN_get_asqtad_deps_links_from_hisq(QOP_FN_FermionLinksHisq *hl);

  /* fermion force routines */

void QOP_FN_hisq_force_multi(QOP_info_t *info,
			     QOP_FN_FermionLinksHisq *flh,
			     QOP_FN_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     QOP_F_Real eps[],
			     QOP_FN_ColorVector *in_pt[],
			     int *n_orders_naik);

void QOP_FN_hisq_deriv_multi_qdp(QOP_info_t *info,
				 QOP_FN_FermionLinksHisq *flh,
				 QDP_FN_ColorMatrix *deriv[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_F_Real eps[],
				 QDP_FN_ColorVector *in_pt[],
				 int *n_orders_naik,
				 int doLastScale);

void QOP_FN_hisq_force_multi_qdp(QOP_info_t *info,
				 QOP_FN_FermionLinksHisq *flh,
				 QDP_FN_ColorMatrix *force[],
				 QOP_hisq_coeffs_t *coef,
				 QOP_F_Real eps[],
				 QDP_FN_ColorVector *in_pt[],
				 int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_create_L_from_raw(int nc, QDP_Lattice *lat, 
                                  QOP_F_Real *links[], QOP_F_Real *clov,
				  QOP_evenodd_t evenodd);

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_FN_GaugeField *gauge);

void QOP_FN_wilson_extract_L_to_raw(QOP_F_Real *links[], QOP_F_Real *clov,
				    QOP_FN_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_FN_wilson_destroy_L(QOP_FN_FermionLinksWilson *field);

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_convert_L_from_raw(int nc, QOP_F_Real *links[],
				   QOP_F_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_FN_wilson_convert_L_to_raw(QOP_F_Real ***links, QOP_F_Real **clov,
				    QOP_FN_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_FN_GaugeField *gauge);

QOP_FN_GaugeField *
  QOP_FN_wilson_convert_L_to_G(QOP_FN_FermionLinksWilson *links);

void QOP_FN_wilson_load_L_from_raw(QOP_FN_FermionLinksWilson *wilson,
				   QOP_F_Real *links[], QOP_F_Real *clov,
				   QOP_evenodd_t evenodd);

void QOP_FN_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_FN_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_FN_GaugeField *gauge);

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_create_L_from_qdp(QDP_FN_ColorMatrix *links[],
				  QDP_FN_DiracPropagator *clov);

void QOP_FN_wilson_extract_L_to_qdp(QDP_FN_ColorMatrix *links[],
				    QDP_FN_DiracPropagator *clov,
				    QOP_FN_FermionLinksWilson *src);

QOP_FN_FermionLinksWilson *
  QOP_FN_wilson_convert_L_from_qdp(QDP_FN_ColorMatrix *links[],
				   QDP_FN_DiracPropagator *clov);

void QOP_FN_wilson_convert_L_to_qdp(QDP_FN_ColorMatrix ***links,
				    QDP_FN_DiracPropagator **clov,
				    QOP_FN_FermionLinksWilson *src);

void QOP_FN_wilson_load_L_from_qdp(QOP_FN_FermionLinksWilson *wilson,
				   QDP_FN_ColorMatrix *links[],
				   QDP_FN_DiracPropagator *clov);


  /* inverter routines */

void QOP_FN_wilson_dslash(QOP_info_t *info,
			  QOP_FN_FermionLinksWilson *flw,
			  QOP_F_Real kappa,
			  int sign,
			  QOP_FN_DiracFermion *out,
			  QOP_FN_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_FN_wilson_diaginv(QOP_info_t *info,
			   QOP_FN_FermionLinksWilson *flw,
			   QOP_F_Real kappa,
			   QOP_FN_DiracFermion *out,
			   QOP_FN_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_FN_wilson_invert(QOP_info_t *info,
			  QOP_FN_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real kappa,
			  QOP_FN_DiracFermion *out_pt,
			  QOP_FN_DiracFermion *in_pt);

void QOP_FN_wilson_invert_multi(QOP_info_t *info,
				QOP_FN_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *kappas[],
				int nkappa[],
				QOP_FN_DiracFermion **out_pt[],
				QOP_FN_DiracFermion *in_pt[],
				int nsrc);

void QOP_FN_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_FN_FermionLinksWilson *flw,
			      QOP_F_Real kappa,
			      int sign,
			      QDP_FN_DiracFermion *out,
			      QDP_FN_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_FN_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_FN_FermionLinksWilson *flw,
			       QOP_F_Real kappa,
			       QDP_FN_DiracFermion *out,
			       QDP_FN_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_FN_wilson_invert_qdp(QOP_info_t *info,
			      QOP_FN_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      QOP_F_Real kappa,
			      QDP_FN_DiracFermion *out_pt,
			      QDP_FN_DiracFermion *in_pt);

void QOP_FN_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_FN_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    QOP_F_Real *kappas[],
				    int nkappa[],
				    QDP_FN_DiracFermion **out_pt[],
				    QDP_FN_DiracFermion *in_pt[],
				    int nsrc);

void QOP_FN_wilson_invert_ne_qdp(QOP_info_t *info,
				 QOP_FN_FermionLinksWilson *flw,
				 QOP_invert_arg_t *inv_arg,
				 QOP_resid_arg_t *res_arg,
				 QOP_F_Real kappa,
				 QDP_FN_DiracFermion *out,
				 QDP_FN_DiracFermion *in);

  /* fermion force routines */

void QOP_FN_wilson_deriv_multi_qdp(QOP_info_t *info,
				   QOP_FN_FermionLinksWilson *flw,
				   QDP_FN_ColorMatrix *deriv[],
				   QOP_F_Real eps[],
				   QDP_FN_DiracFermion *x[],
				   QDP_FN_DiracFermion *y[],
				   int n);

void QOP_FN_wilson_force_multi_qdp(QOP_info_t *info,
				   QOP_FN_FermionLinksWilson *flw,
				   QDP_FN_ColorMatrix *force[],
				   QOP_F_Real eps[],
				   QDP_FN_DiracFermion *x[],
				   QDP_FN_DiracFermion *y[],
				   int n);

void QOP_FN_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
					QOP_FN_FermionLinksWilson *flw,
					QDP_FN_ColorMatrix *deriv[],
					QOP_F_Real kappa[],
					QOP_F_Real eps[],
					QDP_FN_DiracFermion *x[],
					QDP_FN_DiracFermion *y[],
					int n);

void QOP_FN_wilson_force_prec_multi_qdp(QOP_info_t *info,
					QOP_FN_FermionLinksWilson *flw,
					QDP_FN_ColorMatrix *force[],
					QOP_F_Real kappa[],
					QOP_F_Real eps[],
					QDP_FN_DiracFermion *x[],
					QDP_FN_DiracFermion *y[],
					int n);

  // new fermilab action IFLA -- added by bugra --------------- :

void QOP_FN_wilson_ifla_dslash(QOP_info_t *info,
			       QOP_FN_FermionLinksWilson *flw,
			       QOP_F_Real kappa,
			       int sign,
			       QOP_wilson_ifla_coeffs_t *coeffs,
			       QOP_FN_DiracFermion *out,
			       QOP_FN_DiracFermion *in,
			       QOP_evenodd_t eo_out,
			       QOP_evenodd_t eo_in);

void QOP_FN_wilson_ifla_dslash_qdp(QOP_info_t *info,
				   QOP_FN_FermionLinksWilson *flw,
				   QOP_F_Real kappa,
				   int sign,
				   QOP_wilson_ifla_coeffs_t *coeffs,
				   QDP_FN_DiracFermion *out,
				   QDP_FN_DiracFermion *in,
				   QOP_evenodd_t eo_out,
				   QOP_evenodd_t eo_in);

  // MULTIGRID STUFF

#ifndef _QOP_N_MG
#define _QOP_N_MG
typedef struct QOP_N_WilsonMgStruct QOP_N_WilsonMg;
QOP_N_WilsonMg *QOP_N_wilsonMgNew(void);
void QOP_N_wilsonMgFree(QOP_N_WilsonMg *wmg);
void QOP_N_wilsonMgSet(QOP_N_WilsonMg *wmg, int l, char *s, double val);
void QOP_N_wilsonMgSetArray(QOP_N_WilsonMg *wmg, int l, char *s, double *vals, int nval);
void QOP_N_wilsonMgSetup(QOP_N_WilsonMg *wmg);
#endif // _QOP_N_MG

void QOP_FN_wilsonMgSetLinks(QOP_N_WilsonMg *wmg, QOP_FN_FermionLinksWilson *wil);
void QOP_FN_wilsonMgSolve(QOP_info_t *info, QOP_N_WilsonMg *wmg,
			  QOP_FN_FermionLinksWilson *flw,
			  QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg,
			  QLA_Real kappa, QDP_DiracFermion *out, QDP_DiracFermion *in);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_FN_FermionLinksDW *
  QOP_FN_dw_create_L_from_raw(int nc, QDP_Lattice *lat,
                              QOP_F_Real *links[],
			      QOP_evenodd_t evenodd);

QOP_FN_FermionLinksDW *
  QOP_FN_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_FN_GaugeField *gauge);

void QOP_FN_dw_extract_L_to_raw(QOP_F_Real *links[],
				QOP_FN_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_FN_dw_destroy_L(QOP_FN_FermionLinksDW *field);

QOP_FN_FermionLinksDW *
  QOP_FN_dw_convert_L_from_raw(int nc, QOP_F_Real *links[],
			       QOP_evenodd_t evenodd);

void QOP_FN_dw_convert_L_to_raw(QOP_F_Real ***links,
				QOP_FN_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_FN_FermionLinksDW *
  QOP_FN_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_FN_GaugeField *gauge);

QOP_FN_GaugeField *
  QOP_FN_dw_convert_L_to_G(QOP_FN_FermionLinksDW *links);

void QOP_FN_dw_load_L_from_raw(QOP_FN_FermionLinksDW *dw,
			       QOP_F_Real *links[], QOP_evenodd_t evenodd);

void QOP_FN_dw_load_L_from_G(QOP_info_t *info,
			     QOP_FN_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_FN_GaugeField *gauge);

QOP_FN_FermionLinksDW *
  QOP_FN_dw_create_L_from_qdp(QDP_FN_ColorMatrix *links[]);

void QOP_FN_dw_extract_L_to_qdp(QDP_FN_ColorMatrix *links[],
				QOP_FN_FermionLinksDW *src);

QOP_FN_FermionLinksDW *
  QOP_FN_dw_convert_L_from_qdp(QDP_FN_ColorMatrix *links[]);

void QOP_FN_dw_convert_L_to_qdp(QDP_FN_ColorMatrix ***links,
				QOP_FN_FermionLinksDW *src);

void QOP_FN_dw_load_L_from_qdp(QOP_FN_FermionLinksDW *dw,
			       QDP_FN_ColorMatrix *links[]);

  /* inverter routines */

void QOP_FN_dw_dslash(QOP_info_t *info,
		      QOP_FN_FermionLinksDW *links,
		      QOP_F_Real M5,
		      QOP_F_Real m,
		      int sign,
		      QOP_FN_DiracFermion *out_pt[],
		      QOP_FN_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_FN_dw_dslash2(QOP_info_t *info,
		       QOP_FN_FermionLinksDW *links,
		       QOP_F_Real M5,
		       QOP_F_Real m,
		       QOP_FN_DiracFermion *out_pt[],
		       QOP_FN_DiracFermion *in_pt[],
		       int Ls,
		       QOP_evenodd_t eo_out,
		       QOP_evenodd_t eo_in);

void QOP_FN_dw_invert(QOP_info_t *info,
		      QOP_FN_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QOP_F_Real M5,
		      QOP_F_Real m,
		      QOP_FN_DiracFermion *out_pt[],
		      QOP_FN_DiracFermion *in_pt[],
		      int Ls);

void QOP_FN_dw_invert_multi(QOP_info_t *info,
			    QOP_FN_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_F_Real *M5[],
			    QOP_F_Real *m[],
			    int nmass[],
			    QOP_FN_DiracFermion ***out_pt[],
			    QOP_FN_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_FN_dw_dslash_qdp(QOP_info_t *info,
			  QOP_FN_FermionLinksDW *links,
			  QOP_F_Real M5, 
			  QOP_F_Real m,
			  int sign,
			  QDP_FN_DiracFermion *out_pt[],
			  QDP_FN_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_FN_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_FN_FermionLinksDW *links,
			   QOP_F_Real M5,
			   QOP_F_Real m,
			   QDP_FN_DiracFermion *out_pt[],
			   QDP_FN_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);

void QOP_FN_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_FN_FermionLinksDW *fldw,
			   QOP_F_Real M5,
			   QOP_F_Real m,
			   QDP_FN_DiracFermion **out,
			   QDP_FN_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_FN_dw_invert_qdp(QOP_info_t *info,
			  QOP_FN_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QOP_F_Real M5,
			  QOP_F_Real m,
			  QDP_FN_DiracFermion *out[],
			  QDP_FN_DiracFermion *in[],
			  int Ls);

void QOP_FN_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_FN_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				QOP_F_Real *M5[],
				QOP_F_Real *m[],
				int nmass[],
				QDP_FN_DiracFermion ***out[],
				QDP_FN_DiracFermion **in[],
				int nsrc,
				int Ls);

 /* fermion force routines */

void QOP_FN_dw_force(QOP_info_t *info,
		     QOP_FN_GaugeField *gauge,
		     QOP_FN_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     QOP_F_Real eps,
		     QOP_FN_DiracFermion *in_pt);

void QOP_FN_dw_force_multi(QOP_info_t *info,
			   QOP_FN_GaugeField *gauge,
			   QOP_FN_Force *force,
			   QOP_dw_coeffs_t *coef,
			   QOP_F_Real eps[],
			   QOP_FN_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == 'F'
#  if QOP_Colors == 'N'
#    include <qop_fn_generic.h>
#  endif
#  include <qop_fn_precision_generic.h>
#endif
#if QOP_Colors == 'N'
#  include <qop_fn_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_FN_H */
