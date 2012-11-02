// DO NOT EDIT
// generated from qop_pc_internal.h
#ifndef _QOP_F2_INTERNAL_H
#define _QOP_F2_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_F2_ColorVector_struct {
  QDP_F2_ColorVector *cv;
  QOP_F_Real *raw;
};

struct QOP_F2_DiracFermion_struct {
  QDP_F2_DiracFermion *df;
  QOP_F_Real *raw;
};

typedef void (*QOP_F2_gauge_deriv)(QDP_F2_ColorMatrix **d[],
				    QOP_F2_GaugeField *g,
				    QDP_F2_ColorMatrix *c[]);

typedef void (*QOP_F2_gauge_scale)(QDP_F2_ColorMatrix *l[],
				    QOP_F2_GaugeField *g, int inv);

struct QOP_F2_GaugeField_struct {
  QDP_F2_ColorMatrix **links;
  QOP_F_Real **raw;
  QOP_F2_GaugeField **parents;
  QOP_F2_gauge_deriv deriv;
  QOP_F2_gauge_scale scale;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t sign;
  int chained;
  int nparents;
};

struct QOP_F2_Force_struct {
  QDP_F2_ColorMatrix **force;
  QOP_F_Real **raw;
};

typedef struct {
  QDP_F2_ColorVector **u;
  QOP_F_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_F2_eigcg_t_V;

typedef struct {
  QDP_F2_DiracFermion **u;
  QOP_F_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_F2_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_F2_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_F2_ColorMatrix **fatlinks;
  QDP_F2_ColorMatrix **longlinks;
  QDP_F2_ColorMatrix **fwdlinks;
  QDP_F2_ColorMatrix **bcklinks;
  QDP_F2_ColorMatrix **dbllinks;
  QOP_F2_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
};

  /* HISQ datatypes */

struct QOP_F2_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_F2_ColorMatrix **U_links; // gauge links
  QDP_F2_ColorMatrix **V_links; // Fat7 smeared
  QDP_F2_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_F2_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_F2_FermionLinksAsqtad **fn;
  QOP_F2_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_F2_FermionLinksWilson_struct {
  QOP_F_Real clovinvkappa;
  int dblstored;
  QDP_F2_ColorMatrix **links;
  QDP_F2_ColorMatrix **bcklinks;
  QDP_F2_ColorMatrix **dbllinks;
  QOP_F2_GaugeField *qopgf;
  QOP_F2_GaugeField *gauge;
  QDP_F2_DiracPropagator *qdpclov;
  QOP_F_Real *clov, *clovinv;
  QOP_F_Real **rawlinks, *rawclov;
  QOP_F2_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_F2_FermionLinksDW_struct {
  QOP_F2_FermionLinksWilson *flw;
};

/* internal routines */

QOP_F2_FermionLinksAsqtad *QOP_F2_asqtad_create_L_from_L(QOP_F2_FermionLinksAsqtad *fla_src);
QOP_F2_FermionLinksAsqtad *QOP_F2_asqtad_create_L_from_r_times_L(QOP_F_Real *s,
								  QOP_F2_FermionLinksAsqtad *fla_src);
void QOP_F2_asqtad_L_peq_L(QOP_F2_FermionLinksAsqtad *fla, QOP_F2_FermionLinksAsqtad *fla1);
void QOP_F2_qdpM_eq_raw(QDP_F2_ColorMatrix *cm, QOP_F_Real *lnk);
typedef void (QOP_F2_linop_t_V)(QDP_F2_ColorVector *out, QDP_F2_ColorVector *in, QDP_Subset subset);
typedef void (QOP_F2_linop_t_D)(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_F2_linop_t_vD)(QDP_F2_DiracFermion **out, QDP_F2_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cg_V(QOP_F2_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_F2_ColorVector *out,
		    QDP_F2_ColorVector *in,
		    QDP_F2_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cg_D(QOP_F2_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_F2_DiracFermion *out,
		    QDP_F2_DiracFermion *in,
		    QDP_F2_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cg_vD(QOP_F2_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_F2_DiracFermion **out,
		     QDP_F2_DiracFermion **in,
		     QDP_F2_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_F2_invert_cgms_V(QOP_F2_linop_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_F2_ColorVector **out,
		      QDP_F2_ColorVector *in,
		      QDP_F2_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cgms_D(QOP_F2_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_F2_DiracFermion **out,
		      QDP_F2_DiracFermion *in,
		      QDP_F2_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cgms_vD(QOP_F2_linop_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_F_Real *shifts,
		       int nshifts,
		       QDP_F2_DiracFermion ***out,
		       QDP_F2_DiracFermion **in,
		       QDP_F2_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_F2_invert_bicgstab_D(QOP_F2_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_F2_DiracFermion *out,
			  QDP_F2_DiracFermion *in,
			  QDP_F2_DiracFermion *p,
			  QDP_F2_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_F2_invert_eigcg_V(QOP_F2_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_F2_ColorVector *out,
		       QDP_F2_ColorVector *in,
		       QDP_F2_ColorVector *p,
		       QDP_Subset subset,
		       QOP_F2_eigcg_t_V *eigcg);

QOP_status_t
QOP_F2_invert_eigcg_D(QOP_F2_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_F2_DiracFermion *out,
		       QDP_F2_DiracFermion *in,
		       QDP_F2_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_F2_eigcg_t_D *eigcg);


QDP_F2_ColorVector *QOP_F2_asqtad_dslash_get_tmp(QOP_F2_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_F2_DiracFermion *QOP_F2_wilson_dslash_get_tmp(QOP_F2_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);

void QOP_F2_get_mid(QOP_info_t *info, QDP_F2_ColorMatrix *mid[], QDP_Shift shifts[], int ns,
		     QOP_F_Real eps[], QDP_F2_ColorVector *x[], int nterms);

void QOP_F2_asqtad_deriv(QOP_info_t *info, QDP_F2_ColorMatrix *gauge[],
			  QDP_F2_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
			  QDP_F2_ColorMatrix *mid_fat[], QDP_F2_ColorMatrix *mid_naik[]);

void QOP_F2_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QOP_F2_GaugeField *gauge,
					  QDP_F2_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[], QDP_F2_ColorVector *x[], int nsrc);

void QOP_F2_asqtad_force_multi_fnmat_qdp(QOP_info_t *info, QOP_F2_GaugeField *gauge,
					  QDP_F2_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[], QDP_F2_ColorVector *x[], int nterms);

//AB internal operations for HISQ

void
QOP_F2_hisq_force_multi_reunit(QOP_info_t *info,
				QDP_F2_ColorMatrix *gf[4],
				QDP_F2_ColorMatrix *force_accum[4],
				QDP_F2_ColorMatrix *force_accum_old[4]);

void 
QOP_F2_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_F2_FermionLinksHisq *flh,
				       QOP_F2_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_F_Real *epsv,
				       QDP_F2_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_F2_hisq_deriv_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_F2_FermionLinksHisq *flh,
				    QDP_F2_ColorMatrix *deriv[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_F_Real *epsv,
				    QDP_F2_ColorVector *in_pt[], 
				    int *n_orders_naik);

void 
QOP_F2_hisq_force_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_F2_FermionLinksHisq *flh,
				    QDP_F2_ColorMatrix *force[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_F_Real *epsv,
				    QDP_F2_ColorVector *in_pt[], 
				    int *n_orders_naik);

void QOP_F2_u3reunit(QOP_info_t *info, QDP_F2_ColorMatrix *U, QDP_F2_ColorMatrix *V);

void QOP_F2_su3reunit(QOP_info_t *info, QDP_F2_ColorMatrix *U, QDP_F2_ColorMatrix *Ur);

void
QOP_F2_dw_schur2_qdp(QOP_info_t *info, QOP_F2_FermionLinksDW *fldw,
		      QOP_F_Real M5, QOP_F_Real mq,
		      QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
		      int ls, QOP_evenodd_t eo);
void
QOP_F2_dw_schur_qdp(QOP_info_t *info, QOP_F2_FermionLinksDW *fldw,
		     QOP_F_Real M5, QOP_F_Real mq, int sign,
		     QDP_DiracFermion *out[], QDP_F2_DiracFermion *in[],
		     int ls, QOP_evenodd_t eo);
extern void
QOP_F2_dw_EO_project(QOP_F2_FermionLinksDW *fldw,
		      QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
		      QOP_F_Real M5, QOP_F_Real mq, int ls, QOP_evenodd_t eo);
extern void
QOP_F2_dw_EO_reconstruct(QOP_F2_FermionLinksDW *fldw,
			  QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
			  QOP_F_Real M5, QOP_F_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_F2_invert_gcr2_D(QOP_F2_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_F2_DiracFermion *out,
		      QDP_F2_DiracFermion *in,
		      QDP_F2_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_F2_invert_gmres2_D(QOP_F2_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_F2_DiracFermion *out,
			QDP_F2_DiracFermion *in,
			QDP_F2_DiracFermion *r,
			QDP_Subset subset);

QOP_F_Real
QOP_F2_relnorm2_V(QDP_F2_ColorVector **rsd, 
		   QDP_F2_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_F_Real
QOP_F2_relnorm2_D(QDP_F2_DiracFermion **rsd, 
		   QDP_F2_DiracFermion **out, 
		   QDP_Subset subset, int nv);

#if QOP_Precision == 'F'
#  if QOP_Colors == 2
#    include <qop_f2_internal_generic.h>
#  endif
#endif

#endif /* _QOP_F2_INTERNAL_H */
