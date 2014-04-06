// DO NOT EDIT
// generated from qop_pc_internal.h
#ifndef _QOP_DN_INTERNAL_H
#define _QOP_DN_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_DN_ColorVector_struct {
  QDP_DN_ColorVector *cv;
  QOP_D_Real *raw;
};

struct QOP_DN_DiracFermion_struct {
  QDP_DN_DiracFermion *df;
  QOP_D_Real *raw;
};

typedef void (*QOP_DN_gauge_deriv)(QDP_DN_ColorMatrix **d[],
				    QOP_DN_GaugeField *g,
				    QDP_DN_ColorMatrix *c[]);

typedef void (*QOP_DN_gauge_scale)(QDP_DN_ColorMatrix *l[],
				    QOP_DN_GaugeField *g, int inv);

struct QOP_DN_GaugeField_struct {
  QDP_DN_ColorMatrix **links;
  QOP_D_Real **raw;
  QOP_DN_GaugeField **parents;
  QOP_DN_gauge_deriv deriv;
  QOP_DN_gauge_scale scale;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t sign;
  int chained;
  int nparents;
};

struct QOP_DN_Force_struct {
  QDP_DN_ColorMatrix **force;
  QOP_D_Real **raw;
};

typedef struct {
  QDP_DN_ColorVector **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_DN_eigcg_t_V;

typedef struct {
  QDP_DN_DiracFermion **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_DN_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_DN_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_DN_ColorMatrix **fatlinks;
  QDP_DN_ColorMatrix **longlinks;
  QDP_DN_ColorMatrix **fwdlinks;
  QDP_DN_ColorMatrix **bcklinks;
  QDP_DN_ColorMatrix **dbllinks;
  QOP_DN_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
  //QOP_FN_FermionLinksAsqtad *ofla;
};

  /* HISQ datatypes */

struct QOP_DN_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_DN_ColorMatrix **U_links; // gauge links
  QDP_DN_ColorMatrix **V_links; // Fat7 smeared
  QDP_DN_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_DN_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_DN_FermionLinksAsqtad **fn;
  QOP_DN_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_DN_FermionLinksWilson_struct {
  QOP_D_Real clovinvkappa;
  int dblstored;
  QDP_DN_ColorMatrix **links;
  QDP_DN_ColorMatrix **bcklinks;
  QDP_DN_ColorMatrix **dbllinks;
  QOP_DN_GaugeField *qopgf;
  QOP_DN_GaugeField *gauge;
  QDP_DN_DiracPropagator *qdpclov;
  QOP_D_Real *clov, *clovinv;
  QOP_D_Real **rawlinks, *rawclov;
  QOP_DN_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_DN_FermionLinksDW_struct {
  QOP_DN_FermionLinksWilson *flw;
};

/* internal routines */

QOP_DN_FermionLinksAsqtad *QOP_DN_asqtad_create_L_from_L(QOP_DN_FermionLinksAsqtad *fla_src);
QOP_DN_FermionLinksAsqtad *QOP_DN_asqtad_create_L_from_r_times_L(QOP_D_Real s,
								  QOP_DN_FermionLinksAsqtad *fla_src);
void QOP_DN_asqtad_L_peq_L(QOP_DN_FermionLinksAsqtad *fla, QOP_DN_FermionLinksAsqtad *fla1);
void QOP_DN_qdpM_eq_raw(QDP_DN_ColorMatrix *cm, QOP_D_Real *lnk);

typedef void (QOP_DN_linop_t_V)(QDP_DN_ColorVector *out, QDP_DN_ColorVector *in, QDP_Subset subset);
typedef void (QOP_DN_linop_t_D)(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_DN_linop_t_vD)(QDP_DN_DiracFermion **out, QDP_DN_DiracFermion **in, QDP_Subset subset);

typedef QOP_D_Real (QOP_DN_linopn_t_V)(QDP_DN_ColorVector *out, QDP_DN_ColorVector *in, QDP_Subset subset);
typedef QOP_D_Real (QOP_DN_linopn_t_D)(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in, QDP_Subset subset);
typedef QOP_D_Real (QOP_DN_linopn_t_vD)(QDP_DN_DiracFermion **out, QDP_DN_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_DN_invert_cg_V(QOP_DN_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_DN_ColorVector *out,
		    QDP_DN_ColorVector *in,
		    QDP_DN_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_DN_invert_cg_D(QOP_DN_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_DN_DiracFermion *out,
		    QDP_DN_DiracFermion *in,
		    QDP_DN_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_DN_invert_cg_vD(QOP_DN_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_DN_DiracFermion **out,
		     QDP_DN_DiracFermion **in,
		     QDP_DN_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_DN_invert_cgms_V(QOP_DN_linopn_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_DN_ColorVector **out,
		      QDP_DN_ColorVector *in,
		      QDP_DN_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_DN_invert_cgms_D(QOP_DN_linopn_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_DN_DiracFermion **out,
		      QDP_DN_DiracFermion *in,
		      QDP_DN_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_DN_invert_cgms_vD(QOP_DN_linopn_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_D_Real *shifts,
		       int nshifts,
		       QDP_DN_DiracFermion ***out,
		       QDP_DN_DiracFermion **in,
		       QDP_DN_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_DN_invert_bicgstab_D(QOP_DN_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_DN_DiracFermion *out,
			  QDP_DN_DiracFermion *in,
			  QDP_DN_DiracFermion *p,
			  QDP_DN_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_DN_invert_eigcg_V(QOP_DN_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_DN_ColorVector *out,
		       QDP_DN_ColorVector *in,
		       QDP_DN_ColorVector *p,
		       QDP_Subset subset,
		       QOP_DN_eigcg_t_V *eigcg);

QOP_status_t
QOP_DN_invert_eigcg_D(QOP_DN_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_DN_DiracFermion *out,
		       QDP_DN_DiracFermion *in,
		       QDP_DN_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_DN_eigcg_t_D *eigcg);


QDP_DN_ColorVector *QOP_DN_asqtad_dslash_get_tmp(QOP_DN_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_DN_DiracFermion *QOP_DN_wilson_dslash_get_tmp(QOP_DN_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);
QOP_DN_FermionLinksWilson *QOP_DN_wilson_initialize_gauge_L(void);

void QOP_DN_get_mid(QOP_info_t *info, QDP_DN_ColorMatrix *mid[],
		     QDP_Shift shifts[], int ns, QOP_D_Real eps[],
		     QOP_D_Real scale, QDP_DN_ColorVector *x[], int nterms);

void QOP_DN_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QDP_DN_ColorMatrix *links[],
					  QDP_DN_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[], QDP_DN_ColorVector *x[], int nsrc);

void QOP_DN_asqtad_deriv_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_DN_ColorMatrix *links[],
					  QDP_DN_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[],
					  QDP_DN_ColorVector *x[],
					  int nterms);

void QOP_DN_asqtad_force_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_DN_ColorMatrix *links[],
					  QDP_DN_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[],
					  QDP_DN_ColorVector *x[],
					  int nterms);

//AB internal operations for HISQ

void 
QOP_DN_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_DN_FermionLinksHisq *flh,
				       QOP_DN_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_D_Real *epsv,
				       QDP_DN_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_DN_hisq_deriv_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_DN_FermionLinksHisq *flh,
				    QDP_DN_ColorMatrix *deriv[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_D_Real *epsv,
				    QDP_DN_ColorVector *in_pt[], 
				    int *n_orders_naik);

void 
QOP_DN_hisq_force_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_DN_FermionLinksHisq *flh,
				    QDP_DN_ColorMatrix *force[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_D_Real *epsv,
				    QDP_DN_ColorVector *in_pt[], 
				    int *n_orders_naik);

void QOP_DN_dw_schur2_qdp(QOP_info_t *info, QOP_DN_FermionLinksDW *fldw,
			   QOP_D_Real M5, QOP_D_Real mq,
			   QDP_DN_DiracFermion *out[], QDP_DN_DiracFermion *in[],
			   int ls,QOP_evenodd_t eo);
void QOP_DN_dw_schur_qdp(QOP_info_t *info, QOP_DN_FermionLinksDW *fldw,
			  QOP_D_Real M5, QOP_D_Real mq, int sign,
			  QDP_DN_DiracFermion *out[], QDP_DN_DiracFermion *in[],
			  int ls, QOP_evenodd_t eo);
void
QOP_DN_dw_EO_project(QOP_DN_FermionLinksDW *fldw,
		      QDP_DN_DiracFermion *out[], QDP_DN_DiracFermion *in[],
		      QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);
void
QOP_DN_dw_EO_reconstruct(QOP_DN_FermionLinksDW *fldw,
			  QDP_DN_DiracFermion *out[], QDP_DN_DiracFermion *in[],
			  QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_DN_invert_gcr2_D(QOP_DN_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_DN_DiracFermion *out,
		      QDP_DN_DiracFermion *in,
		      QDP_DN_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_DN_invert_gmres2_D(QOP_DN_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_DN_DiracFermion *out,
			QDP_DN_DiracFermion *in,
			QDP_DN_DiracFermion *r,
			QDP_Subset subset);

QOP_D_Real
QOP_DN_relnorm2_V(QDP_DN_ColorVector **rsd, 
		   QDP_DN_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_D_Real
QOP_DN_relnorm2_D(QDP_DN_DiracFermion **rsd, 
		   QDP_DN_DiracFermion **out, 
		   QDP_Subset subset, int nv);

//// MULTIGRID STUFF

typedef struct {
  QOP_DN_FermionLinksWilson *wil;
  QLA_D_Real kappa;
} QOP_DN_WilArgs;

void QOP_DN_wilsonDslash(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in,
			  QOP_DN_FermionLinksWilson *wil, QLA_D_Real kappa,
			  int sign, QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_DN_wilsonDiaginv(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in,
			   QOP_DN_FermionLinksWilson *wil, QLA_D_Real kappa,
			   QOP_evenodd_t pout);
void QOP_DN_wilsonDslashEO(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in,
			    QOP_DN_FermionLinksWilson *wil, QLA_D_Real kappa,
			    int sign, QOP_evenodd_t par);
void QOP_DN_wilsonDslashEOS(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_D_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_DN_wilsonDslashEOH(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_D_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_DN_wilEoProjectD(QDP_DN_DiracFermion *ineo, QDP_DN_DiracFermion *in,
			   QOP_DN_WilArgs *w);
void QOP_DN_wilEoReconstructD(QDP_DN_DiracFermion *out, QDP_DN_DiracFermion *outeo,
			       QDP_DN_DiracFermion *in, QOP_DN_WilArgs *w);

#ifdef HAVE_NCN
#include <qdp_fn.h>
#include <qdp_dn.h>

void QOP_DN_V1eqD(QDP_DN_ColorVector *v[1], QDP_DN_DiracFermion *d, QDP_Subset sub);
void QOP_DN_DeqV1(QDP_DN_DiracFermion *d, QDP_DN_ColorVector *v[1], QDP_Subset sub);
void QOP_DN_V2eqD(QDP_DN_ColorVector *v[2], QDP_DN_DiracFermion *d, QDP_Subset sub);
void QOP_DN_DeqV2(QDP_DN_DiracFermion *d, QDP_DN_ColorVector *v[2], QDP_Subset sub);
void QOP_DN_wilDV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_DN_wilDV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_DN_wilPV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_DN_wilPV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_DN_wilPNEV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_DN_wilEoV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_DN_wilEoV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_DN_wilEoProjectV1(QDP_DN_ColorVector *ineo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_DN_wilEoReconstructV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_DN_wilEoReconstructPV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_DN_wilEoProjectV2(QDP_DN_ColorVector *ineo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_DN_wilEoReconstructV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_DN_wilEoReconstructPV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);

#endif // HAVE_NCN

#ifndef _QOP_N_MG_INTERNAL
#define _QOP_N_MG_INTERNAL

#include <qop_f_internal.h>
#include <qop_d_internal.h>
//#include <qop_mg_internal.h>
//struct QOP_WilMgLevel;

struct QOP_N_WilsonMgStruct {
  QOP_FN_FermionLinksWilson *wilF;
  QOP_DN_FermionLinksWilson *wilD;
  QLA_F_Real kappa;
  QLA_F_Real kappanv;
  QOP_FN_WilArgs vcwaF;
  QOP_FN_WilArgs nvwaF;
  int nlevels;
  struct QOP_WilMgLevel *mg;
  int verbose;
  int profile;
  int itmax;
  QOP_F_Gcr *gcrF;
  QOP_D_Gcr *gcrD;
  int ngcr;
  int nc;
};

#endif // _QOP_N_MG_INTERNAL

#if QOP_Precision == 'D'
#  if QOP_Colors == 'N'
#    include <qop_dn_internal_generic.h>
#  endif
#endif

#endif /* _QOP_DN_INTERNAL_H */
