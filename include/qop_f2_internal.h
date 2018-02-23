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
  //QOP_D2_FermionLinksAsqtad *ofla;
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
QOP_F2_FermionLinksAsqtad *QOP_F2_asqtad_create_L_from_r_times_L(QOP_D_Real s,
								  QOP_F2_FermionLinksAsqtad *fla_src);
void QOP_F2_asqtad_L_peq_L(QOP_F2_FermionLinksAsqtad *fla, QOP_F2_FermionLinksAsqtad *fla1);
void QOP_F2_qdpM_eq_raw(QDP_F2_ColorMatrix *cm, QOP_F_Real *lnk);

typedef void (QOP_F2_linop_t_V)(QDP_F2_ColorVector *out, QDP_F2_ColorVector *in, QDP_Subset subset);
typedef void (QOP_F2_linop_t_D)(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_F2_linop_t_vD)(QDP_F2_DiracFermion **out, QDP_F2_DiracFermion **in, QDP_Subset subset);

typedef QOP_F_Real (QOP_F2_linopn_t_V)(QDP_F2_ColorVector *out, QDP_F2_ColorVector *in, QDP_Subset subset);
typedef QOP_F_Real (QOP_F2_linopn_t_D)(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in, QDP_Subset subset);
typedef QOP_F_Real (QOP_F2_linopn_t_vD)(QDP_F2_DiracFermion **out, QDP_F2_DiracFermion **in, QDP_Subset subset);

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
QOP_F2_invert_cgms_V(QOP_F2_linopn_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_F2_ColorVector **out,
		      QDP_F2_ColorVector *in,
		      QDP_F2_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cgms_D(QOP_F2_linopn_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_F2_DiracFermion **out,
		      QDP_F2_DiracFermion *in,
		      QDP_F2_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_F2_invert_cgms_vD(QOP_F2_linopn_t_vD *linop,
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
QOP_F2_FermionLinksWilson *QOP_F2_wilson_initialize_gauge_L(void);

void QOP_F2_get_mid(QOP_info_t *info, QDP_F2_ColorMatrix *mid[],
		     QDP_Shift shifts[], int ns, QOP_F_Real eps[],
		     QOP_F_Real scale, QDP_F2_ColorVector *x[], int nterms);

void QOP_F2_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QDP_F2_ColorMatrix *links[],
					  QDP_F2_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[], QDP_F2_ColorVector *x[], int nsrc);

void QOP_F2_asqtad_deriv_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_F2_ColorMatrix *links[],
					  QDP_F2_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[],
					  QDP_F2_ColorVector *x[],
					  int nterms);

void QOP_F2_asqtad_force_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_F2_ColorMatrix *links[],
					  QDP_F2_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[],
					  QDP_F2_ColorVector *x[],
					  int nterms);

//AB internal operations for HISQ

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

void QOP_F2_dw_schur2_qdp(QOP_info_t *info, QOP_F2_FermionLinksDW *fldw,
			   QOP_F_Real M5, QOP_F_Real mq,
			   QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
			   int ls,QOP_evenodd_t eo);
void QOP_F2_dw_schur_qdp(QOP_info_t *info, QOP_F2_FermionLinksDW *fldw,
			  QOP_F_Real M5, QOP_F_Real mq, int sign,
			  QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
			  int ls, QOP_evenodd_t eo);
void
QOP_F2_dw_EO_project(QOP_F2_FermionLinksDW *fldw,
		      QDP_F2_DiracFermion *out[], QDP_F2_DiracFermion *in[],
		      QOP_F_Real M5, QOP_F_Real mq, int ls, QOP_evenodd_t eo);
void
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

//// MULTIGRID STUFF

typedef struct {
  QOP_F2_FermionLinksWilson *wil;
  QLA_F_Real kappa;
} QOP_F2_WilArgs;

void QOP_F2_wilsonDslash(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in,
			  QOP_F2_FermionLinksWilson *wil, QLA_F_Real kappa,
			  int sign, QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_F2_wilsonDiaginv(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in,
			   QOP_F2_FermionLinksWilson *wil, QLA_F_Real kappa,
			   QOP_evenodd_t pout);
void QOP_F2_wilsonDslashEO(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in,
			    QOP_F2_FermionLinksWilson *wil, QLA_F_Real kappa,
			    int sign, QOP_evenodd_t par);
void QOP_F2_wilsonDslashEOS(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_F_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_F2_wilsonDslashEOH(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_F_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_F2_wilEoProjectD(QDP_F2_DiracFermion *ineo, QDP_F2_DiracFermion *in,
			   QOP_F2_WilArgs *w);
void QOP_F2_wilEoReconstructD(QDP_F2_DiracFermion *out, QDP_F2_DiracFermion *outeo,
			       QDP_F2_DiracFermion *in, QOP_F2_WilArgs *w);

#ifdef HAVE_NCN
#include <qdp_fn.h>
#include <qdp_dn.h>

void QOP_F2_V1eqD(QDP_FN_ColorVector *v[1], QDP_F2_DiracFermion *d, QDP_Subset sub);
void QOP_F2_DeqV1(QDP_F2_DiracFermion *d, QDP_FN_ColorVector *v[1], QDP_Subset sub);
void QOP_F2_V2eqD(QDP_FN_ColorVector *v[2], QDP_F2_DiracFermion *d, QDP_Subset sub);
void QOP_F2_DeqV2(QDP_F2_DiracFermion *d, QDP_FN_ColorVector *v[2], QDP_Subset sub);
void QOP_F2_wilDV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F2_wilDV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F2_wilPV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F2_wilPV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F2_wilPNEV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F2_wilEoV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_F2_wilEoV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_F2_wilEoProjectV1(QDP_FN_ColorVector *ineo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F2_wilEoReconstructV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F2_wilEoReconstructPV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_F2_wilEoProjectV2(QDP_FN_ColorVector *ineo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_F2_wilEoReconstructV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_F2_wilEoReconstructPV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);

#ifndef _QOP_2_MG_INTERNAL
#define _QOP_2_MG_INTERNAL

#include <qop_f_internal.h>
#include <qop_d_internal.h>
//#include <qop_mg_internal.h>
//struct QOP_WilMgLevel;

struct QOP_2_WilsonMgStruct {
  QOP_F2_FermionLinksWilson *wilF;
  QOP_F2_FermionLinksWilson *wilF_priv;
  QOP_D2_FermionLinksWilson *wilD;
  QLA_F_Real kappa;
  QLA_F_Real kappanv;
  QOP_F2_WilArgs vcwaF;
  QOP_F2_WilArgs nvwaF;
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

#endif // _QOP_2_MG_INTERNAL

#endif // HAVE_NCN

#if QOP_Precision == 'F'
#  if QOP_Colors == 2
#    include <qop_f2_internal_generic.h>
#  endif
#endif

#ifdef HAVE_QLL

void setup_qll_solverF2(QOP_FermionLinksAsqtad *fla);
void free_qll_solverF2(void);
void solve_qllF2(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		    QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg);
void solveMulti_qllF2(QDP_ColorVector *dest[], QDP_ColorVector *src,
			 double ms[], int nm,  QOP_invert_arg_t *invarg,
			 QOP_resid_arg_t *resargs[]);
void * create_qll_gaugeF2(int nc);
void * create_qll_from_gaugeF2(QDP_ColorMatrix *g[]);
void copy_gauge_from_qllF2(QDP_ColorMatrix *g[], void *ff);
void free_qll_gaugeF2(void *ff);
void fat7_qllF2(void *qllfl, void *qllll, QOP_asqtad_coeffs_t *coef,
		   void *qllu, void *qllul);

#endif // HAVE_QLL

#ifdef HAVE_QUDA

void setup_quda_solverF2(QOP_FermionLinksAsqtad *fla);
void free_quda_solverF2(void);
void solve_qudaF2(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		     QOP_invert_arg_t *invarg, QOP_resid_arg_t *resarg,
		     QOP_evenodd_t eo);
void solveMulti_qudaF2(QDP_ColorVector *dest[], QDP_ColorVector *src,
			  double ms[], int nm,  QOP_invert_arg_t *invarg,
			  QOP_resid_arg_t *resargs[]);

#endif // HAVE_QUDA

#endif /* _QOP_F2_INTERNAL_H */
