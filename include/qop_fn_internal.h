// DO NOT EDIT
// generated from qop_pc_internal.h
#ifndef _QOP_FN_INTERNAL_H
#define _QOP_FN_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_FN_ColorVector_struct {
  QDP_FN_ColorVector *cv;
  QOP_F_Real *raw;
};

struct QOP_FN_DiracFermion_struct {
  QDP_FN_DiracFermion *df;
  QOP_F_Real *raw;
};

typedef void (*QOP_FN_gauge_deriv)(QDP_FN_ColorMatrix **d[],
				    QOP_FN_GaugeField *g,
				    QDP_FN_ColorMatrix *c[]);

typedef void (*QOP_FN_gauge_scale)(QDP_FN_ColorMatrix *l[],
				    QOP_FN_GaugeField *g, int inv);

struct QOP_FN_GaugeField_struct {
  QDP_FN_ColorMatrix **links;
  QOP_F_Real **raw;
  QOP_FN_GaugeField **parents;
  QOP_FN_gauge_deriv deriv;
  QOP_FN_gauge_scale scale;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t sign;
  int chained;
  int nparents;
};

struct QOP_FN_Force_struct {
  QDP_FN_ColorMatrix **force;
  QOP_F_Real **raw;
};

typedef struct {
  QDP_FN_ColorVector **u;
  QOP_F_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_FN_eigcg_t_V;

typedef struct {
  QDP_FN_DiracFermion **u;
  QOP_F_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_FN_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_FN_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_FN_ColorMatrix **fatlinks;
  QDP_FN_ColorMatrix **longlinks;
  QDP_FN_ColorMatrix **fwdlinks;
  QDP_FN_ColorMatrix **bcklinks;
  QDP_FN_ColorMatrix **dbllinks;
  QOP_FN_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
  //QOP_DN_FermionLinksAsqtad *ofla;
};

  /* HISQ datatypes */

struct QOP_FN_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_FN_ColorMatrix **U_links; // gauge links
  QDP_FN_ColorMatrix **V_links; // Fat7 smeared
  QDP_FN_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_FN_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_FN_FermionLinksAsqtad **fn;
  QOP_FN_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_FN_FermionLinksWilson_struct {
  QOP_F_Real clovinvkappa;
  int dblstored;
  QDP_FN_ColorMatrix **links;
  QDP_FN_ColorMatrix **bcklinks;
  QDP_FN_ColorMatrix **dbllinks;
  QOP_FN_GaugeField *qopgf;
  QOP_FN_GaugeField *gauge;
  QDP_FN_DiracPropagator *qdpclov;
  QOP_F_Real *clov, *clovinv;
  QOP_F_Real **rawlinks, *rawclov;
  QOP_FN_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_FN_FermionLinksDW_struct {
  QOP_FN_FermionLinksWilson *flw;
};

/* internal routines */

QOP_FN_FermionLinksAsqtad *QOP_FN_asqtad_create_L_from_L(QOP_FN_FermionLinksAsqtad *fla_src);
QOP_FN_FermionLinksAsqtad *QOP_FN_asqtad_create_L_from_r_times_L(QOP_D_Real s,
								  QOP_FN_FermionLinksAsqtad *fla_src);
void QOP_FN_asqtad_L_peq_L(QOP_FN_FermionLinksAsqtad *fla, QOP_FN_FermionLinksAsqtad *fla1);
void QOP_FN_qdpM_eq_raw(QDP_FN_ColorMatrix *cm, QOP_F_Real *lnk);

typedef void (QOP_FN_linop_t_V)(QDP_FN_ColorVector *out, QDP_FN_ColorVector *in, QDP_Subset subset);
typedef void (QOP_FN_linop_t_D)(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_FN_linop_t_vD)(QDP_FN_DiracFermion **out, QDP_FN_DiracFermion **in, QDP_Subset subset);

typedef QOP_F_Real (QOP_FN_linopn_t_V)(QDP_FN_ColorVector *out, QDP_FN_ColorVector *in, QDP_Subset subset);
typedef QOP_F_Real (QOP_FN_linopn_t_D)(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in, QDP_Subset subset);
typedef QOP_F_Real (QOP_FN_linopn_t_vD)(QDP_FN_DiracFermion **out, QDP_FN_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_FN_invert_cg_V(QOP_FN_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_FN_ColorVector *out,
		    QDP_FN_ColorVector *in,
		    QDP_FN_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_FN_invert_cg_D(QOP_FN_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_FN_DiracFermion *out,
		    QDP_FN_DiracFermion *in,
		    QDP_FN_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_FN_invert_cg_vD(QOP_FN_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_FN_DiracFermion **out,
		     QDP_FN_DiracFermion **in,
		     QDP_FN_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_FN_invert_cgms_V(QOP_FN_linopn_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_FN_ColorVector **out,
		      QDP_FN_ColorVector *in,
		      QDP_FN_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_FN_invert_cgms_D(QOP_FN_linopn_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_F_Real *shifts,
		      int nshifts,
		      QDP_FN_DiracFermion **out,
		      QDP_FN_DiracFermion *in,
		      QDP_FN_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_FN_invert_cgms_vD(QOP_FN_linopn_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_F_Real *shifts,
		       int nshifts,
		       QDP_FN_DiracFermion ***out,
		       QDP_FN_DiracFermion **in,
		       QDP_FN_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_FN_invert_bicgstab_D(QOP_FN_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_FN_DiracFermion *out,
			  QDP_FN_DiracFermion *in,
			  QDP_FN_DiracFermion *p,
			  QDP_FN_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_FN_invert_eigcg_V(QOP_FN_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_FN_ColorVector *out,
		       QDP_FN_ColorVector *in,
		       QDP_FN_ColorVector *p,
		       QDP_Subset subset,
		       QOP_FN_eigcg_t_V *eigcg);

QOP_status_t
QOP_FN_invert_eigcg_D(QOP_FN_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_FN_DiracFermion *out,
		       QDP_FN_DiracFermion *in,
		       QDP_FN_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_FN_eigcg_t_D *eigcg);


QDP_FN_ColorVector *QOP_FN_asqtad_dslash_get_tmp(QOP_FN_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_FN_DiracFermion *QOP_FN_wilson_dslash_get_tmp(QOP_FN_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);
QOP_FN_FermionLinksWilson *QOP_FN_wilson_initialize_gauge_L(void);

void QOP_FN_get_mid(QOP_info_t *info, QDP_FN_ColorMatrix *mid[],
		     QDP_Shift shifts[], int ns, QOP_F_Real eps[],
		     QOP_F_Real scale, QDP_FN_ColorVector *x[], int nterms);

void QOP_FN_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QDP_FN_ColorMatrix *links[],
					  QDP_FN_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[], QDP_FN_ColorVector *x[], int nsrc);

void QOP_FN_asqtad_deriv_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_FN_ColorMatrix *links[],
					  QDP_FN_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[],
					  QDP_FN_ColorVector *x[],
					  int nterms);

void QOP_FN_asqtad_force_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_FN_ColorMatrix *links[],
					  QDP_FN_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_F_Real eps[],
					  QDP_FN_ColorVector *x[],
					  int nterms);

//AB internal operations for HISQ

void 
QOP_FN_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_FN_FermionLinksHisq *flh,
				       QOP_FN_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_F_Real *epsv,
				       QDP_FN_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_FN_hisq_deriv_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_FN_FermionLinksHisq *flh,
				    QDP_FN_ColorMatrix *deriv[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_F_Real *epsv,
				    QDP_FN_ColorVector *in_pt[], 
				    int *n_orders_naik);

void 
QOP_FN_hisq_force_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_FN_FermionLinksHisq *flh,
				    QDP_FN_ColorMatrix *force[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_F_Real *epsv,
				    QDP_FN_ColorVector *in_pt[], 
				    int *n_orders_naik);

void QOP_FN_dw_schur2_qdp(QOP_info_t *info, QOP_FN_FermionLinksDW *fldw,
			   QOP_F_Real M5, QOP_F_Real mq,
			   QDP_FN_DiracFermion *out[], QDP_FN_DiracFermion *in[],
			   int ls,QOP_evenodd_t eo);
void QOP_FN_dw_schur_qdp(QOP_info_t *info, QOP_FN_FermionLinksDW *fldw,
			  QOP_F_Real M5, QOP_F_Real mq, int sign,
			  QDP_FN_DiracFermion *out[], QDP_FN_DiracFermion *in[],
			  int ls, QOP_evenodd_t eo);
void
QOP_FN_dw_EO_project(QOP_FN_FermionLinksDW *fldw,
		      QDP_FN_DiracFermion *out[], QDP_FN_DiracFermion *in[],
		      QOP_F_Real M5, QOP_F_Real mq, int ls, QOP_evenodd_t eo);
void
QOP_FN_dw_EO_reconstruct(QOP_FN_FermionLinksDW *fldw,
			  QDP_FN_DiracFermion *out[], QDP_FN_DiracFermion *in[],
			  QOP_F_Real M5, QOP_F_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_FN_invert_gcr2_D(QOP_FN_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_FN_DiracFermion *out,
		      QDP_FN_DiracFermion *in,
		      QDP_FN_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_FN_invert_gmres2_D(QOP_FN_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_FN_DiracFermion *out,
			QDP_FN_DiracFermion *in,
			QDP_FN_DiracFermion *r,
			QDP_Subset subset);

QOP_F_Real
QOP_FN_relnorm2_V(QDP_FN_ColorVector **rsd, 
		   QDP_FN_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_F_Real
QOP_FN_relnorm2_D(QDP_FN_DiracFermion **rsd, 
		   QDP_FN_DiracFermion **out, 
		   QDP_Subset subset, int nv);

//// MULTIGRID STUFF

typedef struct {
  QOP_FN_FermionLinksWilson *wil;
  QLA_F_Real kappa;
} QOP_FN_WilArgs;

void QOP_FN_wilsonDslash(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in,
			  QOP_FN_FermionLinksWilson *wil, QLA_F_Real kappa,
			  int sign, QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_FN_wilsonDiaginv(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in,
			   QOP_FN_FermionLinksWilson *wil, QLA_F_Real kappa,
			   QOP_evenodd_t pout);
void QOP_FN_wilsonDslashEO(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in,
			    QOP_FN_FermionLinksWilson *wil, QLA_F_Real kappa,
			    int sign, QOP_evenodd_t par);
void QOP_FN_wilsonDslashEOS(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_F_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_FN_wilsonDslashEOH(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_F_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_FN_wilEoProjectD(QDP_FN_DiracFermion *ineo, QDP_FN_DiracFermion *in,
			   QOP_FN_WilArgs *w);
void QOP_FN_wilEoReconstructD(QDP_FN_DiracFermion *out, QDP_FN_DiracFermion *outeo,
			       QDP_FN_DiracFermion *in, QOP_FN_WilArgs *w);

#ifdef HAVE_NCN
#include <qdp_fn.h>
#include <qdp_dn.h>

void QOP_FN_V1eqD(QDP_FN_ColorVector *v[1], QDP_FN_DiracFermion *d, QDP_Subset sub);
void QOP_FN_DeqV1(QDP_FN_DiracFermion *d, QDP_FN_ColorVector *v[1], QDP_Subset sub);
void QOP_FN_V2eqD(QDP_FN_ColorVector *v[2], QDP_FN_DiracFermion *d, QDP_Subset sub);
void QOP_FN_DeqV2(QDP_FN_DiracFermion *d, QDP_FN_ColorVector *v[2], QDP_Subset sub);
void QOP_FN_wilDV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_FN_wilDV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_FN_wilPV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_FN_wilPV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_FN_wilPNEV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_FN_wilEoV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *in[1], int sign, void *args);
void QOP_FN_wilEoV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *in[2], int sign, void *args);
void QOP_FN_wilEoProjectV1(QDP_FN_ColorVector *ineo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_FN_wilEoReconstructV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_FN_wilEoReconstructPV1(QDP_FN_ColorVector *out[1], QDP_FN_ColorVector *outeo[1], QDP_FN_ColorVector *in[1], void *args);
void QOP_FN_wilEoProjectV2(QDP_FN_ColorVector *ineo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_FN_wilEoReconstructV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);
void QOP_FN_wilEoReconstructPV2(QDP_FN_ColorVector *out[2], QDP_FN_ColorVector *outeo[2], QDP_FN_ColorVector *in[2], void *args);

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

#endif // HAVE_NCN

#if QOP_Precision == 'F'
#  if QOP_Colors == 'N'
#    include <qop_fn_internal_generic.h>
#  endif
#endif

#ifdef HAVE_QLL

void setup_qll_solverFN(QOP_FermionLinksAsqtad *fla);
void free_qll_solverFN(void);
void solve_qllFN(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		    QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg);
void solveMulti_qllFN(QDP_ColorVector *dest[], QDP_ColorVector *src,
			 double ms[], int nm,  QOP_invert_arg_t *invarg,
			 QOP_resid_arg_t *resargs[]);
void * create_qll_gaugeFN(int nc);
void * create_qll_from_gaugeFN(QDP_ColorMatrix *g[]);
void copy_gauge_from_qllFN(QDP_ColorMatrix *g[], void *ff);
void free_qll_gaugeFN(void *ff);
void fat7_qllFN(void *qllfl, void *qllll, QOP_asqtad_coeffs_t *coef,
		   void *qllu, void *qllul);

#endif // HAVE_QLL

#endif /* _QOP_FN_INTERNAL_H */
