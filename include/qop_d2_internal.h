// DO NOT EDIT
// generated from qop_pc_internal.h
#ifndef _QOP_D2_INTERNAL_H
#define _QOP_D2_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_D2_ColorVector_struct {
  QDP_D2_ColorVector *cv;
  QOP_D_Real *raw;
};

struct QOP_D2_DiracFermion_struct {
  QDP_D2_DiracFermion *df;
  QOP_D_Real *raw;
};

typedef void (*QOP_D2_gauge_deriv)(QDP_D2_ColorMatrix **d[],
				    QOP_D2_GaugeField *g,
				    QDP_D2_ColorMatrix *c[]);

typedef void (*QOP_D2_gauge_scale)(QDP_D2_ColorMatrix *l[],
				    QOP_D2_GaugeField *g, int inv);

struct QOP_D2_GaugeField_struct {
  QDP_D2_ColorMatrix **links;
  QOP_D_Real **raw;
  QOP_D2_GaugeField **parents;
  QOP_D2_gauge_deriv deriv;
  QOP_D2_gauge_scale scale;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t sign;
  int chained;
  int nparents;
};

struct QOP_D2_Force_struct {
  QDP_D2_ColorMatrix **force;
  QOP_D_Real **raw;
};

typedef struct {
  QDP_D2_ColorVector **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_D2_eigcg_t_V;

typedef struct {
  QDP_D2_DiracFermion **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_D2_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_D2_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_D2_ColorMatrix **fatlinks;
  QDP_D2_ColorMatrix **longlinks;
  QDP_D2_ColorMatrix **fwdlinks;
  QDP_D2_ColorMatrix **bcklinks;
  QDP_D2_ColorMatrix **dbllinks;
  QOP_D2_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
  //QOP_F2_FermionLinksAsqtad *ofla;
};

  /* HISQ datatypes */

struct QOP_D2_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_D2_ColorMatrix **U_links; // gauge links
  QDP_D2_ColorMatrix **V_links; // Fat7 smeared
  QDP_D2_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_D2_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_D2_FermionLinksAsqtad **fn;
  QOP_D2_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_D2_FermionLinksWilson_struct {
  QOP_D_Real clovinvkappa;
  int dblstored;
  QDP_D2_ColorMatrix **links;
  QDP_D2_ColorMatrix **bcklinks;
  QDP_D2_ColorMatrix **dbllinks;
  QOP_D2_GaugeField *qopgf;
  QOP_D2_GaugeField *gauge;
  QDP_D2_DiracPropagator *qdpclov;
  QOP_D_Real *clov, *clovinv;
  QOP_D_Real **rawlinks, *rawclov;
  QOP_D2_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_D2_FermionLinksDW_struct {
  QOP_D2_FermionLinksWilson *flw;
};

/* internal routines */

QOP_D2_FermionLinksAsqtad *QOP_D2_asqtad_create_L_from_L(QOP_D2_FermionLinksAsqtad *fla_src);
QOP_D2_FermionLinksAsqtad *QOP_D2_asqtad_create_L_from_r_times_L(QOP_D_Real s,
								  QOP_D2_FermionLinksAsqtad *fla_src);
void QOP_D2_asqtad_L_peq_L(QOP_D2_FermionLinksAsqtad *fla, QOP_D2_FermionLinksAsqtad *fla1);
void QOP_D2_qdpM_eq_raw(QDP_D2_ColorMatrix *cm, QOP_D_Real *lnk);

typedef void (QOP_D2_linop_t_V)(QDP_D2_ColorVector *out, QDP_D2_ColorVector *in, QDP_Subset subset);
typedef void (QOP_D2_linop_t_D)(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_D2_linop_t_vD)(QDP_D2_DiracFermion **out, QDP_D2_DiracFermion **in, QDP_Subset subset);

typedef QOP_D_Real (QOP_D2_linopn_t_V)(QDP_D2_ColorVector *out, QDP_D2_ColorVector *in, QDP_Subset subset);
typedef QOP_D_Real (QOP_D2_linopn_t_D)(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in, QDP_Subset subset);
typedef QOP_D_Real (QOP_D2_linopn_t_vD)(QDP_D2_DiracFermion **out, QDP_D2_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_D2_invert_cg_V(QOP_D2_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_D2_ColorVector *out,
		    QDP_D2_ColorVector *in,
		    QDP_D2_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_D2_invert_cg_D(QOP_D2_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_D2_DiracFermion *out,
		    QDP_D2_DiracFermion *in,
		    QDP_D2_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_D2_invert_cg_vD(QOP_D2_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_D2_DiracFermion **out,
		     QDP_D2_DiracFermion **in,
		     QDP_D2_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_D2_invert_cgms_V(QOP_D2_linopn_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_D2_ColorVector **out,
		      QDP_D2_ColorVector *in,
		      QDP_D2_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_D2_invert_cgms_D(QOP_D2_linopn_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_D2_DiracFermion **out,
		      QDP_D2_DiracFermion *in,
		      QDP_D2_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_D2_invert_cgms_vD(QOP_D2_linopn_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_D_Real *shifts,
		       int nshifts,
		       QDP_D2_DiracFermion ***out,
		       QDP_D2_DiracFermion **in,
		       QDP_D2_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_D2_invert_bicgstab_D(QOP_D2_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_D2_DiracFermion *out,
			  QDP_D2_DiracFermion *in,
			  QDP_D2_DiracFermion *p,
			  QDP_D2_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_D2_invert_eigcg_V(QOP_D2_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_D2_ColorVector *out,
		       QDP_D2_ColorVector *in,
		       QDP_D2_ColorVector *p,
		       QDP_Subset subset,
		       QOP_D2_eigcg_t_V *eigcg);

QOP_status_t
QOP_D2_invert_eigcg_D(QOP_D2_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_D2_DiracFermion *out,
		       QDP_D2_DiracFermion *in,
		       QDP_D2_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_D2_eigcg_t_D *eigcg);


QDP_D2_ColorVector *QOP_D2_asqtad_dslash_get_tmp(QOP_D2_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_D2_DiracFermion *QOP_D2_wilson_dslash_get_tmp(QOP_D2_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);
QOP_D2_FermionLinksWilson *QOP_D2_wilson_initialize_gauge_L(void);

void QOP_D2_get_mid(QOP_info_t *info, QDP_D2_ColorMatrix *mid[],
		     QDP_Shift shifts[], int ns, QOP_D_Real eps[],
		     QOP_D_Real scale, QDP_D2_ColorVector *x[], int nterms);

void QOP_D2_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QDP_D2_ColorMatrix *links[],
					  QDP_D2_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[], QDP_D2_ColorVector *x[], int nsrc);

void QOP_D2_asqtad_deriv_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_D2_ColorMatrix *links[],
					  QDP_D2_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[],
					  QDP_D2_ColorVector *x[],
					  int nterms);

void QOP_D2_asqtad_force_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_D2_ColorMatrix *links[],
					  QDP_D2_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_D_Real eps[],
					  QDP_D2_ColorVector *x[],
					  int nterms);

//AB internal operations for HISQ

void 
QOP_D2_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_D2_FermionLinksHisq *flh,
				       QOP_D2_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_D_Real *epsv,
				       QDP_D2_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_D2_hisq_deriv_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_D2_FermionLinksHisq *flh,
				    QDP_D2_ColorMatrix *deriv[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_D_Real *epsv,
				    QDP_D2_ColorVector *in_pt[], 
				    int *n_orders_naik);

void 
QOP_D2_hisq_force_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_D2_FermionLinksHisq *flh,
				    QDP_D2_ColorMatrix *force[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_D_Real *epsv,
				    QDP_D2_ColorVector *in_pt[], 
				    int *n_orders_naik);

void QOP_D2_dw_schur2_qdp(QOP_info_t *info, QOP_D2_FermionLinksDW *fldw,
			   QOP_D_Real M5, QOP_D_Real mq,
			   QDP_D2_DiracFermion *out[], QDP_D2_DiracFermion *in[],
			   int ls,QOP_evenodd_t eo);
void QOP_D2_dw_schur_qdp(QOP_info_t *info, QOP_D2_FermionLinksDW *fldw,
			  QOP_D_Real M5, QOP_D_Real mq, int sign,
			  QDP_D2_DiracFermion *out[], QDP_D2_DiracFermion *in[],
			  int ls, QOP_evenodd_t eo);
void
QOP_D2_dw_EO_project(QOP_D2_FermionLinksDW *fldw,
		      QDP_D2_DiracFermion *out[], QDP_D2_DiracFermion *in[],
		      QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);
void
QOP_D2_dw_EO_reconstruct(QOP_D2_FermionLinksDW *fldw,
			  QDP_D2_DiracFermion *out[], QDP_D2_DiracFermion *in[],
			  QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_D2_invert_gcr2_D(QOP_D2_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_D2_DiracFermion *out,
		      QDP_D2_DiracFermion *in,
		      QDP_D2_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_D2_invert_gmres2_D(QOP_D2_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_D2_DiracFermion *out,
			QDP_D2_DiracFermion *in,
			QDP_D2_DiracFermion *r,
			QDP_Subset subset);

QOP_D_Real
QOP_D2_relnorm2_V(QDP_D2_ColorVector **rsd, 
		   QDP_D2_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_D_Real
QOP_D2_relnorm2_D(QDP_D2_DiracFermion **rsd, 
		   QDP_D2_DiracFermion **out, 
		   QDP_Subset subset, int nv);

//// MULTIGRID STUFF

typedef struct {
  QOP_D2_FermionLinksWilson *wil;
  QLA_D_Real kappa;
} QOP_D2_WilArgs;

void QOP_D2_wilsonDslash(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in,
			  QOP_D2_FermionLinksWilson *wil, QLA_D_Real kappa,
			  int sign, QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_D2_wilsonDiaginv(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in,
			   QOP_D2_FermionLinksWilson *wil, QLA_D_Real kappa,
			   QOP_evenodd_t pout);
void QOP_D2_wilsonDslashEO(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in,
			    QOP_D2_FermionLinksWilson *wil, QLA_D_Real kappa,
			    int sign, QOP_evenodd_t par);
void QOP_D2_wilsonDslashEOS(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_D_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_D2_wilsonDslashEOH(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_D_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_D2_wilEoProjectD(QDP_D2_DiracFermion *ineo, QDP_D2_DiracFermion *in,
			   QOP_D2_WilArgs *w);
void QOP_D2_wilEoReconstructD(QDP_D2_DiracFermion *out, QDP_D2_DiracFermion *outeo,
			       QDP_D2_DiracFermion *in, QOP_D2_WilArgs *w);

#ifdef HAVE_NCN
#include <qdp_fn.h>
#include <qdp_dn.h>

void QOP_D2_V1eqD(QDP_DN_ColorVector *v[1], QDP_D2_DiracFermion *d, QDP_Subset sub);
void QOP_D2_DeqV1(QDP_D2_DiracFermion *d, QDP_DN_ColorVector *v[1], QDP_Subset sub);
void QOP_D2_V2eqD(QDP_DN_ColorVector *v[2], QDP_D2_DiracFermion *d, QDP_Subset sub);
void QOP_D2_DeqV2(QDP_D2_DiracFermion *d, QDP_DN_ColorVector *v[2], QDP_Subset sub);
void QOP_D2_wilDV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D2_wilDV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D2_wilPV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D2_wilPV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D2_wilPNEV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D2_wilEoV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *in[1], int sign, void *args);
void QOP_D2_wilEoV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *in[2], int sign, void *args);
void QOP_D2_wilEoProjectV1(QDP_DN_ColorVector *ineo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D2_wilEoReconstructV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D2_wilEoReconstructPV1(QDP_DN_ColorVector *out[1], QDP_DN_ColorVector *outeo[1], QDP_DN_ColorVector *in[1], void *args);
void QOP_D2_wilEoProjectV2(QDP_DN_ColorVector *ineo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_D2_wilEoReconstructV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);
void QOP_D2_wilEoReconstructPV2(QDP_DN_ColorVector *out[2], QDP_DN_ColorVector *outeo[2], QDP_DN_ColorVector *in[2], void *args);

#ifndef _QOP_2_MG_INTERNAL
#define _QOP_2_MG_INTERNAL

#include <qop_f_internal.h>
#include <qop_d_internal.h>
//#include <qop_mg_internal.h>
//struct QOP_WilMgLevel;

struct QOP_2_WilsonMgStruct {
  QOP_F2_FermionLinksWilson *wilF;
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

#if QOP_Precision == 'D'
#  if QOP_Colors == 2
#    include <qop_d2_internal_generic.h>
#  endif
#endif

#ifdef HAVE_QLL

void setup_qll_solverD2(QOP_FermionLinksAsqtad *fla);
void free_qll_solverD2(void);
void solve_qllD2(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		    QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg);
void solveMulti_qllD2(QDP_ColorVector *dest[], QDP_ColorVector *src,
			 double ms[], int nm,  QOP_invert_arg_t *invarg,
			 QOP_resid_arg_t *resargs[]);
void * create_qll_gaugeD2(int nc);
void * create_qll_from_gaugeD2(QDP_ColorMatrix *g[]);
void copy_gauge_from_qllD2(QDP_ColorMatrix *g[], void *ff);
void free_qll_gaugeD2(void *ff);
void fat7_qllD2(void *qllfl, void *qllll, QOP_asqtad_coeffs_t *coef,
		   void *qllu, void *qllul);

#endif // HAVE_QLL

#ifdef HAVE_QUDA

void setup_quda_solverD2(QOP_FermionLinksAsqtad *fla);
void free_quda_solverD2(void);
void solve_qudaD2(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		     QOP_invert_arg_t *invarg, QOP_resid_arg_t *resarg,
		     QOP_evenodd_t eo);
void solveMulti_qudaD2(QDP_ColorVector *dest[], QDP_ColorVector *src,
			  double ms[], int nm,  QOP_invert_arg_t *invarg,
			  QOP_resid_arg_t *resargs[]);

#endif // HAVE_QUDA

#endif /* _QOP_D2_INTERNAL_H */
