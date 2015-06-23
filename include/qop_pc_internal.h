#ifndef _QOP_PC_INTERNAL_H
#define _QOP_PC_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_IPC_ColorVector_struct {
  QDP_PC_ColorVector *cv;
  QOP_P_Real *raw;
};

struct QOP_IPC_DiracFermion_struct {
  QDP_PC_DiracFermion *df;
  QOP_P_Real *raw;
};

typedef void (*QOP_IPC_gauge_deriv)(QDP_PC_ColorMatrix **d[],
				    QOP_PC_GaugeField *g,
				    QDP_PC_ColorMatrix *c[]);

typedef void (*QOP_IPC_gauge_scale)(QDP_PC_ColorMatrix *l[],
				    QOP_PC_GaugeField *g, int inv);

struct QOP_IPC_GaugeField_struct {
  QDP_PC_ColorMatrix **links;
  QOP_P_Real **raw;
  QOP_PC_GaugeField **parents;
  QOP_IPC_gauge_deriv deriv;
  QOP_IPC_gauge_scale scale;
  int *r0;
  QOP_bc_t bc;
  QOP_staggered_sign_t sign;
  int chained;
  int nparents;
};

struct QOP_IPC_Force_struct {
  QDP_PC_ColorMatrix **force;
  QOP_P_Real **raw;
};

typedef struct {
  QDP_PC_ColorVector **u;
  QOP_P_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_IPC_eigcg_t_V;

typedef struct {
  QDP_PC_DiracFermion **u;
  QOP_P_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_IPC_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_IPC_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_PC_ColorMatrix **fatlinks;
  QDP_PC_ColorMatrix **longlinks;
  QDP_PC_ColorMatrix **fwdlinks;
  QDP_PC_ColorMatrix **bcklinks;
  QDP_PC_ColorMatrix **dbllinks;
  QOP_IPC_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
  //QOP_OC_FermionLinksAsqtad *ofla;
};

  /* HISQ datatypes */

struct QOP_IPC_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_PC_ColorMatrix **U_links; // gauge links
  QDP_PC_ColorMatrix **V_links; // Fat7 smeared
  QDP_PC_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_PC_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_PC_FermionLinksAsqtad **fn;
  QOP_PC_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_IPC_FermionLinksWilson_struct {
  QOP_P_Real clovinvkappa;
  int dblstored;
  QDP_PC_ColorMatrix **links;
  QDP_PC_ColorMatrix **bcklinks;
  QDP_PC_ColorMatrix **dbllinks;
  QOP_PC_GaugeField *qopgf;
  QOP_PC_GaugeField *gauge;
  QDP_PC_DiracPropagator *qdpclov;
  QOP_P_Real *clov, *clovinv;
  QOP_P_Real **rawlinks, *rawclov;
  QOP_IPC_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_IPC_FermionLinksDW_struct {
  QOP_PC_FermionLinksWilson *flw;
};

/* internal routines */

QOP_PC_FermionLinksAsqtad *QOP_IPC_asqtad_create_L_from_L(QOP_PC_FermionLinksAsqtad *fla_src);
QOP_PC_FermionLinksAsqtad *QOP_IPC_asqtad_create_L_from_r_times_L(QOP_D_Real s,
								  QOP_PC_FermionLinksAsqtad *fla_src);
void QOP_IPC_asqtad_L_peq_L(QOP_PC_FermionLinksAsqtad *fla, QOP_PC_FermionLinksAsqtad *fla1);
void QOP_IPC_qdpM_eq_raw(QDP_PC_ColorMatrix *cm, QOP_P_Real *lnk);

typedef void (QOP_IPC_linop_t_V)(QDP_PC_ColorVector *out, QDP_PC_ColorVector *in, QDP_Subset subset);
typedef void (QOP_IPC_linop_t_D)(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_IPC_linop_t_vD)(QDP_PC_DiracFermion **out, QDP_PC_DiracFermion **in, QDP_Subset subset);

typedef QOP_P_Real (QOP_IPC_linopn_t_V)(QDP_PC_ColorVector *out, QDP_PC_ColorVector *in, QDP_Subset subset);
typedef QOP_P_Real (QOP_IPC_linopn_t_D)(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in, QDP_Subset subset);
typedef QOP_P_Real (QOP_IPC_linopn_t_vD)(QDP_PC_DiracFermion **out, QDP_PC_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_cg_V(QOP_IPC_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_PC_ColorVector *out,
		    QDP_PC_ColorVector *in,
		    QDP_PC_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_cg_D(QOP_IPC_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_PC_DiracFermion *out,
		    QDP_PC_DiracFermion *in,
		    QDP_PC_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_cg_vD(QOP_IPC_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_PC_DiracFermion **out,
		     QDP_PC_DiracFermion **in,
		     QDP_PC_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_IPC_invert_cgms_V(QOP_IPC_linopn_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_P_Real *shifts,
		      int nshifts,
		      QDP_PC_ColorVector **out,
		      QDP_PC_ColorVector *in,
		      QDP_PC_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_cgms_D(QOP_IPC_linopn_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_P_Real *shifts,
		      int nshifts,
		      QDP_PC_DiracFermion **out,
		      QDP_PC_DiracFermion *in,
		      QDP_PC_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_cgms_vD(QOP_IPC_linopn_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_P_Real *shifts,
		       int nshifts,
		       QDP_PC_DiracFermion ***out,
		       QDP_PC_DiracFermion **in,
		       QDP_PC_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_IPC_invert_bicgstab_D(QOP_IPC_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_PC_DiracFermion *out,
			  QDP_PC_DiracFermion *in,
			  QDP_PC_DiracFermion *p,
			  QDP_PC_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_eigcg_V(QOP_IPC_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_PC_ColorVector *out,
		       QDP_PC_ColorVector *in,
		       QDP_PC_ColorVector *p,
		       QDP_Subset subset,
		       QOP_IPC_eigcg_t_V *eigcg);

QOP_status_t
QOP_IPC_invert_eigcg_D(QOP_IPC_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_PC_DiracFermion *out,
		       QDP_PC_DiracFermion *in,
		       QDP_PC_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_IPC_eigcg_t_D *eigcg);


QDP_PC_ColorVector *QOP_IPC_asqtad_dslash_get_tmp(QOP_PC_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_PC_DiracFermion *QOP_IPC_wilson_dslash_get_tmp(QOP_PC_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);
QOP_PC_FermionLinksWilson *QOP_IPC_wilson_initialize_gauge_L(void);

void QOP_IPC_get_mid(QOP_info_t *info, QDP_PC_ColorMatrix *mid[],
		     QDP_Shift shifts[], int ns, QOP_P_Real eps[],
		     QOP_P_Real scale, QDP_PC_ColorVector *x[], int nterms);

void QOP_IPC_asqtad_force_multi_asvec_qdp(QOP_info_t *info, QDP_PC_ColorMatrix *links[],
					  QDP_PC_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
					  QOP_P_Real eps[], QDP_PC_ColorVector *x[], int nsrc);

void QOP_IPC_asqtad_deriv_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_PC_ColorMatrix *links[],
					  QDP_PC_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_P_Real eps[],
					  QDP_PC_ColorVector *x[],
					  int nterms);

void QOP_IPC_asqtad_force_multi_fnmat_qdp(QOP_info_t *info,
					  QDP_PC_ColorMatrix *links[],
					  QDP_PC_ColorMatrix *force[],
					  QOP_asqtad_coeffs_t *coef,
					  QOP_P_Real eps[],
					  QDP_PC_ColorVector *x[],
					  int nterms);

//AB internal operations for HISQ

void 
QOP_IPC_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_PC_FermionLinksHisq *flh,
				       QOP_PC_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_P_Real *epsv,
				       QDP_PC_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_IPC_hisq_deriv_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_PC_FermionLinksHisq *flh,
				    QDP_PC_ColorMatrix *deriv[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_P_Real *epsv,
				    QDP_PC_ColorVector *in_pt[], 
				    int *n_orders_naik);

void 
QOP_IPC_hisq_force_multi_fnmat2_qdp(QOP_info_t *info,  
				    QOP_PC_FermionLinksHisq *flh,
				    QDP_PC_ColorMatrix *force[],
				    QOP_hisq_coeffs_t *hisq_coeff,
				    QOP_P_Real *epsv,
				    QDP_PC_ColorVector *in_pt[], 
				    int *n_orders_naik);

void QOP_IPC_dw_schur2_qdp(QOP_info_t *info, QOP_PC_FermionLinksDW *fldw,
			   QOP_P_Real M5, QOP_P_Real mq,
			   QDP_PC_DiracFermion *out[], QDP_PC_DiracFermion *in[],
			   int ls,QOP_evenodd_t eo);
void QOP_IPC_dw_schur_qdp(QOP_info_t *info, QOP_PC_FermionLinksDW *fldw,
			  QOP_P_Real M5, QOP_P_Real mq, int sign,
			  QDP_PC_DiracFermion *out[], QDP_PC_DiracFermion *in[],
			  int ls, QOP_evenodd_t eo);
void
QOP_IPC_dw_EO_project(QOP_PC_FermionLinksDW *fldw,
		      QDP_PC_DiracFermion *out[], QDP_PC_DiracFermion *in[],
		      QOP_P_Real M5, QOP_P_Real mq, int ls, QOP_evenodd_t eo);
void
QOP_IPC_dw_EO_reconstruct(QOP_PC_FermionLinksDW *fldw,
			  QDP_PC_DiracFermion *out[], QDP_PC_DiracFermion *in[],
			  QOP_P_Real M5, QOP_P_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_IPC_invert_gcr2_D(QOP_IPC_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_PC_DiracFermion *out,
		      QDP_PC_DiracFermion *in,
		      QDP_PC_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_IPC_invert_gmres2_D(QOP_IPC_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_PC_DiracFermion *out,
			QDP_PC_DiracFermion *in,
			QDP_PC_DiracFermion *r,
			QDP_Subset subset);

QOP_P_Real
QOP_IPC_relnorm2_V(QDP_PC_ColorVector **rsd, 
		   QDP_PC_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_P_Real
QOP_IPC_relnorm2_D(QDP_PC_DiracFermion **rsd, 
		   QDP_PC_DiracFermion **out, 
		   QDP_Subset subset, int nv);

//// MULTIGRID STUFF

typedef struct {
  QOP_PC_FermionLinksWilson *wil;
  QLA_P_Real kappa;
} QOP_IPC_WilArgs;

void QOP_IPC_wilsonDslash(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in,
			  QOP_PC_FermionLinksWilson *wil, QLA_P_Real kappa,
			  int sign, QOP_evenodd_t pout, QOP_evenodd_t pin);
void QOP_IPC_wilsonDiaginv(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in,
			   QOP_PC_FermionLinksWilson *wil, QLA_P_Real kappa,
			   QOP_evenodd_t pout);
void QOP_IPC_wilsonDslashEO(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in,
			    QOP_PC_FermionLinksWilson *wil, QLA_P_Real kappa,
			    int sign, QOP_evenodd_t par);
void QOP_IPC_wilsonDslashEOS(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_P_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_IPC_wilsonDslashEOH(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *in,
			     QOP_FermionLinksWilson *wil, QLA_P_Real kappa,
			     int sign, QOP_evenodd_t par);
void QOP_IPC_wilEoProjectD(QDP_PC_DiracFermion *ineo, QDP_PC_DiracFermion *in,
			   QOP_IPC_WilArgs *w);
void QOP_IPC_wilEoReconstructD(QDP_PC_DiracFermion *out, QDP_PC_DiracFermion *outeo,
			       QDP_PC_DiracFermion *in, QOP_IPC_WilArgs *w);

#ifdef HAVE_NCN
#include <qdp_fn.h>
#include <qdp_dn.h>

void QOP_IPC_V1eqD(QDP_N_ColorVector *v[1], QDP_PC_DiracFermion *d, QDP_Subset sub);
void QOP_IPC_DeqV1(QDP_PC_DiracFermion *d, QDP_N_ColorVector *v[1], QDP_Subset sub);
void QOP_IPC_V2eqD(QDP_N_ColorVector *v[2], QDP_PC_DiracFermion *d, QDP_Subset sub);
void QOP_IPC_DeqV2(QDP_PC_DiracFermion *d, QDP_N_ColorVector *v[2], QDP_Subset sub);
void QOP_IPC_wilDV1(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1], int sign, void *args);
void QOP_IPC_wilDV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2], int sign, void *args);
void QOP_IPC_wilPV1(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1], int sign, void *args);
void QOP_IPC_wilPV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2], int sign, void *args);
void QOP_IPC_wilPNEV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2], int sign, void *args);
void QOP_IPC_wilEoV1(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1], int sign, void *args);
void QOP_IPC_wilEoV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2], int sign, void *args);
void QOP_IPC_wilEoProjectV1(QDP_N_ColorVector *ineo[1], QDP_N_ColorVector *in[1], void *args);
void QOP_IPC_wilEoReconstructV1(QDP_N_ColorVector *out[1], QDP_N_ColorVector *outeo[1], QDP_N_ColorVector *in[1], void *args);
void QOP_IPC_wilEoReconstructPV1(QDP_N_ColorVector *out[1], QDP_N_ColorVector *outeo[1], QDP_N_ColorVector *in[1], void *args);
void QOP_IPC_wilEoProjectV2(QDP_N_ColorVector *ineo[2], QDP_N_ColorVector *in[2], void *args);
void QOP_IPC_wilEoReconstructV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *outeo[2], QDP_N_ColorVector *in[2], void *args);
void QOP_IPC_wilEoReconstructPV2(QDP_N_ColorVector *out[2], QDP_N_ColorVector *outeo[2], QDP_N_ColorVector *in[2], void *args);

#endif // HAVE_NCN

#ifndef _QOP_C_MG_INTERNAL
#define _QOP_C_MG_INTERNAL

#include <qop_f_internal.h>
#include <qop_d_internal.h>
//#include <qop_mg_internal.h>
//struct QOP_WilMgLevel;

struct QOP_C_WilsonMgStruct {
  QDP_Lattice *qdp_lattice;
  QOP_FC_FermionLinksWilson *wilF;
  QOP_DC_FermionLinksWilson *wilD;
  QLA_F_Real kappa;
  QLA_F_Real kappanv;
  QOP_FC_WilArgs vcwaF;
  QOP_FC_WilArgs nvwaF;
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

#endif // _QOP_C_MG_INTERNAL

#if QOP_Precision == _QOP_Precision
#  if QOP_Colors == _QOP_Colors
#    include <qop_pc_internal_generic.h>
#  endif
#endif

#ifdef HAVE_QLL

void IPC(setup_qll_solver)(QOP_FermionLinksAsqtad *fla);
void IPC(free_qll_solver)(void);
void IPC(solve_qll)(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
		    QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg);
void IPC(solveMulti_qll)(QDP_ColorVector *dest[], QDP_ColorVector *src,
			 double ms[], int nm,  QOP_invert_arg_t *invarg,
			 QOP_resid_arg_t *resargs[]);
void * IPC(create_qll_gauge)(int nc);
void * IPC(create_qll_from_gauge)(QDP_ColorMatrix *g[]);
void IPC(copy_gauge_from_qll)(QDP_ColorMatrix *g[], void *ff);
void IPC(free_qll_gauge)(void *ff);
void IPC(fat7_qll)(void *qllfl, void *qllll, QOP_asqtad_coeffs_t *coef,
		   void *qllu, void *qllul);

#endif // HAVE_QLL

#endif /* _QOP_PC_INTERNAL_H */
