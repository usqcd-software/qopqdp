// DO NOT EDIT
// generated from qop_pc_internal.h
#ifndef _QOP_D3_INTERNAL_H
#define _QOP_D3_INTERNAL_H

//typedef struct {
//  int tmp;
//} QOPPC(common_t);
//extern QOPPC(common_t) QOPPC(common);

struct QOP_D3_ColorVector_struct {
  QDP_D3_ColorVector *cv;
  QOP_D_Real *raw;
};

struct QOP_D3_DiracFermion_struct {
  QDP_D3_DiracFermion *df;
  QOP_D_Real *raw;
};

struct QOP_D3_GaugeField_struct {
  QDP_D3_ColorMatrix **links;
  QOP_D_Real **raw;
};

struct QOP_D3_Force_struct {
  QDP_D3_ColorMatrix **force;
  QOP_D_Real **raw;
};

typedef struct {
  QDP_D3_ColorVector **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_D3_eigcg_t_V;

typedef struct {
  QDP_D3_DiracFermion **u;
  QOP_D_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_D3_eigcg_t_D;

  /* Asqtad datatypes */

struct QOP_D3_FermionLinksAsqtad_struct {
  int dblstored, nlinks;
  QDP_D3_ColorMatrix **fatlinks;
  QDP_D3_ColorMatrix **longlinks;
  QDP_D3_ColorMatrix **fwdlinks;
  QDP_D3_ColorMatrix **bcklinks;
  QDP_D3_ColorMatrix **dbllinks;
  QOP_D3_eigcg_t_V eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
};

  /* HISQ datatypes */

struct QOP_D3_FermionLinksHisq_struct {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDP_D3_ColorMatrix **U_links; // gauge links
  QDP_D3_ColorMatrix **V_links; // Fat7 smeared
  QDP_D3_ColorMatrix **Y_unitlinks; // projected to U(3),
  QDP_D3_ColorMatrix **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOP_D3_FermionLinksAsqtad **fn;
  QOP_D3_FermionLinksAsqtad *fn_deps;
};

  /* Wilson datatypes */

struct QOP_D3_FermionLinksWilson_struct {
  QOP_D_Real clovinvkappa;
  int dblstored;
  QDP_D3_ColorMatrix **links;
  QDP_D3_ColorMatrix **bcklinks;
  QDP_D3_ColorMatrix **dbllinks;
  QOP_D3_GaugeField *qopgf;
  QDP_D3_DiracPropagator *qdpclov;
  QOP_D_Real *clov, *clovinv;
  QOP_D_Real **rawlinks, *rawclov;
  QOP_D3_eigcg_t_D eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOP_D3_FermionLinksDW_struct {
  QOP_D3_FermionLinksWilson *flw;
};

/* internal routines */

QOP_D3_FermionLinksAsqtad *QOP_D3_asqtad_create_L_from_L(QOP_D3_FermionLinksAsqtad *fla_src);
QOP_D3_FermionLinksAsqtad *QOP_D3_asqtad_create_L_from_r_times_L(QOP_D_Real *s,
								  QOP_D3_FermionLinksAsqtad *fla_src);
void QOP_D3_asqtad_L_peq_L(QOP_D3_FermionLinksAsqtad *fla, QOP_D3_FermionLinksAsqtad *fla1);
void QOP_D3_qdpM_eq_raw(QDP_D3_ColorMatrix *cm, QOP_D_Real *lnk);
typedef void (QOP_D3_linop_t_V)(QDP_D3_ColorVector *out, QDP_D3_ColorVector *in, QDP_Subset subset);
typedef void (QOP_D3_linop_t_D)(QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in, QDP_Subset subset);
typedef void (QOP_D3_linop_t_vD)(QDP_D3_DiracFermion **out, QDP_D3_DiracFermion **in, QDP_Subset subset);

QOP_status_t
QOP_D3_invert_cg_V(QOP_D3_linop_t_V *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_D3_ColorVector *out,
		    QDP_D3_ColorVector *in,
		    QDP_D3_ColorVector *p,
		    QDP_Subset subset);

QOP_status_t
QOP_D3_invert_cg_D(QOP_D3_linop_t_D *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_D3_DiracFermion *out,
		    QDP_D3_DiracFermion *in,
		    QDP_D3_DiracFermion *p,
		    QDP_Subset subset);

QOP_status_t
QOP_D3_invert_cg_vD(QOP_D3_linop_t_vD *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_D3_DiracFermion **out,
		     QDP_D3_DiracFermion **in,
		     QDP_D3_DiracFermion **p,
		     QDP_Subset subset,
		     int _n);

QOP_status_t
QOP_D3_invert_cgms_V(QOP_D3_linop_t_V *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_D3_ColorVector **out,
		      QDP_D3_ColorVector *in,
		      QDP_D3_ColorVector *p,
		      QDP_Subset subset);

QOP_status_t
QOP_D3_invert_cgms_D(QOP_D3_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QOP_D_Real *shifts,
		      int nshifts,
		      QDP_D3_DiracFermion **out,
		      QDP_D3_DiracFermion *in,
		      QDP_D3_DiracFermion *p,
		      QDP_Subset subset);

QOP_status_t
QOP_D3_invert_cgms_vD(QOP_D3_linop_t_vD *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t **res_arg,
		       QOP_D_Real *shifts,
		       int nshifts,
		       QDP_D3_DiracFermion ***out,
		       QDP_D3_DiracFermion **in,
		       QDP_D3_DiracFermion **p,
		       QDP_Subset subset,
		       int _n);

QOP_status_t
QOP_D3_invert_bicgstab_D(QOP_D3_linop_t_D *linop,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  QDP_D3_DiracFermion *out,
			  QDP_D3_DiracFermion *in,
			  QDP_D3_DiracFermion *p,
			  QDP_D3_DiracFermion *r,
			  QDP_Subset subset);

QOP_status_t
QOP_D3_invert_eigcg_V(QOP_D3_linop_t_V *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_D3_ColorVector *out,
		       QDP_D3_ColorVector *in,
		       QDP_D3_ColorVector *p,
		       QDP_Subset subset,
		       QOP_D3_eigcg_t_V *eigcg);

QOP_status_t
QOP_D3_invert_eigcg_D(QOP_D3_linop_t_D *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_D3_DiracFermion *out,
		       QDP_D3_DiracFermion *in,
		       QDP_D3_DiracFermion *p,
		       QDP_Subset subset,
		       QOP_D3_eigcg_t_D *eigcg);


QDP_D3_ColorVector *QOP_D3_asqtad_dslash_get_tmp(QOP_D3_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);
QDP_D3_DiracFermion *QOP_D3_wilson_dslash_get_tmp(QOP_D3_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);

void QOP_D3_get_mid(QOP_info_t *info, QDP_D3_ColorMatrix *mid[], QDP_Shift shifts[], int ns,
		     QOP_D_Real eps[], QDP_D3_ColorVector *x[], int nterms);

void QOP_D3_asqtad_deriv(QOP_info_t *info, QDP_D3_ColorMatrix *gauge[],
			  QDP_D3_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
			  QDP_D3_ColorMatrix *mid_fat[], QDP_D3_ColorMatrix *mid_naik[]);

void QOP_D3_asqtad_force_multi_fnmat(QOP_info_t *info, QOP_D3_GaugeField *gauge,
				      QOP_D3_Force *force, QOP_asqtad_coeffs_t *coef,
				      QOP_D_Real eps[], QDP_D3_ColorVector *x[], int nterms);

//AB internal operations for HISQ

void
QOP_D3_hisq_force_multi_reunit(QOP_info_t *info,
				QDP_D3_ColorMatrix *gf[4],
				QDP_D3_ColorMatrix *force_accum[4],
				QDP_D3_ColorMatrix *force_accum_old[4]);

void 
QOP_D3_hisq_force_multi_wrapper_fnmat(QOP_info_t *info,  
				       QOP_D3_FermionLinksHisq *flh,
				       QOP_D3_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       QOP_D_Real *epsv,
				       QDP_D3_ColorVector *in_pt[], 
				       int *n_orders_naik);

void 
QOP_D3_hisq_deriv_multi_wrapper_fnmat2(QOP_info_t *info,  
					QOP_D3_FermionLinksHisq *flh,
					QOP_D3_Force *Force, 
					QOP_hisq_coeffs_t *hisq_coeff,
					QOP_D_Real *epsv,
					QDP_D3_ColorVector *in_pt[], 
					int *n_orders_naik);

void 
QOP_D3_hisq_force_multi_wrapper_fnmat2(QOP_info_t *info,  
					QOP_D3_FermionLinksHisq *flh,
					QOP_D3_Force *Force, 
					QOP_hisq_coeffs_t *hisq_coeff,
					QOP_D_Real *epsv,
					QDP_D3_ColorVector *in_pt[], 
					int *n_orders_naik);

void QOP_D3_u3reunit(QOP_info_t *info, QDP_D3_ColorMatrix *U, QDP_D3_ColorMatrix *V);

void QOP_D3_su3reunit(QOP_info_t *info, QDP_D3_ColorMatrix *U, QDP_D3_ColorMatrix *Ur);

void
QOP_D3_dw_schur2_qdp(QOP_info_t *info, QOP_D3_FermionLinksDW *fldw,
		      QOP_D_Real M5, QOP_D_Real mq,
		      QDP_D3_DiracFermion *out[], QDP_D3_DiracFermion *in[],
		      int ls, QOP_evenodd_t eo);
void
QOP_D3_dw_schur_qdp(QOP_info_t *info, QOP_D3_FermionLinksDW *fldw,
		     QOP_D_Real M5, QOP_D_Real mq, int sign,
		     QDP_DiracFermion *out[], QDP_D3_DiracFermion *in[],
		     int ls, QOP_evenodd_t eo);
extern void
QOP_D3_dw_EO_project(QOP_D3_FermionLinksDW *fldw,
		      QDP_D3_DiracFermion *out[], QDP_D3_DiracFermion *in[],
		      QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);
extern void
QOP_D3_dw_EO_reconstruct(QOP_D3_FermionLinksDW *fldw,
			  QDP_D3_DiracFermion *out[], QDP_D3_DiracFermion *in[],
			  QOP_D_Real M5, QOP_D_Real mq, int ls, QOP_evenodd_t eo);

QOP_status_t
QOP_D3_invert_gcr2_D(QOP_D3_linop_t_D *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_D3_DiracFermion *out,
		      QDP_D3_DiracFermion *in,
		      QDP_D3_DiracFermion *r,
		      QDP_Subset subset);

QOP_status_t
QOP_D3_invert_gmres2_D(QOP_D3_linop_t_D *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			QDP_D3_DiracFermion *out,
			QDP_D3_DiracFermion *in,
			QDP_D3_DiracFermion *r,
			QDP_Subset subset);

QOP_D_Real
QOP_D3_relnorm2_V(QDP_D3_ColorVector **rsd, 
		   QDP_D3_ColorVector **out, 
		   QDP_Subset subset, int nv);

QOP_D_Real
QOP_D3_relnorm2_D(QDP_D3_DiracFermion **rsd, 
		   QDP_D3_DiracFermion **out, 
		   QDP_Subset subset, int nv);

#if QOP_Precision == 'D'
#  if QOP_Colors == 3
#    include <qop_d3_internal_generic.h>
#  endif
#endif

#endif /* _QOP_D3_INTERNAL_H */
