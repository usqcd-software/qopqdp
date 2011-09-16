#ifndef _QOP_INTERNALP_H
#define _QOP_INTERNALP_H

#if QOP_Precision == 1
#define QOPP(x) QOP_F_##x
#define QOPPC(x) QOP_F3_##x
#define QDPP(x) QDP_F_##x
#define QDPPC(x) QDP_F3_##x
#define REAL float
#else
#define QOPP(x) QOP_D_##x
#define QOPPC(x) QOP_D3_##x
#define QDPP(x) QDP_D_##x
#define QDPPC(x) QDP_D3_##x
#define REAL double
#endif

//AB HISQ derivatives, rank-4 tensor
#if QOP_Precision == 1
#define QLA_ColorTensor4 QLA_F3_ColorTensor4
#else
#define QLA_ColorTensor4 QLA_D3_ColorTensor4
#endif


#define QLATYPE_V QLA_ColorVector
#define QLATYPE_D QLA_DiracFermion
#define QLATYPE_M QLA_ColorMatrix

#define QOP_qdp_eq_raw(abbr, qdp, raw, evenodd) \
  QDPPC(insert_##abbr)(qdp, (QLATYPE_##abbr *)raw, QDP_all)

#define QOP_raw_eq_qdp(abbr, raw, qdp, evenodd)		      \
  if(evenodd==QOP_EVEN) {				      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_even); \
  } else if(evenodd==QOP_ODD) {				      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_odd);  \
  } else {						      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_all);  \
  }

typedef struct {
  int tmp;
} QOPPC(common_t);
extern QOPPC(common_t) QOPPC(common);

struct QOPPC(ColorVector_struct) {
  QDPPC(ColorVector) *cv;
  REAL *raw;
};

struct QOPPC(DiracFermion_struct) {
  QDPPC(DiracFermion) *df;
  REAL *raw;
};

struct QOPPC(GaugeField_struct) {
  QDPPC(ColorMatrix) **links;
  REAL **raw;
};

struct QOPPC(Force_struct) {
  QDPPC(ColorMatrix) **force;
  REAL **raw;
};

typedef struct {
  QDP_ColorVector **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOPPC(eigcg_t_V);

typedef struct {
  QDP_DiracFermion **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOPPC(eigcg_t_D);

  /* Asqtad datatypes */

struct QOPPC(FermionLinksAsqtad_struct) {
  int dblstored, nlinks;
  QDPPC(ColorMatrix) **fatlinks;
  QDPPC(ColorMatrix) **longlinks;
  QDPPC(ColorMatrix) **fwdlinks;
  QDPPC(ColorMatrix) **bcklinks;
  QDPPC(ColorMatrix) **dbllinks;
  QOPPC(eigcg_t_V) eigcg;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
};

  /* HISQ datatypes */

struct QOPPC(FermionLinksHisq_struct) {
  //  int dblstored, nlinks;
  int n_naiks, WeqY;
  //AB intermediate links
  QDPPC(ColorMatrix) **U_links; // gauge links
  QDPPC(ColorMatrix) **V_links; // Fat7 smeared
  QDPPC(ColorMatrix) **Y_unitlinks; // projected to U(3),
  QDPPC(ColorMatrix) **W_unitlinks; // projected to SU(3)
  // normally we project only to U(3) and W_unitlink is a pointer to Y_unitlink
  //AB actual array where extra index distinguishes
  //   different epsilon corrections to 1-link and Naik terms
  QOPPC(FermionLinksAsqtad) **fn;
  QOPPC(FermionLinksAsqtad) *fn_deps;
//  QOPPC(eigcg_t_V) eigcg;
//  QDP_Shift shifts[8];
//  QDP_Shift shifts_dbl[16];
//  QDP_ShiftDir shiftdirs_dbl[16];
};

  /* Wilson datatypes */

struct QOPPC(FermionLinksWilson_struct) {
  REAL clovinvkappa;
  int dblstored;
  QDPPC(ColorMatrix) **links;
  QDPPC(ColorMatrix) **bcklinks;
  QDPPC(ColorMatrix) **dbllinks;
  QOPPC(GaugeField) *qopgf;
  QDPPC(DiracPropagator) *qdpclov;
  REAL *clov, *clovinv;
  REAL **rawlinks, *rawclov;
  QOPPC(eigcg_t_D) eigcg;
};

  /* Domain Wall datatypes */

// Current DWF implementation explicitly calls Wilson op
struct QOPPC(FermionLinksDW_struct) {
  QOPPC(FermionLinksWilson) *flw;
};

/* internal routines */

QOPPC(FermionLinksAsqtad) *
  QOPPC(asqtad_create_L_from_L)(QOPPC(FermionLinksAsqtad) *fla_src);

QOPPC(FermionLinksAsqtad) *QOPPC(asqtad_create_L_from_r_times_L)(QLA_Real *s,
	      QOPPC(FermionLinksAsqtad) *fla_src);

void QOPPC(asqtad_L_peq_L)(QOPPC(FermionLinksAsqtad) *fla,
			   QOPPC(FermionLinksAsqtad) *fla1);

void QOPPC(qdpM_eq_raw)(QDP_ColorMatrix *cm, REAL *lnk);


typedef void (QOPPC(linop_t_V))(QDP_ColorVector *out, QDP_ColorVector *in,
				QDP_Subset subset);

typedef void (QOPPC(linop_t_D))(QDP_DiracFermion *out, QDP_DiracFermion *in,
				QDP_Subset subset);

typedef void (QOPPC(linop_t_vD))(QDP_DiracFermion **out, QDP_DiracFermion **in,
				 QDP_Subset subset);

QOP_status_t
QOPPC(invert_cg_V)(QOPPC(linop_t_V) *linop,
		   QOP_invert_arg_t *inv_arg,
		   QOP_resid_arg_t *res_arg,
		   QDP_ColorVector *out,
		   QDP_ColorVector *in,
		   QDP_ColorVector *p,
		   QDP_Subset subset);

QOP_status_t
QOPPC(invert_cg_D)(QOPPC(linop_t_D) *linop,
		   QOP_invert_arg_t *inv_arg,
		   QOP_resid_arg_t *res_arg,
		   QDP_DiracFermion *out,
		   QDP_DiracFermion *in,
		   QDP_DiracFermion *p,
		   QDP_Subset subset);

QOP_status_t
QOPPC(invert_cg_vD)(QOPPC(linop_t_vD) *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t *res_arg,
		    QDP_DiracFermion **out,
		    QDP_DiracFermion **in,
		    QDP_DiracFermion **p,
		    QDP_Subset subset,
		    int _n);

QOP_status_t
QOPPC(invert_cgms_V)(QOPPC(linop_t_V) *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t **res_arg,
		     QLA_Real *shifts,
		     int nshifts,
		     QDP_ColorVector **out,
		     QDP_ColorVector *in,
		     QDP_ColorVector *p,
		     QDP_Subset subset);

QOP_status_t
QOPPC(invert_cgms_D)(QOPPC(linop_t_D) *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t **res_arg,
		     QLA_Real *shifts,
		     int nshifts,
		     QDP_DiracFermion **out,
		     QDP_DiracFermion *in,
		     QDP_DiracFermion *p,
		     QDP_Subset subset);

QOP_status_t
QOPPC(invert_cgms_vD)(QOPPC(linop_t_vD) *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t **res_arg,
		      QLA_Real *shifts,
		      int nshifts,
		      QDP_DiracFermion ***out,
		      QDP_DiracFermion **in,
		      QDP_DiracFermion **p,
		      QDP_Subset subset,
		      int _n);

QOP_status_t
QOPPC(invert_bicgstab_D)(QOPPC(linop_t_D) *linop,
			 QOP_invert_arg_t *inv_arg,
			 QOP_resid_arg_t *res_arg,
			 QDP_DiracFermion *out,
			 QDP_DiracFermion *in,
			 QDP_DiracFermion *p,
			 QDP_DiracFermion *r,
			 QDP_Subset subset);

QOP_status_t
QOPPC(invert_eigcg_V)(QOPPC(linop_t_V) *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_ColorVector *out,
		      QDP_ColorVector *in,
		      QDP_ColorVector *p,
		      QDP_Subset subset,
		      QOPPC(eigcg_t_V) *eigcg);

QOP_status_t
QOPPC(invert_eigcg_D)(QOPPC(linop_t_D) *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      QDP_DiracFermion *out,
		      QDP_DiracFermion *in,
		      QDP_DiracFermion *p,
		      QDP_Subset subset,
		      QOPPC(eigcg_t_D) *eigcg);

//AB internal operations for HISQ

void
QOPPC(hisq_force_multi_reunit)(QOP_info_t *info,
				     QDP_ColorMatrix *gf[4],
				     QDP_ColorMatrix *force_accum[4],
				     QDP_ColorMatrix *force_accum_old[4]);


void 
QOPPC(hisq_force_multi_wrapper_fnmat)(QOP_info_t *info,  
				      QOP_FermionLinksHisq *flh,
				      QOP_Force *Force, 
				      QOP_hisq_coeffs_t *hisq_coeff,
				      REAL *epsv,
				      QDP_ColorVector *in_pt[], 
				      int *n_orders_naik);

void 
QOPPC(hisq_force_multi_wrapper_fnmat2)(QOP_info_t *info,  
				       QOP_FermionLinksHisq *flh,
				       QOP_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       REAL *epsv,
				       QDP_ColorVector *in_pt[], 
				       int *n_orders_naik);

void
QOPPC(u3reunit)(QDP_ColorMatrix *U, QDP_ColorMatrix *V);

void
QOPPC(su3reunit)(QDP_ColorMatrix *U,QDP_ColorMatrix *Ur);

QOP_status_t
QOPPC(invert_gcr2_D)(QOPPC(linop_t_D) *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     QDP_DiracFermion *out,
		     QDP_DiracFermion *in,
		     QDP_DiracFermion *r,
		     QDP_Subset subset);

QOP_status_t
QOPPC(invert_gmres2_D)(QOPPC(linop_t_D) *linop,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QDP_DiracFermion *out,
		       QDP_DiracFermion *in,
		       QDP_DiracFermion *r,
		       QDP_Subset subset);

QLA_Real
QOPPC(relnorm2_V)(QDP_ColorVector **rsd, 
		  QDP_ColorVector **out, 
		  QDP_Subset subset, int nv);

QLA_Real
QOPPC(relnorm2_D)(QDP_DiracFermion **rsd, 
		  QDP_DiracFermion **out, 
		  QDP_Subset subset, int nv);


#endif /* _QOP_INTERNALP_H */
