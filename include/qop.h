#ifndef _QOP_H
#define _QOP_H

#ifdef __cplusplus
extern "C" {
#endif

  /* Enumerate in order of increasing verbosity */
#define QOP_VERB_OFF    0
#define QOP_VERB_LOW    1
#define QOP_VERB_MED    2
#define QOP_VERB_HI     3
#define QOP_VERB_DEBUG  4

/* Maximum number of Naik terms */
/* Must match MAX_NAIK in MILC code */
#define QOP_MAX_NAIK 10

typedef enum {
  QOP_SUCCESS = 0,
  QOP_FAIL = 1,
  QOP_BC_ERROR,  /* inconsistent boundary condition */
  QOP_MEM_ERROR, /* memory errors such as failure to allocate, etc */
  QOP_TOP_ERROR  /* to be used as the counter for numeber of errors? */
} QOP_status_t;

typedef enum {
  QOP_EVEN = 0,
  QOP_ODD = 1,
  QOP_EVENODD = 2
} QOP_evenodd_t;

typedef enum {
  QOP_UNITARIZE_U3 = 0,
  QOP_UNITARIZE_SU3 = 1
} QOP_hisq_unitarize_group_t;

/* So far, only rational and none are supported */
typedef enum {
  QOP_UNITARIZE_NONE = 0,
  //  QOP_UNITARIZE_APE = 1,
  //  QOP_UNITARIZE_ROOT = 2,
  QOP_UNITARIZE_RATIONAL = 3,
  //  QOP_UNITARIZE_HISQ = 4,
  QOP_UNITARIZE_ANALYTIC = 5,
  //  QOP_UNITARIZE_STOUT = 6
} QOP_hisq_unitarize_method_t;

typedef struct {
  int (*node_number)(const int coords[]); /* node no for given latt coord */
  int (*node_index)(const int coords[]);  /* site rank for given latt coord */
  int latdim;                             /* number of lattice dimensions */
  int *latsize;                           /* physical lattice lengths */
  int machdim;                            /* number of logical machine dims */
  int *machsize;                          /* logical grid lengths */
  int this_node;                          /* lexicographic node number */
  int sites_on_node;
} QOP_layout_t;

typedef struct {
  double re;
  double im;
} QOP_Complex;

  /* User must allocate a length 'latdim' array for the boundary phase in
     each direction.  These will multiply the last link in each direction. */
typedef struct {
  QOP_Complex *phase;
} QOP_bc_t;

  /* User must allocate a length 'latdim' array for the staggered sign on
     the links in each direction.  The int field is a bitmask where a 1
     bit means that the phase will depend on that direction.
     0 = no phase, 1 = (-1)^(x0), 5 = (-1)^(x0+x2), etc. */
typedef struct {
  int *signmask;
} QOP_staggered_sign_t;

typedef struct {
  char *tag;
  double value;
  void *extra;
} QOP_opt_t;

typedef struct {
  double final_sec;        /* (out) total number of seconds used */
  double final_flop;       /* (out) total number of flops performed */
  QOP_status_t status;     /* (out) error status */
  int count1, count2;      /* (out) generic counters */
} QOP_info_t;

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  int max_iter;            /* (in) max number of iterations */
  int restart;             /* (in) number of iterations before restart */
  int max_restarts;        /* (in) number of restarts allowed */
  QOP_evenodd_t evenodd;   /* (in) subset of source vector */
} QOP_invert_arg_t;
#define QOP_INVERT_ARG_DEFAULT ((QOP_invert_arg_t){2000,1000,5,QOP_EVENODD})

  /* these are quantities that vary for each mass in the multi inverter */
typedef struct {
  double rsqmin;           /* (in) desired squared residual. Ignored if 0. */
  double final_rsq;        /* (out) actual squared residual */
  double relmin;           /* (in) desired squared relative norm Ignored if 0. */
  double final_rel;        /* (out) actual squared relative norm. */
  int final_iter;          /* (out) number of iterations done */
  int final_restart;       /* (out) number of restarts done */
} QOP_resid_arg_t;
#define QOP_RESID_ARG_DEFAULT ((QOP_resid_arg_t){1e-6,0,0,0,0,0})

typedef struct QOP_F3_ColorVector_struct  QOP_F3_ColorVector;
typedef struct QOP_F3_DiracFermion_struct QOP_F3_DiracFermion;
typedef struct QOP_F3_GaugeField_struct   QOP_F3_GaugeField;
typedef struct QOP_F3_Force_struct        QOP_F3_Force;

typedef struct QOP_D3_ColorVector_struct  QOP_D3_ColorVector;
typedef struct QOP_D3_DiracFermion_struct QOP_D3_DiracFermion;
typedef struct QOP_D3_GaugeField_struct   QOP_D3_GaugeField;
typedef struct QOP_D3_Force_struct        QOP_D3_Force;

  /* Improved gauge datatypes */

typedef struct {
   double plaquette;
   double rectangle;
   double parallelogram;
} QOP_gauge_coeffs_t;

  /* Asqtad datatypes */
  /* to follow the convetion defined in Orginos, et al (PRD 60, 045403) */
typedef struct {
  double one_link;
  double three_staple;
  double five_staple;
  double seven_staple;
  double lepage;
  double naik;
} QOP_asqtad_coeffs_t;

typedef struct QOP_F3_FermionLinksAsqtad_struct QOP_F3_FermionLinksAsqtad;
typedef struct QOP_D3_FermionLinksAsqtad_struct QOP_D3_FermionLinksAsqtad;

  /* HISQ datatypes*/
typedef struct {
  int n_naiks;
  double eps_naik[QOP_MAX_NAIK];
  QOP_hisq_unitarize_group_t ugroup;
  QOP_hisq_unitarize_method_t umethod;
  double fat7_one_link;
  double fat7_three_staple;
  double fat7_five_staple;
  double fat7_seven_staple;
  double fat7_lepage;
  double asqtad_one_link;
  double asqtad_three_staple;
  double asqtad_five_staple;
  double asqtad_seven_staple;
  double asqtad_lepage;
  double asqtad_naik;
  double difference_one_link;
  double difference_naik;
} QOP_hisq_coeffs_t;

typedef struct QOP_F3_FermionLinksHisq_struct QOP_F3_FermionLinksHisq;
typedef struct QOP_D3_FermionLinksHisq_struct QOP_D3_FermionLinksHisq;


  /* Wilson datatypes */

typedef struct {
  double clov_s;
  double clov_t;
  double aniso;
} QOP_wilson_coeffs_t;

typedef struct QOP_F3_FermionLinksWilson_struct QOP_F3_FermionLinksWilson;
typedef struct QOP_D3_FermionLinksWilson_struct QOP_D3_FermionLinksWilson;

  /* Domain Wall datatypes */

typedef struct {
  // need to add Mobius parameters here
} QOP_dw_coeffs_t;

typedef struct QOP_F3_FermionLinksDW_struct QOP_F3_FermionLinksDW;
typedef struct QOP_D3_FermionLinksDW_struct QOP_D3_FermionLinksDW;


  /**********************/
  /*  General routines  */
  /**********************/

QOP_status_t QOP_init(QOP_layout_t *layout);
QOP_status_t QOP_finalize(void);
int QOP_is_initialized(void);

int QOP_verbose(int level);
int QOP_profcontrol(int level);

int QOP_node_number_raw (int coords[]);
int QOP_node_index_raw_V(int coords[], QOP_evenodd_t evenodd);
int QOP_node_index_raw_D(int coords[], QOP_evenodd_t evenodd);
int QOP_node_index_raw_G(int coords[], QOP_evenodd_t evenodd);
int QOP_node_index_raw_F(int coords[], QOP_evenodd_t evenodd);
int QOP_sites_on_node_raw_V(QOP_evenodd_t evenodd);
int QOP_sites_on_node_raw_D(QOP_evenodd_t evenodd);
int QOP_sites_on_node_raw_G(QOP_evenodd_t evenodd);
int QOP_sites_on_node_raw_F(QOP_evenodd_t evenodd);

#define QOP_qla_type_F3_V QLA_F3_ColorVector
#define QOP_qla_type_F3_D QLA_F3_DiracFermion
#define QOP_qla_type_F3_M QLA_F3_ColorMatrix
#define QOP_qla_type_D3_V QLA_D3_ColorVector
#define QOP_qla_type_D3_D QLA_D3_DiracFermion
#define QOP_qla_type_D3_M QLA_D3_ColorMatrix
#define QOP_raw_size(P, T) (QDP_sites_on_node*sizeof(QOP_qla_type_##P##3_##T))
#define QOP_F3_raw_size_V(evenodd) QOP_raw_size(F, V)
#define QOP_F3_raw_size_D(evenodd) QOP_raw_size(F, D)
#define QOP_F3_raw_size_G(evenodd) QOP_raw_size(F, M)
#define QOP_F3_raw_size_F(evenodd) QOP_raw_size(F, M)
#define QOP_D3_raw_size_V(evenodd) QOP_raw_size(D, V)
#define QOP_D3_raw_size_D(evenodd) QOP_raw_size(D, D)
#define QOP_D3_raw_size_G(evenodd) QOP_raw_size(D, M)
#define QOP_D3_raw_size_F(evenodd) QOP_raw_size(D, M)

#define QOP_elem(P, T, raw, i, ...) QLA_##P##3_elem_##T(((QOP_qla_type##P##3_##T *)raw)[i], __ARGV__)
#define QOP_set(P, T, raw, i, re, im, ...) QLA_c_eq_r_plus_ir(QOP_elem(P, T, raw, i, __ARGV__), re, im)
#define QOP_F3_raw_set_V(raw, evenodd, i, ic, re, im) QOP_set(F, V, raw, i, re, im, ic)
#define QOP_F3_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_set(F, D, raw, i, re, im, ic, is)
#define QOP_F3_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_set(F, M, raw, i, re, im, ic, jc)
#define QOP_F3_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_set(F, M, raw, i, re, im, ic, jc)
#define QOP_D3_raw_set_V(raw, evenodd, i, ic, re, im) QOP_set(D, V, raw, i, re, im, ic)
#define QOP_D3_raw_set_D(raw, evenodd, i, ic, is, re, im) QOP_set(D, D, raw, i, re, im, ic, is)
#define QOP_D3_raw_set_G(raw, evenodd, i, ic, jc, re, im) QOP_set(D, M, raw, i, re, im, ic, jc)
#define QOP_D3_raw_set_F(raw, evenodd, i, ic, jc, re, im) QOP_set(D, M, raw, i, re, im, ic, jc)
#define QOP_get(P, T, re, im, raw, i, ...) { QLA_##P##_Complex _c = QOP_elem(P, T, raw, i, __ARGV__); re = QLA_real(_c); im = QLA_imag(_c); }
#define QOP_F3_raw_get_V(re, im, raw, evenodd, i, ic) QOP_set(F, V, re, im, raw, i, ic)
#define QOP_F3_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_set(F, D, re, im, raw, i, ic, is)
#define QOP_F3_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_set(F, M, re, im, raw, i, ic, jc)
#define QOP_F3_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_set(F, M, re, im, raw, i, ic, jc)
#define QOP_D3_raw_get_V(re, im, raw, evenodd, i, ic) QOP_set(D, V, re, im, raw, i, ic)
#define QOP_D3_raw_get_D(re, im, raw, evenodd, i, ic, is) QOP_set(D, D, re, im, raw, i, ic, is)
#define QOP_D3_raw_get_M(re, im, raw, evenodd, i, ic, jc) QOP_set(D, M, re, im, raw, i, ic, jc)
#define QOP_D3_raw_get_G(re, im, raw, evenodd, i, ic, jc) QOP_set(D, M, re, im, raw, i, ic, jc)


  /* single precision */

/* create a QOP field with a copy of the raw source field */
QOP_F3_ColorVector  *QOP_F3_create_V_from_raw(float *src,
					      QOP_evenodd_t evenodd);
QOP_F3_DiracFermion *QOP_F3_create_D_from_raw(float *src,
					      QOP_evenodd_t evenodd);
QOP_F3_GaugeField   *QOP_F3_create_G_from_raw(float *links[],
					      QOP_evenodd_t evenodd);
QOP_F3_Force        *QOP_F3_create_F_from_raw(float *force[],
					      QOP_evenodd_t evenodd);

/* copy QOP field into a raw field */
void QOP_F3_extract_V_to_raw(float *dest, QOP_F3_ColorVector *src,
			     QOP_evenodd_t evenodd);
void QOP_F3_extract_D_to_raw(float *dest, QOP_F3_DiracFermion *src,
			     QOP_evenodd_t evenodd);
void QOP_F3_extract_G_to_raw(float *dest[], QOP_F3_GaugeField *src,
			     QOP_evenodd_t evenodd);
void QOP_F3_extract_F_to_raw(float *dest[], QOP_F3_Force *src,
			     QOP_evenodd_t evenodd);

/* destroy a QOP field */
/* if the QOP field was created with a convert from raw function then
   the user must still free the original raw field themself */
void QOP_F3_destroy_V(QOP_F3_ColorVector *field);
void QOP_F3_destroy_D(QOP_F3_DiracFermion *field);
void QOP_F3_destroy_G(QOP_F3_GaugeField *field);
void QOP_F3_destroy_F(QOP_F3_Force *field);

/* create a QOP field using the raw source field */
/* the raw source is not freed and the user must not change or free it until
   the QOP field has been converted back to raw or destroyed */
QOP_F3_ColorVector  *QOP_F3_convert_V_from_raw(float *src,
					       QOP_evenodd_t evenodd);
QOP_F3_DiracFermion *QOP_F3_convert_D_from_raw(float *src,
					       QOP_evenodd_t evenodd);
QOP_F3_GaugeField   *QOP_F3_convert_G_from_raw(float *links[],
					       QOP_evenodd_t evenodd);
QOP_F3_Force        *QOP_F3_convert_F_from_raw(float *force[],
					       QOP_evenodd_t evenodd);

/* create a raw field from the data in the QOP field and destroy it */
/* if the QOP field was created with a convert from raw function then
   this will return the same raw source used as input */
float  *QOP_F3_convert_V_to_raw(QOP_F3_ColorVector *src,
				QOP_evenodd_t evenodd);
float  *QOP_F3_convert_D_to_raw(QOP_F3_DiracFermion *src,
				QOP_evenodd_t evenodd);
float **QOP_F3_convert_G_to_raw(QOP_F3_GaugeField *src,
				QOP_evenodd_t evenodd);
float **QOP_F3_convert_F_to_raw(QOP_F3_Force *src,
				QOP_evenodd_t evenodd);

  /* double precision */

QOP_D3_ColorVector  *QOP_D3_create_V_from_raw(double *src,
					      QOP_evenodd_t evenodd);
QOP_D3_DiracFermion *QOP_D3_create_D_from_raw(double *src,
					      QOP_evenodd_t evenodd);
QOP_D3_GaugeField   *QOP_D3_create_G_from_raw(double *links[],
					      QOP_evenodd_t evenodd);
QOP_D3_Force        *QOP_D3_create_F_from_raw(double *force[],
					      QOP_evenodd_t evenodd);

void QOP_D3_extract_V_to_raw(double *dest, QOP_D3_ColorVector *src,
			     QOP_evenodd_t evenodd);
void QOP_D3_extract_D_to_raw(double *dest, QOP_D3_DiracFermion *src,
			     QOP_evenodd_t evenodd);
void QOP_D3_extract_G_to_raw(double *dest[], QOP_D3_GaugeField *src,
			     QOP_evenodd_t evenodd);
void QOP_D3_extract_F_to_raw(double *dest[], QOP_D3_Force *src,
			     QOP_evenodd_t evenodd);

void QOP_D3_destroy_V(QOP_D3_ColorVector *field);
void QOP_D3_destroy_D(QOP_D3_DiracFermion *field);
void QOP_D3_destroy_G(QOP_D3_GaugeField *field);
void QOP_D3_destroy_F(QOP_D3_Force *field);

QOP_D3_ColorVector  *QOP_D3_convert_V_from_raw(double *src,
					       QOP_evenodd_t evenodd);
QOP_D3_DiracFermion *QOP_D3_convert_D_from_raw(double *src,
					       QOP_evenodd_t evenodd);
QOP_D3_GaugeField   *QOP_D3_convert_G_from_raw(double *links[],
					       QOP_evenodd_t evenodd);
QOP_D3_Force        *QOP_D3_convert_F_from_raw(double *force[],
					       QOP_evenodd_t evenodd);

double  *QOP_D3_convert_V_to_raw(QOP_D3_ColorVector *src,
				 QOP_evenodd_t evenodd);
double  *QOP_D3_convert_D_to_raw(QOP_D3_DiracFermion *src,
				 QOP_evenodd_t evenodd);
double **QOP_D3_convert_G_to_raw(QOP_D3_GaugeField *src,
				 QOP_evenodd_t evenodd);
double **QOP_D3_convert_F_to_raw(QOP_D3_Force *src,
				 QOP_evenodd_t evenodd);

  /* puts in boundary condition and staggered phases in place */
  /* if either bc or ksphase is NULL those phases are ignored */
  /* see the corresponding structure definitions for conventions */

void QOP_F3_rephase_G(QOP_F3_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_D3_rephase_G(QOP_D3_GaugeField *links,
		      int *r0,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);


  /********************/
  /*  Gauge routines  */
  /********************/

void QOP_F3_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_F3_GaugeField *gauge,
					float *acts, float *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_D3_symanzik_1loop_gauge_action(QOP_info_t *info,
					QOP_D3_GaugeField *gauge,
					double *acts, double *actt,
					QOP_gauge_coeffs_t *coeffs);

void QOP_F3_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_F3_GaugeField *gauge, 
				       QOP_F3_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       float eps);

void QOP_D3_symanzik_1loop_gauge_force(QOP_info_t *info, 
				       QOP_D3_GaugeField *gauge, 
				       QOP_D3_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       double eps);

void QOP_F3_symanzik_1loop_gauge_deriv(QOP_info_t *info, 
				       QOP_F3_GaugeField *gauge, 
				       QOP_F3_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       float eps);

void QOP_D3_symanzik_1loop_gauge_deriv(QOP_info_t *info, 
				       QOP_D3_GaugeField *gauge, 
				       QOP_D3_Force *force,
				       QOP_gauge_coeffs_t *coeffs,
				       double eps);


  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* fermion matrix link routines */

  /* single precision */

QOP_F3_FermionLinksAsqtad *
  QOP_F3_asqtad_create_L_from_raw(float *fatlinks[], float *longlinks[],
				  QOP_evenodd_t evenodd);

  /* CJ: This is smearing routine in effect */
QOP_F3_FermionLinksAsqtad *
  QOP_F3_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_F3_GaugeField *gauge);

void QOP_F3_asqtad_extract_L_to_raw(float *fatlinks[], float *longlinks[],
				    QOP_F3_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_F3_asqtad_destroy_L(QOP_F3_FermionLinksAsqtad *field);

QOP_F3_FermionLinksAsqtad *
  QOP_F3_asqtad_convert_L_from_raw(float *fatlinks[], float *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_F3_asqtad_convert_L_to_raw(float ***fatlinks, float ***longlinks,
				    QOP_F3_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_F3_asqtad_load_L_from_raw(QOP_F3_FermionLinksAsqtad *asqtad,
				   float *fatlinks[], float *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_F3_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_F3_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_F3_GaugeField *gauge);

void QOP_F3_asqtad_rephase_L(QOP_F3_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);


  /* double precision */

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_raw(double *fatlinks[], double *longlinks[],
				  QOP_evenodd_t evenodd);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_G(QOP_info_t *info,
				QOP_asqtad_coeffs_t *coeffs,
				QOP_D3_GaugeField *gauge);

void QOP_D3_asqtad_extract_L_to_raw(double *fatlinks[], double *longlinks[],
				    QOP_D3_FermionLinksAsqtad *src,
				    QOP_evenodd_t evenodd);

void QOP_D3_asqtad_destroy_L(QOP_D3_FermionLinksAsqtad *field);

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_convert_L_from_raw(double *fatlinks[], double *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_D3_asqtad_convert_L_to_raw(double ***fatlinks, double ***longlinks,
				    QOP_D3_FermionLinksAsqtad *,
				    QOP_evenodd_t evenodd);

void QOP_D3_asqtad_load_L_from_raw(QOP_D3_FermionLinksAsqtad *asqtad,
				   double *fatlinks[], double *longlinks[],
				   QOP_evenodd_t evenodd);

void QOP_D3_asqtad_load_L_from_G(QOP_info_t *info,
				 QOP_D3_FermionLinksAsqtad *asqtad,
				 QOP_asqtad_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

void QOP_D3_asqtad_rephase_L(QOP_D3_FermionLinksAsqtad *fla,
			     int *r0,
			     QOP_bc_t *bc,
			     QOP_staggered_sign_t *sign);

  /* inverter routines */

QOP_status_t QOP_asqtad_invert_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_asqtad_dslash(QOP_info_t *info,
			  QOP_F3_FermionLinksAsqtad *asqtad,
			  float mass,
			  QOP_F3_ColorVector *out,
			  QOP_F3_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_F3_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_F3_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_F3_ColorVector *out,
			      QOP_F3_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_F3_asqtad_diaginv(QOP_info_t *info,
			   QOP_F3_FermionLinksAsqtad *asqtad,
			   float mass,
			   QOP_F3_ColorVector *out,
			   QOP_F3_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_F3_asqtad_invert(QOP_info_t *info,
			  QOP_F3_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  float mass,
			  QOP_F3_ColorVector *out_pt,
			  QOP_F3_ColorVector *in_pt);

void QOP_F3_asqtad_invert_multi(QOP_info_t *info,
				QOP_F3_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				float *masses[],
				int nmass[],
				QOP_F3_ColorVector **out_pt[],
				QOP_F3_ColorVector *in_pt[],
				int nsrc);

void QOP_D3_asqtad_dslash(QOP_info_t *info,
			  QOP_D3_FermionLinksAsqtad *asqtad,
			  double mass,
			  QOP_D3_ColorVector *out,
			  QOP_D3_ColorVector *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_D3_asqtad_dslash_dir(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      int dir, int fb,
			      double wtfat, double wtlong,
			      QOP_D3_ColorVector *out,
			      QOP_D3_ColorVector *in,
			      QOP_evenodd_t eo_out);

void QOP_D3_asqtad_diaginv(QOP_info_t *info,
			   QOP_D3_FermionLinksAsqtad *asqtad,
			   double mass,
			   QOP_D3_ColorVector *out,
			   QOP_D3_ColorVector *in,
			   QOP_evenodd_t eo);

void QOP_D3_asqtad_invert(QOP_info_t *info,
			  QOP_D3_FermionLinksAsqtad *asqtad,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  double mass,
			  QOP_D3_ColorVector *out_pt,
			  QOP_D3_ColorVector *in_pt);

void QOP_D3_asqtad_invert_multi(QOP_info_t *info,
				QOP_D3_FermionLinksAsqtad *asqtad,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				double *masses[],
				int nmass[],
				QOP_D3_ColorVector **out_pt[],
				QOP_D3_ColorVector *in_pt[],
				int nsrc);

  /* fermion force routines */

QOP_status_t QOP_asqtad_force_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_asqtad_force(QOP_info_t *info,
			 QOP_F3_GaugeField *gauge,
			 QOP_F3_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 float eps,
			 QOP_F3_ColorVector *in_pt);

void QOP_F3_asqtad_force_multi(QOP_info_t *info,
			       QOP_F3_GaugeField *gauge,
			       QOP_F3_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       float eps[],
			       QOP_F3_ColorVector *in_pt[],
			       int nsrc);

void QOP_D3_asqtad_force(QOP_info_t *info,
			 QOP_D3_GaugeField *gauge,
			 QOP_D3_Force *force,
			 QOP_asqtad_coeffs_t *coeffs,
			 double eps,
			 QOP_D3_ColorVector *in_pt);

void QOP_D3_asqtad_force_multi(QOP_info_t *info,
			       QOP_D3_GaugeField *gauge,
			       QOP_D3_Force *force,
			       QOP_asqtad_coeffs_t *coef,
			       double eps[],
			       QOP_D3_ColorVector *in_pt[],
			       int nsrc);


  /*********************/
  /*  HISQ routines  */
  /*********************/

  /* QOP_info_t counters as they are used in HISQ routines */
  /* The user must initialize them and reset them */
#define QOP_info_hisq_svd_counter(info) ((info)->count1)
#define QOP_info_hisq_force_filter_counter(info) ((info)->count2) 

  /* fermion matrix link routines */

QOP_status_t QOP_hisq_links_set_opts(QOP_opt_t opts[], int nopts);

  /* single precision */

QOP_F3_FermionLinksHisq *
  QOP_F3_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_F3_GaugeField *gauge);

void QOP_F3_hisq_destroy_L(QOP_F3_FermionLinksHisq *field);

QOP_F3_FermionLinksAsqtad **
  QOP_F3_get_asqtad_links_from_hisq(QOP_F3_FermionLinksHisq *hl);
  
QOP_F3_FermionLinksAsqtad *
  QOP_F3_get_asqtad_deps_links_from_hisq(QOP_F3_FermionLinksHisq *hl);

  /* double precision */

QOP_D3_FermionLinksHisq *
  QOP_D3_hisq_create_L_from_G(QOP_info_t *info,
			      QOP_hisq_coeffs_t *coeffs,
			      QOP_D3_GaugeField *gauge);

void QOP_D3_hisq_destroy_L(QOP_D3_FermionLinksHisq *field);


QOP_D3_FermionLinksAsqtad **
  QOP_D3_get_asqtad_links_from_hisq(QOP_D3_FermionLinksHisq *hl);
  
QOP_D3_FermionLinksAsqtad *
  QOP_D3_get_asqtad_deps_links_from_hisq(QOP_D3_FermionLinksHisq *hl);


  /* fermion force routines */

QOP_status_t QOP_hisq_force_set_opts(QOP_opt_t opts[], int nopts);


void QOP_F3_hisq_force_multi(QOP_info_t *info,
			     QOP_F3_FermionLinksHisq *flh,
			     QOP_F3_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     float eps[],
			     QOP_F3_ColorVector *in_pt[],
			     int *n_orders_naik);
  
void QOP_D3_hisq_force_multi(QOP_info_t *info,
			     QOP_D3_FermionLinksHisq *flh,
			     QOP_D3_Force *force,
			     QOP_hisq_coeffs_t *coef,
			     double eps[],
			     QOP_D3_ColorVector *in_pt[],
			     int *n_orders_naik);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_raw(float *links[], float *clov,
				  QOP_evenodd_t evenodd);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_F3_GaugeField *gauge);

void QOP_F3_wilson_extract_L_to_raw(float *links[], float *clov,
				    QOP_F3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_F3_wilson_destroy_L(QOP_F3_FermionLinksWilson *field);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_raw(float *links[], float *clov,
				   QOP_evenodd_t evenodd);

void QOP_F3_wilson_convert_L_to_raw(float ***links, float **clov,
				    QOP_F3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_F3_GaugeField *gauge);

QOP_F3_GaugeField *
  QOP_F3_wilson_convert_L_to_G(QOP_F3_FermionLinksWilson *links);

void QOP_F3_wilson_load_L_from_raw(QOP_F3_FermionLinksWilson *wilson,
				   float *links[], float *clov,
				   QOP_evenodd_t evenodd);

void QOP_F3_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_F3_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_F3_GaugeField *gauge);

  /* double precision */

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_raw(double *links[], double *clov,
				  QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_D3_GaugeField *gauge);

void QOP_D3_wilson_extract_L_to_raw(double *links[], double *clov,
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_D3_wilson_destroy_L(QOP_D3_FermionLinksWilson *field);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_raw(double *links[], double *clov,
				   QOP_evenodd_t evenodd);

void QOP_D3_wilson_convert_L_to_raw(double ***links, double **clov,
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_wilson_convert_L_to_G(QOP_D3_FermionLinksWilson *links);

void QOP_D3_wilson_load_L_from_raw(QOP_D3_FermionLinksWilson *wilson,
				   double *links[], double *clov,
				   QOP_evenodd_t evenodd);

void QOP_D3_wilson_load_L_from_G(QOP_info_t *info,
				 QOP_D3_FermionLinksWilson *wilson,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

  /* create single-precision QOP objects from double-precision */

QOP_F3_GaugeField *
QOP_FD3_create_G_from_G(QOP_D3_GaugeField *qopgf_double);

QOP_F3_FermionLinksWilson *
QOP_FD3_wilson_create_L_from_L(QOP_D3_FermionLinksWilson *flw_double);

QOP_F3_FermionLinksAsqtad *
QOP_FD3_asqtad_create_L_from_L(QOP_D3_FermionLinksAsqtad *fla_src);

  /* inverter routines */

QOP_status_t QOP_wilson_invert_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_wilson_dslash(QOP_info_t *info,
			  QOP_F3_FermionLinksWilson *flw,
			  float kappa,
			  int sign,
			  QOP_F3_DiracFermion *out,
			  QOP_F3_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_F3_wilson_diaginv(QOP_info_t *info,
			   QOP_F3_FermionLinksWilson *flw,
			   float kappa,
			   QOP_F3_DiracFermion *out,
			   QOP_F3_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_F3_wilson_invert(QOP_info_t *info,
			  QOP_F3_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  float kappa,
			  QOP_F3_DiracFermion *out_pt,
			  QOP_F3_DiracFermion *in_pt);

void QOP_F3_wilson_invert_multi(QOP_info_t *info,
				QOP_F3_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				float *kappas[],
				int nkappa[],
				QOP_F3_DiracFermion **out_pt[],
				QOP_F3_DiracFermion *in_pt[],
				int nsrc);

void QOP_D3_wilson_dslash(QOP_info_t *info,
			  QOP_D3_FermionLinksWilson *flw,
			  double kappa,
			  int sign,
			  QOP_D3_DiracFermion *out,
			  QOP_D3_DiracFermion *in,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);

void QOP_D3_wilson_diaginv(QOP_info_t *info,
			   QOP_D3_FermionLinksWilson *flw,
			   double kappa,
			   QOP_D3_DiracFermion *out,
			   QOP_D3_DiracFermion *in,
			   QOP_evenodd_t eo);

void QOP_D3_wilson_invert(QOP_info_t *info,
			  QOP_D3_FermionLinksWilson *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  double kappa,
			  QOP_D3_DiracFermion *out_pt,
			  QOP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_invert_multi(QOP_info_t *info,
				QOP_D3_FermionLinksWilson *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				double *kappas[],
				int nkappa[],
				QOP_D3_DiracFermion **out_pt[],
				QOP_D3_DiracFermion *in_pt[],
				int nsrc);

  /* fermion force routines */

QOP_status_t QOP_wilson_force_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_wilson_force(QOP_info_t *info,
			 QOP_F3_GaugeField *gauge,
			 QOP_F3_Force *force,
			 QOP_wilson_coeffs_t *coeffs,
			 float eps,
			 QOP_F3_DiracFermion *in_pt);

void QOP_F3_wilson_force_multi(QOP_info_t *info,
			       QOP_F3_GaugeField *gauge,
			       QOP_F3_Force *force,
			       QOP_wilson_coeffs_t *coef,
			       float eps[],
			       QOP_F3_DiracFermion *in_pt[],
			       int nsrc);

void QOP_D3_wilson_force(QOP_info_t *info,
			 QOP_D3_GaugeField *gauge,
			 QOP_D3_Force *force,
			 QOP_wilson_coeffs_t *coeffs,
			 double eps,
			 QOP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_force_multi(QOP_info_t *info,
			       QOP_D3_GaugeField *gauge,
			       QOP_D3_Force *force,
			       QOP_wilson_coeffs_t *coef,
			       double eps[],
			       QOP_D3_DiracFermion *in_pt[],
			       int nsrc);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

  /* fermion matrix link routines */

QOP_F3_FermionLinksDW *
  QOP_F3_dw_create_L_from_raw(float *links[], QOP_evenodd_t evenodd);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_F3_GaugeField *gauge);

void QOP_F3_dw_extract_L_to_raw(float *links[],
				QOP_F3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_F3_dw_destroy_L(QOP_F3_FermionLinksDW *field);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_convert_L_from_raw(float *links[], QOP_evenodd_t evenodd);

void QOP_F3_dw_convert_L_to_raw(float ***links,
				QOP_F3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_F3_GaugeField *gauge);

QOP_F3_GaugeField *
  QOP_F3_dw_convert_L_to_G(QOP_F3_FermionLinksDW *links);

void QOP_F3_dw_load_L_from_raw(QOP_F3_FermionLinksDW *dw,
			       float *links[], QOP_evenodd_t evenodd);

void QOP_F3_dw_load_L_from_G(QOP_info_t *info,
			     QOP_F3_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_F3_GaugeField *gauge);


  /* double precision */

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_raw(double *links[], QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_D3_GaugeField *gauge);

void QOP_D3_dw_extract_L_to_raw(double *links[],
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_D3_dw_destroy_L(QOP_D3_FermionLinksDW *field);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_raw(double *links[], QOP_evenodd_t evenodd);

void QOP_D3_dw_convert_L_to_raw(double ***links,
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_dw_convert_L_to_G(QOP_D3_FermionLinksDW *links);

void QOP_D3_dw_load_L_from_raw(QOP_D3_FermionLinksDW *dw,
			       double *links[],
			       QOP_evenodd_t evenodd);

void QOP_D3_dw_load_L_from_G(QOP_info_t *info,
			     QOP_D3_FermionLinksDW *dw,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_D3_GaugeField *gauge);

  /* inverter routines */

QOP_status_t QOP_dw_invert_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_dw_dslash(QOP_info_t *info,
		      QOP_F3_FermionLinksDW *links,
		      float M5,
		      float m,
		      int sign,
		      QOP_F3_DiracFermion *out_pt[],
		      QOP_F3_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_F3_dw_invert(QOP_info_t *info,
		      QOP_F3_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      float M5,
		      float m,
		      QOP_F3_DiracFermion *out_pt[],
		      QOP_F3_DiracFermion *in_pt[],
		      int Ls);

void QOP_F3_dw_invert_multi(QOP_info_t *info,
			    QOP_F3_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    float *M5[],
			    float *m[],
			    int nmass[],
			    QOP_F3_DiracFermion ***out_pt[],
			    QOP_F3_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

void QOP_D3_dw_dslash(QOP_info_t *info,
		      QOP_D3_FermionLinksDW *links,
		      double M5,
		      double m,
		      int sign,
		      QOP_D3_DiracFermion *out_pt[],
		      QOP_D3_DiracFermion *in_pt[],
		      int Ls,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in);

void QOP_D3_dw_invert(QOP_info_t *info,
		      QOP_D3_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      double M5,
		      double m,
		      QOP_D3_DiracFermion *out_pt[],
		      QOP_D3_DiracFermion *in_pt[],
		      int Ls);

void QOP_D3_dw_invert_multi(QOP_info_t *info,
			    QOP_D3_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    double *M5[],
			    double *m[],
			    int nmass[],
			    QOP_D3_DiracFermion ***out_pt[],
			    QOP_D3_DiracFermion **in_pt[],
			    int nsrc,
			    int Ls);

 /* fermion force routines */

QOP_status_t QOP_dw_force_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_dw_force(QOP_info_t *info,
		     QOP_F3_GaugeField *gauge,
		     QOP_F3_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     float eps,
		     QOP_F3_DiracFermion *in_pt);

void QOP_F3_dw_force_multi(QOP_info_t *info,
			   QOP_F3_GaugeField *gauge,
			   QOP_F3_Force *force,
			   QOP_dw_coeffs_t *coef,
			   float eps[],
			   QOP_F3_DiracFermion *in_pt[],
			   int nsrc);

void QOP_D3_dw_force(QOP_info_t *info,
		     QOP_D3_GaugeField *gauge,
		     QOP_D3_Force *force,
		     QOP_dw_coeffs_t *coeffs,
		     double eps,
		     QOP_D3_DiracFermion *in_pt);

void QOP_D3_dw_force_multi(QOP_info_t *info,
			   QOP_D3_GaugeField *gauge,
			   QOP_D3_Force *force,
			   QOP_dw_coeffs_t *coef,
			   double eps[],
			   QOP_D3_DiracFermion *in_pt[],
			   int nsrc);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#ifndef QOP_Precision
#define QOP_Precision 1
#endif

#if QOP_Precision == 1

#define QOP_ColorVector  QOP_F3_ColorVector
#define QOP_DiracFermion QOP_F3_DiracFermion
#define QOP_GaugeField   QOP_F3_GaugeField
#define QOP_Force        QOP_F3_Force

#define QOP_create_V_from_raw QOP_F3_create_V_from_raw
#define QOP_create_D_from_raw QOP_F3_create_D_from_raw
#define QOP_create_G_from_raw QOP_F3_create_G_from_raw
#define QOP_create_F_from_raw QOP_F3_create_F_from_raw

#define QOP_extract_V_to_raw QOP_F3_extract_V_to_raw
#define QOP_extract_D_to_raw QOP_F3_extract_D_to_raw
#define QOP_extract_G_to_raw QOP_F3_extract_G_to_raw
#define QOP_extract_F_to_raw QOP_F3_extract_F_to_raw

#define QOP_destroy_V QOP_F3_destroy_V
#define QOP_destroy_D QOP_F3_destroy_D
#define QOP_destroy_G QOP_F3_destroy_G
#define QOP_destroy_F QOP_F3_destroy_F

#define QOP_convert_V_from_raw QOP_F3_convert_V_from_raw
#define QOP_convert_D_from_raw QOP_F3_convert_D_from_raw
#define QOP_convert_G_from_raw QOP_F3_convert_G_from_raw
#define QOP_convert_F_from_raw QOP_F3_convert_F_from_raw

#define QOP_convert_V_to_raw QOP_F3_convert_V_to_raw
#define QOP_convert_D_to_raw QOP_F3_convert_D_to_raw
#define QOP_convert_G_to_raw QOP_F3_convert_G_to_raw
#define QOP_convert_F_to_raw QOP_F3_convert_F_to_raw

#define QOP_rephase_G QOP_F3_rephase_G

#define QOP_symanzik_1loop_gauge_action  QOP_F3_symanzik_1loop_gauge_action
#define QOP_symanzik_1loop_gauge_force   QOP_F3_symanzik_1loop_gauge_force
#define QOP_symanzik_1loop_gauge_deriv   QOP_F3_symanzik_1loop_gauge_deriv

#define QOP_FermionLinksAsqtad        QOP_F3_FermionLinksAsqtad
#define QOP_asqtad_create_L_from_raw  QOP_F3_asqtad_create_L_from_raw
#define QOP_asqtad_create_L_from_G    QOP_F3_asqtad_create_L_from_G
#define QOP_asqtad_extract_L_to_raw   QOP_F3_asqtad_extract_L_to_raw
#define QOP_asqtad_destroy_L          QOP_F3_asqtad_destroy_L
#define QOP_asqtad_convert_L_from_raw QOP_F3_asqtad_convert_L_from_raw
#define QOP_asqtad_convert_L_to_raw   QOP_F3_asqtad_convert_L_to_raw
#define QOP_asqtad_load_L_from_raw    QOP_F3_asqtad_load_L_from_raw
#define QOP_asqtad_load_L_from_G      QOP_F3_asqtad_load_L_from_G
#define QOP_asqtad_rephase_L          QOP_F3_asqtad_rephase_L

#define QOP_asqtad_dslash       QOP_F3_asqtad_dslash
#define QOP_asqtad_dslash_dir   QOP_F3_asqtad_dslash_dir
#define QOP_asqtad_diaginv      QOP_F3_asqtad_diaginv
#define QOP_asqtad_invert       QOP_F3_asqtad_invert
#define QOP_asqtad_invert_multi QOP_F3_asqtad_invert_multi
#define QOP_asqtad_force        QOP_F3_asqtad_force
#define QOP_asqtad_force_multi  QOP_F3_asqtad_force_multi


#define QOP_FermionLinksHisq        QOP_F3_FermionLinksHisq
#define QOP_hisq_create_L_from_G    QOP_F3_hisq_create_L_from_G
#define QOP_hisq_destroy_L          QOP_F3_hisq_destroy_L
#define QOP_get_asqtad_links_from_hisq QOP_F3_get_asqtad_links_from_hisq
#define QOP_get_asqtad_deps_links_from_hisq QOP_F3_get_asqtad_deps_links_from_hisq

#define QOP_hisq_invert       QOP_F3_hisq_invert
#define QOP_hisq_invert_multi QOP_F3_hisq_invert_multi

#define QOP_hisq_force        QOP_F3_hisq_force
#define QOP_hisq_force_multi  QOP_F3_hisq_force_multi

#define QOP_FermionLinksWilson        QOP_F3_FermionLinksWilson
#define QOP_wilson_create_L_from_raw  QOP_F3_wilson_create_L_from_raw
#define QOP_wilson_create_L_from_G    QOP_F3_wilson_create_L_from_G
#define QOP_wilson_extract_L_to_raw   QOP_F3_wilson_extract_L_to_raw
#define QOP_wilson_destroy_L          QOP_F3_wilson_destroy_L
#define QOP_wilson_convert_L_from_raw QOP_F3_wilson_convert_L_from_raw
#define QOP_wilson_convert_L_to_raw   QOP_F3_wilson_convert_L_to_raw
#define QOP_wilson_convert_L_from_G   QOP_F3_wilson_convert_L_from_G
#define QOP_wilson_convert_L_to_G     QOP_F3_wilson_convert_L_to_G
#define QOP_wilson_load_L_from_raw    QOP_F3_wilson_load_L_from_raw
#define QOP_wilson_load_L_from_G      QOP_F3_wilson_load_L_from_G

#define QOP_wilson_dslash       QOP_F3_wilson_dslash
#define QOP_wilson_invert       QOP_F3_wilson_invert
#define QOP_wilson_invert_multi QOP_F3_wilson_invert_multi
#define QOP_wilson_force        QOP_F3_wilson_force
#define QOP_wilson_force_multi  QOP_F3_wilson_force_multi

#define QOP_FermionLinksDW        QOP_F3_FermionLinksDW
#define QOP_dw_create_L_from_raw  QOP_F3_dw_create_L_from_raw
#define QOP_dw_create_L_from_G    QOP_F3_dw_create_L_from_G
#define QOP_dw_extract_L_to_raw   QOP_F3_dw_extract_L_to_raw
#define QOP_dw_destroy_L          QOP_F3_dw_destroy_L
#define QOP_dw_convert_L_from_raw QOP_F3_dw_convert_L_from_raw
#define QOP_dw_convert_L_to_raw   QOP_F3_dw_convert_L_to_raw
#define QOP_dw_convert_L_from_G   QOP_F3_dw_convert_L_from_G
#define QOP_dw_convert_L_to_G     QOP_F3_dw_convert_L_to_G
#define QOP_dw_load_L_from_raw    QOP_F3_dw_load_L_from_raw
#define QOP_dw_load_L_from_G      QOP_F3_dw_load_L_from_G

#define QOP_dw_dslash       QOP_F3_dw_dslash
#define QOP_dw_invert       QOP_F3_dw_invert
#define QOP_dw_invert_multi QOP_F3_dw_invert_multi
#define QOP_dw_force        QOP_F3_dw_force
#define QOP_dw_force_multi  QOP_F3_dw_force_multi

#else

#define QOP_ColorVector  QOP_D3_ColorVector
#define QOP_DiracFermion QOP_D3_DiracFermion
#define QOP_GaugeField   QOP_D3_GaugeField
#define QOP_Force        QOP_D3_Force

#define QOP_create_V_from_raw QOP_D3_create_V_from_raw
#define QOP_create_D_from_raw QOP_D3_create_D_from_raw
#define QOP_create_G_from_raw QOP_D3_create_G_from_raw
#define QOP_create_F_from_raw QOP_D3_create_F_from_raw

#define QOP_extract_V_to_raw QOP_D3_extract_V_to_raw
#define QOP_extract_D_to_raw QOP_D3_extract_D_to_raw
#define QOP_extract_G_to_raw QOP_D3_extract_G_to_raw
#define QOP_extract_F_to_raw QOP_D3_extract_F_to_raw

#define QOP_destroy_V QOP_D3_destroy_V
#define QOP_destroy_D QOP_D3_destroy_D
#define QOP_destroy_G QOP_D3_destroy_G
#define QOP_destroy_F QOP_D3_destroy_F

#define QOP_convert_V_from_raw QOP_D3_convert_V_from_raw
#define QOP_convert_D_from_raw QOP_D3_convert_D_from_raw
#define QOP_convert_G_from_raw QOP_D3_convert_G_from_raw
#define QOP_convert_F_from_raw QOP_D3_convert_F_from_raw

#define QOP_convert_V_to_raw QOP_D3_convert_V_to_raw
#define QOP_convert_D_to_raw QOP_D3_convert_D_to_raw
#define QOP_convert_G_to_raw QOP_D3_convert_G_to_raw
#define QOP_convert_F_to_raw QOP_D3_convert_F_to_raw

#define QOP_rephase_G QOP_D3_rephase_G

#define QOP_symanzik_1loop_gauge_action  QOP_D3_symanzik_1loop_gauge_action
#define QOP_symanzik_1loop_gauge_force   QOP_D3_symanzik_1loop_gauge_force
#define QOP_symanzik_1loop_gauge_deriv   QOP_D3_symanzik_1loop_gauge_deriv

#define QOP_FermionLinksAsqtad       QOP_D3_FermionLinksAsqtad
#define QOP_asqtad_create_L_from_raw QOP_D3_asqtad_create_L_from_raw
#define QOP_asqtad_create_L_from_G   QOP_D3_asqtad_create_L_from_G
#define QOP_asqtad_extract_L_to_raw  QOP_D3_asqtad_extract_L_to_raw
#define QOP_asqtad_destroy_L         QOP_D3_asqtad_destroy_L
#define QOP_asqtad_convert_L_from_raw QOP_D3_asqtad_convert_L_from_raw
#define QOP_asqtad_convert_L_to_raw   QOP_D3_asqtad_convert_L_to_raw
#define QOP_asqtad_load_L_from_raw    QOP_D3_asqtad_load_L_from_raw
#define QOP_asqtad_load_L_from_G      QOP_D3_asqtad_load_L_from_G
#define QOP_asqtad_rephase_L          QOP_D3_asqtad_rephase_L

#define QOP_asqtad_dslash       QOP_D3_asqtad_dslash
#define QOP_asqtad_dslash_dir   QOP_D3_asqtad_dslash_dir
#define QOP_asqtad_diaginv      QOP_D3_asqtad_diaginv
#define QOP_asqtad_invert       QOP_D3_asqtad_invert
#define QOP_asqtad_invert_multi QOP_D3_asqtad_invert_multi
#define QOP_asqtad_force        QOP_D3_asqtad_force
#define QOP_asqtad_force_multi  QOP_D3_asqtad_force_multi

#define QOP_FermionLinksHisq       QOP_D3_FermionLinksHisq
#define QOP_hisq_create_L_from_G   QOP_D3_hisq_create_L_from_G
#define QOP_hisq_destroy_L         QOP_D3_hisq_destroy_L
#define QOP_get_asqtad_links_from_hisq QOP_D3_get_asqtad_links_from_hisq
#define QOP_get_asqtad_deps_links_from_hisq QOP_D3_get_asqtad_deps_links_from_hisq

#define QOP_hisq_invert       QOP_D3_hisq_invert
#define QOP_hisq_invert_multi QOP_D3_hisq_invert_multi

#define QOP_hisq_force        QOP_D3_hisq_force
#define QOP_hisq_force_multi  QOP_D3_hisq_force_multi

#define QOP_FermionLinksWilson        QOP_D3_FermionLinksWilson
#define QOP_wilson_create_L_from_raw  QOP_D3_wilson_create_L_from_raw
#define QOP_wilson_create_L_from_G    QOP_D3_wilson_create_L_from_G
#define QOP_wilson_extract_L_to_raw   QOP_D3_wilson_extract_L_to_raw
#define QOP_wilson_destroy_L          QOP_D3_wilson_destroy_L
#define QOP_wilson_convert_L_from_raw QOP_D3_wilson_convert_L_from_raw
#define QOP_wilson_convert_L_to_raw   QOP_D3_wilson_convert_L_to_raw
#define QOP_wilson_convert_L_from_G   QOP_D3_wilson_convert_L_from_G
#define QOP_wilson_convert_L_to_G     QOP_D3_wilson_convert_L_to_G
#define QOP_wilson_load_L_from_raw    QOP_D3_wilson_load_L_from_raw
#define QOP_wilson_load_L_from_G      QOP_D3_wilson_load_L_from_G

#define QOP_wilson_dslash       QOP_D3_wilson_dslash
#define QOP_wilson_diaginv      QOP_D3_wilson_diaginv
#define QOP_wilson_invert       QOP_D3_wilson_invert
#define QOP_wilson_invert_multi QOP_D3_wilson_invert_multi
#define QOP_wilson_force        QOP_D3_wilson_force
#define QOP_wilson_force_multi  QOP_D3_wilson_force_multi

#define QOP_FermionLinksDW        QOP_D3_FermionLinksDW
#define QOP_dw_create_L_from_raw  QOP_D3_dw_create_L_from_raw
#define QOP_dw_create_L_from_G    QOP_D3_dw_create_L_from_G
#define QOP_dw_extract_L_to_raw   QOP_D3_dw_extract_L_to_raw
#define QOP_dw_destroy_L          QOP_D3_dw_destroy_L
#define QOP_dw_convert_L_from_raw QOP_D3_dw_convert_L_from_raw
#define QOP_dw_convert_L_to_raw   QOP_D3_dw_convert_L_to_raw
#define QOP_dw_convert_L_from_G   QOP_D3_dw_convert_L_from_G
#define QOP_dw_convert_L_to_G     QOP_D3_dw_convert_L_to_G
#define QOP_dw_load_L_from_raw    QOP_D3_dw_load_L_from_raw
#define QOP_dw_load_L_from_G      QOP_D3_dw_load_L_from_G

#define QOP_dw_dslash       QOP_D3_dw_dslash
#define QOP_dw_invert       QOP_D3_dw_invert
#define QOP_dw_invert_multi QOP_D3_dw_invert_multi
#define QOP_dw_force        QOP_D3_dw_force
#define QOP_dw_force_multi  QOP_D3_dw_force_multi

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_H */
