#ifndef _QOP_INT_H
#define _QOP_INT_H

#ifdef __cplusplus
extern "C" {
#endif

  /* Enumerate in order of increasing verbosity */
#define QOP_VERB_OFF    0
#define QOP_VERB_LOW    1
#define QOP_VERB_MED    2
#define QOP_VERB_HI     3
#define QOP_VERB_DEBUG  4

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

typedef float QOP_F_Real;
typedef double QOP_D_Real;

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
#define QOP_LAYOUT_ZERO ((QOP_layout_t){NULL,NULL,0,NULL,0,NULL,0,0})

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
#define QOP_INFO_ZERO ((QOP_info_t){0,0,QOP_SUCCESS,0,0})
#define QOP_info_zero(a) do { \
 (a).final_sec = 0; \
 (a).final_flop = 0; \
 (a).status = QOP_SUCCESS; \
 (a).count1 = 0; \
 (a).count2 = 0; \
} while(0)
#define QOP_info_copy(a,b) do { \
 (a).final_sec = (b).final_sec; \
 (a).final_flop = (b).final_flop; \
 (a).status = (b).status; \
 (a).count1 = (b).count1; \
 (a).count2 = (b).count2; \
} while(0)
#define QOP_info_addfrom(a,b) do { \
 (a).final_sec += (b).final_sec; \
 (a).final_flop += (b).final_flop; \
 if((b).status>(a).status) (a).status = (b).status; \
 (a).count1 += (b).count1; \
 (a).count2 += (b).count2; \
} while(0)
  /* QOP_info_t counters as they are used in HISQ routines */
  /* The user must initialize them and reset them */
#define QOP_info_hisq_svd_counter(info) ((info)->count1)
#define QOP_info_hisq_force_filter_counter(info) ((info)->count2) 

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  double mixed_rsq;        // rsq above which to use mixed precision
  int max_iter;            /* (in) max number of iterations */
  int restart;             /* (in) number of iterations before restart */
  int max_restarts;        /* (in) number of restarts allowed */
  QOP_evenodd_t evenodd;   /* (in) subset of source vector */
} QOP_invert_arg_t;
#define QOP_INVERT_ARG_DEFAULT ((QOP_invert_arg_t){1e9,2000,1000,5,QOP_EVENODD})

  /* these are quantities that vary for each mass in the multi inverter */
typedef struct {
  double rsqmin;         /* (in) desired squared residual. Ignored if 0 */
  double final_rsq;      /* (out) actual squared residual */
  double relmin;         /* (in) desired squared relative norm. Ignored if 0 */
  double final_rel;      /* (out) actual squared relative norm */
  int final_iter;        /* (out) number of iterations done */
  int final_restart;     /* (out) number of restarts done */
} QOP_resid_arg_t;
#define QOP_RESID_ARG_DEFAULT ((QOP_resid_arg_t){1e-6,0,0,0,0,0})

  /* Improved gauge datatypes */

typedef struct {
   double plaquette;
   double rectangle;
   double parallelogram;
   double adjoint_plaquette;
} QOP_gauge_coeffs_t;
#define QOP_GAUGE_COEFFS_ZERO ((QOP_gauge_coeffs_t){0,0,0,0})

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
#define QOP_ASQTAD_COEFFS_ZERO ((QOP_asqtad_coeffs_t){0,0,0,0,0,0})

/* Maximum number of Naik terms */
#define QOP_MAX_NAIK 20
#define QOP_EPS_NAIK_ZERO {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}

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
#define QOP_HISQ_COEFFS_ZERO \
  ((QOP_hisq_coeffs_t){1, QOP_EPS_NAIK_ZERO, QOP_UNITARIZE_U3, \
      QOP_UNITARIZE_RATIONAL, 0,0,0,0,0, 0,0,0,0,0,0, 0,0})

  /* Wilson datatypes */
typedef struct {
  double clov_s;
  double clov_t;
  double aniso;
} QOP_wilson_coeffs_t;
#define QOP_WILSON_COEFFS_ZERO ((QOP_wilson_coeffs_t){0,0,1})

  /* Oktay-Kronfeld OK action datatypes */
typedef struct {
  double kapifla ;
  double kappa_s ;
  double kappa_t ;
  double r_s     ;
  double r_t     ;
  double zeta    ;
  double c_E     ; 
  double c_B     ; 
  double c_1     ;
  double c_2     ;
  double c_3     ;
  double c_4     ;
  double c_5     ;
  double c_EE    ; 
  double u0      ;
} QOP_wilson_ifla_coeffs_t;
#define QOP_WILSON_IFLA_COEFFS_ZERO \
  ((QOP_wilson_ifla_coeffs_t){0,0,0, 0,0, 0, 0,0,0,0,0,0,0,0, 0})

  /* Domain Wall datatypes */
typedef struct {
  // need to add Mobius parameters here
  double dummy;
} QOP_dw_coeffs_t;
#define QOP_DW_COEFFS_ZERO ((QOP_dw_coeffs_t){0})


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


QOP_status_t QOP_asqtad_invert_set_opts(QOP_opt_t opts[], int nopts);
QOP_status_t QOP_asqtad_force_set_opts(QOP_opt_t opts[], int nopts);

QOP_status_t QOP_hisq_links_set_opts(QOP_opt_t opts[], int nopts);
QOP_status_t QOP_hisq_force_set_opts(QOP_opt_t opts[], int nopts);

QOP_status_t QOP_wilson_invert_set_opts(QOP_opt_t opts[], int nopts);
QOP_status_t QOP_wilson_force_set_opts(QOP_opt_t opts[], int nopts);

QOP_status_t QOP_dw_invert_set_opts(QOP_opt_t opts[], int nopts);
QOP_status_t QOP_dw_force_set_opts(QOP_opt_t opts[], int nopts);

#ifdef __cplusplus
}
#endif

#endif /* _QOP_INT_H */
