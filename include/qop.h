#ifndef _QOP_H
#define _QOP_H

#ifdef __cplusplus
extern "C" {
#endif

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
} QOP_info_t;

  /* these are quantities that apply to all masses in the multi inverter */
typedef struct {
  int max_iter;            /* (in) max number of iterations before restart */
  int restart;             /* (in) number of restarts allowed */
  QOP_evenodd_t evenodd;   /* (in) subset of source vector */
} QOP_invert_arg_t;

  /* these are quantities that vary for each mass in the multi inverter */
typedef struct {
  double rsqmin;           /* (in) desired squared residual */
  double final_rsq;        /* (out) actual squared residual */
  int final_iter;          /* (out) number of iterations done */
} QOP_resid_arg_t;

typedef struct QOP_F3_ColorVector_struct  QOP_F3_ColorVector;
typedef struct QOP_F3_DiracFermion_struct QOP_F3_DiracFermion;
typedef struct QOP_F3_GaugeField_struct   QOP_F3_GaugeField;
typedef struct QOP_F3_Force_struct        QOP_F3_Force;

typedef struct QOP_D3_ColorVector_struct  QOP_D3_ColorVector;
typedef struct QOP_D3_DiracFermion_struct QOP_D3_DiracFermion;
typedef struct QOP_D3_GaugeField_struct   QOP_D3_GaugeField;
typedef struct QOP_D3_Force_struct        QOP_D3_Force;

   /* Imp gauge */
typedef struct{  
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

  /* Wilson datatypes */

typedef struct {
  double clov_c;
} QOP_wilson_coeffs_t;

typedef struct QOP_F3_FermionLinksWilson_struct QOP_F3_FermionLinksWilson;
typedef struct QOP_D3_FermionLinksWilson_struct QOP_D3_FermionLinksWilson;

  /* Domain Wall datatypes */

typedef struct {
  double clov_c;
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
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);

void QOP_D3_rephase_G(QOP_D3_GaugeField *links,
		      QOP_bc_t *bc,
		      QOP_staggered_sign_t *sign);


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

  /* inverter routines */

QOP_status_t QOP_asqtad_invert_set_opts(QOP_opt_t opts[], int nopts);

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
/*  gauge force routines */

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


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* fermion matrix link routines */

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_raw(float *links[], float *clov[],
				  QOP_evenodd_t evenodd);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_F3_GaugeField *gauge);

void QOP_F3_wilson_extract_L_to_raw(float *links[], float *clov[],
				    QOP_F3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_F3_wilson_destroy_L(QOP_F3_FermionLinksWilson *field);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_raw(float *links[], float *clov[],
				   QOP_evenodd_t evenodd);

void QOP_F3_wilson_convert_L_to_raw(float ***links, float ***clov,
				    QOP_F3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_F3_GaugeField *gauge);

QOP_F3_GaugeField *
  QOP_F3_wilson_convert_L_to_G(QOP_F3_FermionLinksWilson *links);

  /* double precision */

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_raw(double *links[], double *clov[],
				  QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_G(QOP_info_t *info,
				QOP_wilson_coeffs_t *coeffs,
				QOP_D3_GaugeField *gauge);

void QOP_D3_wilson_extract_L_to_raw(double *links[], double *clov[],
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

void QOP_D3_wilson_destroy_L(QOP_D3_FermionLinksWilson *field);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_raw(double *links[], double *clov[],
				   QOP_evenodd_t evenodd);

void QOP_D3_wilson_convert_L_to_raw(double ***links, double ***clov,
				    QOP_D3_FermionLinksWilson *src,
				    QOP_evenodd_t evenodd);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_G(QOP_info_t *info,
				 QOP_wilson_coeffs_t *coeffs,
				 QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_wilson_convert_L_to_G(QOP_D3_FermionLinksWilson *links);

  /* inverter routines */

QOP_status_t QOP_wilson_invert_set_opts(QOP_opt_t opts[], int nopts);

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
  QOP_F3_dw_create_L_from_raw(float *links[], float *clov[],
			      QOP_evenodd_t evenodd);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_F3_GaugeField *gauge);

void QOP_F3_dw_extract_L_to_raw(float *links[], float *clov[],
				QOP_F3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_F3_dw_destroy_L(QOP_F3_FermionLinksDW *field);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_convert_L_from_raw(float *links[], float *clov[],
			       QOP_evenodd_t evenodd);

void QOP_F3_dw_convert_L_to_raw(float ***links, float ***clov,
				QOP_F3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_F3_GaugeField *gauge);

QOP_F3_GaugeField *
  QOP_F3_dw_convert_L_to_G(QOP_F3_FermionLinksDW *links);

  /* double precision */

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_raw(double *links[], double *clov[],
			      QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_G(QOP_info_t *info,
			    QOP_dw_coeffs_t *coeffs,
			    QOP_D3_GaugeField *gauge);

void QOP_D3_dw_extract_L_to_raw(double *links[], double *clov[],
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

void QOP_D3_dw_destroy_L(QOP_D3_FermionLinksDW *field);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_raw(double *links[], double *clov[],
			       QOP_evenodd_t evenodd);

void QOP_D3_dw_convert_L_to_raw(double ***links, double ***clov,
				QOP_D3_FermionLinksDW *src,
				QOP_evenodd_t evenodd);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_G(QOP_info_t *info,
			     QOP_dw_coeffs_t *coeffs,
			     QOP_D3_GaugeField *gauge);

QOP_D3_GaugeField *
  QOP_D3_dw_convert_L_to_G(QOP_D3_FermionLinksDW *links);

  /* inverter routines */

QOP_status_t QOP_dw_invert_set_opts(QOP_opt_t opts[], int nopts);

void QOP_F3_dw_invert(QOP_info_t *info,
		      QOP_F3_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      float m0,
		      float M,
		      QOP_F3_DiracFermion *out_pt[],
		      QOP_F3_DiracFermion *in_pt[],
		      int Ls);

void QOP_F3_dw_invert_multi(QOP_info_t *info,
			    QOP_F3_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    float *m0[],
			    float *M[],
			    int nmass[],
			    QOP_F3_DiracFermion ***out_pt[],
			    QOP_F3_DiracFermion **in_pt[],
			    int Ls,
			    int nsrc);

void QOP_D3_dw_invert(QOP_info_t *info,
		      QOP_D3_FermionLinksDW *links,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      double m0,
		      double M,
		      QOP_D3_DiracFermion *out_pt[],
		      QOP_D3_DiracFermion *in_pt[],
		      int Ls);

void QOP_D3_dw_invert_multi(QOP_info_t *info,
			    QOP_D3_FermionLinksDW *links,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    double *m0[],
			    double *M[],
			    int nmass[],
			    QOP_D3_DiracFermion ***out_pt[],
			    QOP_D3_DiracFermion **in_pt[],
			    int Ls,
			    int nsrc);

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

#define QOP_FermionLinksAsqtad        QOP_F3_FermionLinksAsqtad
#define QOP_asqtad_create_L_from_raw  QOP_F3_asqtad_create_L_from_raw
#define QOP_asqtad_create_L_from_G    QOP_F3_asqtad_create_L_from_G
#define QOP_asqtad_extract_L_to_raw   QOP_F3_asqtad_extract_L_to_raw
#define QOP_asqtad_destroy_L          QOP_F3_asqtad_destroy_L
#define QOP_asqtad_convert_L_from_raw QOP_F3_asqtad_convert_L_from_raw
#define QOP_asqtad_convert_L_to_raw   QOP_F3_asqtad_convert_L_to_raw

#define QOP_asqtad_invert       QOP_F3_asqtad_invert
#define QOP_asqtad_invert_multi QOP_F3_asqtad_invert_multi
#define QOP_asqtad_force        QOP_F3_asqtad_force
#define QOP_asqtad_force_multi  QOP_F3_asqtad_force_multi
#define QOP_symanzik_1loop_gauge_force     QOP_F3_symanzik_1loop_gauge_force

#define QOP_FermionLinksWilson        QOP_F3_FermionLinksWilson
#define QOP_wilson_create_L_from_raw  QOP_F3_wilson_create_L_from_raw
#define QOP_wilson_create_L_from_G    QOP_F3_wilson_create_L_from_G
#define QOP_wilson_extract_L_to_raw   QOP_F3_wilson_extract_L_to_raw
#define QOP_wilson_destroy_L          QOP_F3_wilson_destroy_L
#define QOP_wilson_convert_L_from_raw QOP_F3_wilson_convert_L_from_raw
#define QOP_wilson_convert_L_to_raw   QOP_F3_wilson_convert_L_to_raw
#define QOP_wilson_convert_L_from_G   QOP_F3_wilson_convert_L_from_G
#define QOP_wilson_convert_L_to_G     QOP_F3_wilson_convert_L_to_G

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

#define QOP_FermionLinksAsqtad       QOP_D3_FermionLinksAsqtad
#define QOP_asqtad_create_L_from_raw QOP_D3_asqtad_create_L_from_raw
#define QOP_asqtad_create_L_from_G   QOP_D3_asqtad_create_L_from_G
#define QOP_asqtad_extract_L_to_raw  QOP_D3_asqtad_extract_L_to_raw
#define QOP_asqtad_destroy_L         QOP_D3_asqtad_destroy_L
#define QOP_asqtad_convert_L_from_raw QOP_D3_asqtad_convert_L_from_raw
#define QOP_asqtad_convert_L_to_raw   QOP_D3_asqtad_convert_L_to_raw

#define QOP_asqtad_invert       QOP_D3_asqtad_invert
#define QOP_asqtad_invert_multi QOP_D3_asqtad_invert_multi
#define QOP_asqtad_force        QOP_D3_asqtad_force
#define QOP_asqtad_force_multi  QOP_D3_asqtad_force_multi
#define QOP_symanzik_1loop_gauge_force     QOP_D3_symanzik_1loop_gauge_force

#define QOP_FermionLinksWilson        QOP_D3_FermionLinksWilson
#define QOP_wilson_create_L_from_raw  QOP_D3_wilson_create_L_from_raw
#define QOP_wilson_create_L_from_G    QOP_D3_wilson_create_L_from_G
#define QOP_wilson_extract_L_to_raw   QOP_D3_wilson_extract_L_to_raw
#define QOP_wilson_destroy_L          QOP_D3_wilson_destroy_L
#define QOP_wilson_convert_L_from_raw QOP_D3_wilson_convert_L_from_raw
#define QOP_wilson_convert_L_to_raw   QOP_D3_wilson_convert_L_to_raw
#define QOP_wilson_convert_L_from_G   QOP_D3_wilson_convert_L_from_G
#define QOP_wilson_convert_L_to_G     QOP_D3_wilson_convert_L_to_G

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

#define QOP_dw_invert       QOP_D3_dw_invert
#define QOP_dw_invert_multi QOP_D3_dw_invert_multi
#define QOP_dw_force        QOP_D3_dw_force
#define QOP_dw_force_multi  QOP_D3_dw_force_multi

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_H */
