#ifndef _QOP_H
#define _QOP_H

#ifdef __cplusplus
extern "C" {
#endif

#define QOP_MAX_DIM 6

typedef enum QOP_status {
  QOP_SUCCESS = 0,
  QOP_FAIL = 1,
  QOP_BC_ERROR,
  QOP_MEM_ERROR,
  QOP_TOP_ERROR
} QOP_status_t;

typedef enum QOP_evenodd {
  QOP_EVEN = 0,
  QOP_ODD = 1,
  QOP_EVENODD = 2
} QOP_evenodd_t;

typedef struct QOP_layout_i {
  int ndims;
  int sites[QOP_MAX_DIM];
  int bc[QOP_MAX_DIM];
} QOP_layout;

typedef struct QOP_invert_arg_i {
  double mass;
  double rsqmin;
  double final_rsq;
  double final_sec;
  double final_flop;
  int max_iter;
  int restart;
  int final_iter;
  QOP_evenodd_t evenodd;
} QOP_invert_arg;


  /*  Asqtad routines  */

QOP_status_t QOP_asqtad_invert_init(QOP_layout *layout);
QOP_status_t QOP_asqtad_invert_finalize(void);
QOP_status_t QOP_asqtad_invert_set_opt(char *tag, double value);

QOP_status_t QOP_F_asqtad_invert_load_links_raw(float *fatlinks[],
                                                float *longlinks[]);
int QOP_F_asqtad_invert_raw(QOP_invert_arg *inv_arg,
			    float *out_pt, float *in_pt);
QOP_status_t QOP_F_asqtad_invert_unload_links(void);

QOP_status_t QOP_D_asqtad_invert_load_links_raw(double *fatlinks[],
                                                double *longlinks[]);
int QOP_D_asqtad_invert_raw(QOP_invert_arg *inv_arg,
			    double *out_pt, double *in_pt);
QOP_status_t QOP_D_asqtad_invert_unload_links(void);


  /*  Wilson routines  */

QOP_status_t QOP_wilson_invert_init(QOP_layout *layout);
QOP_status_t QOP_wilson_invert_finalize(void);
QOP_status_t QOP_wilson_invert_set_opt(char *tag, double value);

QOP_status_t QOP_F_wilson_invert_load_links_raw(float *links[]);
int QOP_F_wilson_invert_raw(QOP_invert_arg *inv_arg,
			    float *out_pt, float *in_pt);
QOP_status_t QOP_F_wilson_invert_unload_links(void);

QOP_status_t QOP_D_wilson_invert_load_links_raw(double *links[]);
int QOP_D_wilson_invert_raw(QOP_invert_arg *inv_arg,
			    double *out_pt, double *in_pt);
QOP_status_t QOP_D_wilson_invert_unload_links(void);


  /* Mapping of generic names to specific precision */

#ifndef QOP_Precision
#define QOP_Precision 1
#endif

#if QOP_Precision == 1

#define QOP_asqtad_invert_load_links_raw QOP_F_asqtad_invert_load_links_raw
#define QOP_asqtad_invert_raw QOP_F_asqtad_invert_raw
#define QOP_asqtad_invert_unload_links QOP_F_asqtad_invert_unload_links

#define QOP_wilson_invert_load_links_raw QOP_F_wilson_invert_load_links_raw
#define QOP_wilson_invert_raw QOP_F_wilson_invert_raw
#define QOP_wilson_invert_unload_links QOP_F_wilson_invert_unload_links

#else

#define QOP_asqtad_invert_load_links_raw QOP_D_asqtad_invert_load_links_raw
#define QOP_asqtad_invert_raw QOP_D_asqtad_invert_raw
#define QOP_asqtad_invert_unload_links QOP_D_asqtad_invert_unload_links

#define QOP_wilson_invert_load_links_raw QOP_D_wilson_invert_load_links_raw
#define QOP_wilson_invert_raw QOP_D_wilson_invert_raw
#define QOP_wilson_invert_unload_links QOP_D_wilson_invert_unload_links

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_H */
