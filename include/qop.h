#ifndef _QOPQDP_H
#define _QOPQDP_H

#include <qdp.h>

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
QOP_status_t QOP_asqtad_set_opt(char *tag, double value);

QOP_status_t QOP_F_asqtad_invert_load_links_raw(float *fatlinks[],
                                                float *longlinks[]);
QOP_status_t QOP_F_asqtad_invert_load_links_qdp(QDP_F3_ColorMatrix *fatlnk[],
                                                QDP_F3_ColorMatrix *longlnk[]);
int QOP_F_asqtad_inv_raw(QOP_invert_arg *inv_arg,
                         float *out_pt, float *in_pt);
int QOP_F_asqtad_inv_qdp(QOP_invert_arg *inv_arg,
                         QDP_F3_ColorVector *out, QDP_F3_ColorVector *in);
QOP_status_t QOP_F_asqtad_invert_unload_links(void);

QOP_status_t QOP_D_asqtad_invert_load_links_raw(double *fatlinks[],
                                                double *longlinks[]);
QOP_status_t QOP_D_asqtad_invert_load_links_qdp(QDP_D3_ColorMatrix *fatlnk[],
                                                QDP_D3_ColorMatrix *longlnk[]);
int QOP_D_asqtad_inv_raw(QOP_invert_arg *inv_arg,
                         double *out_pt, double *in_pt);
int QOP_D_asqtad_inv_qdp(QOP_invert_arg *inv_arg,
                         QDP_D3_ColorVector *out, QDP_D3_ColorVector *in);
QOP_status_t QOP_D_asqtad_invert_unload_links(void);


  /*  Wilson routines  */

QOP_status_t QOP_wilson_invert_init(QOP_layout *layout);
QOP_status_t QOP_wilson_invert_finalize(void);
QOP_status_t QOP_wilson_set_opt(char *tag, double value);

QOP_status_t QOP_F_wilson_invert_load_links_raw(float *links[]);
QOP_status_t QOP_F_wilson_invert_load_links_qdp(QDP_F3_ColorMatrix *links[]);
int QOP_F_wilson_inv_raw(QOP_invert_arg *inv_arg, float *out_pt, float *in_pt);
int QOP_F_wilson_inv_qdp(QOP_invert_arg *inv_arg,
                         QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in);
QOP_status_t QOP_F_wilson_invert_unload_links(void);

QOP_status_t QOP_D_wilson_invert_load_links_raw(double *links[]);
QOP_status_t QOP_D_wilson_invert_load_links_qdp(QDP_D3_ColorMatrix *links[]);
int QOP_D_wilson_inv_raw(QOP_invert_arg *inv_arg,
			 double *out_pt, double *in_pt);
int QOP_D_wilson_inv_qdp(QOP_invert_arg *inv_arg,
                         QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in);
QOP_status_t QOP_D_wilson_invert_unload_links(void);


#if QOP_Precision == 1

#define QOP_asqtad_invert_load_links_raw QOP_F_asqtad_invert_load_links_raw
#define QOP_asqtad_invert_load_links_qdp QOP_F_asqtad_invert_load_links_qdp
#define QOP_asqtad_inv_raw QOP_F_asqtad_inv_raw
#define QOP_asqtad_inv_qdp QOP_F_asqtad_inv_qdp
#define QOP_asqtad_invert_unload_links QOP_F_asqtad_invert_unload_links

#define QOP_wilson_invert_load_links_raw QOP_F_wilson_invert_load_links_raw
#define QOP_wilson_invert_load_links_qdp QOP_F_wilson_invert_load_links_qdp
#define QOP_wilson_inv_raw QOP_F_wilson_inv_raw
#define QOP_wilson_inv_qdp QOP_F_wilson_inv_qdp
#define QOP_wilson_invert_unload_links QOP_F_wilson_invert_unload_links

#else

#define QOP_asqtad_invert_load_links_raw QOP_D_asqtad_invert_load_links_raw
#define QOP_asqtad_invert_load_links_qdp QOP_D_asqtad_invert_load_links_qdp
#define QOP_asqtad_inv_raw QOP_D_asqtad_inv_raw
#define QOP_asqtad_inv_qdp QOP_D_asqtad_inv_qdp
#define QOP_asqtad_invert_unload_links QOP_D_asqtad_invert_unload_links

#define QOP_wilson_invert_load_links_raw QOP_D_wilson_invert_load_links_raw
#define QOP_wilson_invert_load_links_qdp QOP_D_wilson_invert_load_links_qdp
#define QOP_wilson_inv_raw QOP_D_wilson_inv_raw
#define QOP_wilson_inv_qdp QOP_D_wilson_inv_qdp
#define QOP_wilson_invert_unload_links QOP_D_wilson_invert_unload_links

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOPQDP_H */
