#ifndef _QOP_INTERNAL_H
#define _QOP_INTERNAL_H

#include <qop.h>
#include <qop_qdp.h>
#include <qop_internal_p.h>
#include <qmp.h>

#define oppsub(eo) ((4-(eo))%3)
#define qdpsub(eo) ((eo)==2 ? QDP_all : QDP_even_and_odd[eo])

#define CHECK_INIT							\
  if(!QOP_common.inited) {						\
    QMP_error("Error: QOP not initialized in function %s\n", __func__); \
    exit(1);								\
  }
#define CHECK_NOT_INIT						       \
  if(QOP_common.inited) {					       \
    QMP_error("Error: QOP is initialized in function %s\n", __func__); \
    exit(1);							       \
  }

#define QOP_malloc(var, type, num)					\
  (var) = (type *) malloc(num*sizeof(type));				\
  if(!(var)) {								\
    QMP_error("Error: QOP ran out of memory in function %s\n", __func__); \
    exit(1);								\
  }

#define QOP_printf0 if(QDP_this_node==0) printf

#define QOP_error(str) \
  fprintf(stderr, "QOP error: %s\n", str); \
  QDP_abort();

typedef struct {
  int inited;
  int verbosity;
  int proflevel;
  int we_inited_qdp;
  int ndim;
  QDP_Shift neighbor3[4];
  QDP_ShiftDir shiftfwd[8], shiftbck[8];
} QOP_common_t;
extern QOP_common_t QOP_common;

typedef struct {
  int inited;
  int style;
  int nsvec;
  int nvec;
  QDP_Shift shifts[8];
  QDP_Shift shifts_dbl[16];
  QDP_ShiftDir shiftdirs_dbl[16];
} QOP_asqtad_t;
extern QOP_asqtad_t QOP_asqtad;

typedef struct {
  int inited;
  int fnmat_src_min;
  int veclength;
} QOP_asqtad_force_t;
extern QOP_asqtad_force_t QOP_asqtad_ff;

#ifdef __cplusplus
extern "C" {
#endif

double QOP_time(void);
QOP_status_t QOP_asqtad_invert_init(void);

#ifdef __cplusplus
}
#endif

#endif /* _QOP_INTERNAL_H */
