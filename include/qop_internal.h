#ifndef _QOP_INTERNAL_H
#define _QOP_INTERNAL_H

#include <qop.h>
#include <qop_qdp.h>
#include <qop_config.h>
//AB 4-tensor definition, should be moved to QLA
typedef struct { QLA_F_Complex t4[3][3][3][3]; } QLA_F3_ColorTensor4;
typedef struct { QLA_D_Complex t4[3][3][3][3]; } QLA_D3_ColorTensor4;
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

#define ASQTAD_DSLASH_BEGIN CHECK_INIT
#define ASQTAD_DSLASH_END
#define ASQTAD_FORCE_BEGIN CHECK_INIT
#define ASQTAD_FORCE_END
#define ASQTAD_INVERT_BEGIN CHECK_INIT
#define ASQTAD_INVERT_END

#define HISQ_LINKS_BEGIN CHECK_INIT
#define HISQ_LINKS_END
#define HISQ_FORCE_BEGIN CHECK_INIT
#define HISQ_FORCE_END

#define WILSON_DSLASH_BEGIN CHECK_INIT
#define WILSON_DSLASH_END
#define WILSON_INVERT_BEGIN CHECK_INIT
#define WILSON_INVERT_END

#define DW_INVERT_BEGIN CHECK_INIT
#define DW_INVERT_END

#define QOP_malloc(var, type, num)					\
  (var) = (type *) malloc(num*sizeof(type));				\
  if(!(var)) {								\
    QMP_error("Error: QOP ran out of memory in function %s\n", __func__); \
    exit(1);								\
  }

#define QOP_printf0 if(QDP_this_node==0) printf
#define VERB(PRI, ...) if(QOP_common.verbosity>=QOP_VERB_##PRI) QOP_printf0(__VA_ARGS__)

#ifdef DO_TRACE
#define TRACE QOP_printf0("%s %s %i\n", __FILE__, __func__, __LINE__);
#else
#define TRACE
#endif

#define QOP_error(str) \
  fprintf(stderr, "QOP error: %s\n", str); \
  QDP_abort(1);


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
  int optnum;
  int cgtype;
  int eigcg_nev;
  int eigcg_m;
  int eigcg_numax;
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

/* Currently unused */
typedef struct {
  int inited;
  int fnmat_src_min;
  int veclength;
  double force_filter;
} QOP_hisq_force_t;
extern QOP_hisq_force_t QOP_hisq_ff;

typedef struct {
  int inited;
  int want_deps;
  int want_aux;
  QOP_hisq_unitarize_method_t umethod;
  int reunit_allow_svd;
  int reunit_svd_only;
  double reunit_svd_rel_error;
  double reunit_svd_abs_error;
  int svd_values_info;
} QOP_hisq_links_t;
extern QOP_hisq_links_t QOP_hisq_links;

#ifdef __cplusplus
extern "C" {
#endif

double QOP_time(void);
QOP_status_t QOP_asqtad_invert_init(void);

int
QOP_F3_u3_un_analytic( QOP_info_t *info,
		       QLA_F3_ColorMatrix *V, QLA_F3_ColorMatrix *W );
  
int
QOP_F3_u3_un_der_analytic( QOP_info_t *info, 
			   QLA_F3_ColorMatrix *V, QLA_F3_ColorTensor4 *dwdv, 
			   QLA_F3_ColorTensor4 *dwdagdv );
QLA_F_Complex 
QOP_F3_su3_mat_det( QLA_F3_ColorMatrix *U) ;

int
QOP_D3_u3_un_analytic( QOP_info_t *info,
		       QLA_D3_ColorMatrix *V, QLA_D3_ColorMatrix *W );

int
QOP_D3_u3_un_der_analytic( QOP_info_t *info,
			    QLA_D3_ColorMatrix *V, QLA_D3_ColorTensor4 *dwdv, 
			    QLA_D3_ColorTensor4 *dwdagdv );
QLA_D_Complex 
QOP_D3_su3_mat_det( QLA_D3_ColorMatrix *U) ;



#ifdef __cplusplus
}
#endif

#endif /* _QOP_INTERNAL_H */
