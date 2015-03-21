#ifndef _QOP_INTERNAL_H
#define _QOP_INTERNAL_H

#include <qop.h>
#include <qop_qdp.h>
#include <qop_config.h>

#if QOP_Colors == 'N'

#undef QOP_Nc
#define QOP_Nc NC
#undef QDP_Nc
#define QDP_Nc QOP_Nc
#undef QLA_Nc
#define QLA_Nc QOP_Nc
#define NCPROT int NC,
#define NCPROT1 int NC,
#define NCPROTVOID int NC
#define NCARG NC,
#define NCARG1 NC,
#define NCARGVOID NC

#ifndef QLA_ColorMatrix
#if QOP_Precision == 'F'
#define QLA_ColorVector(x)  QLA_FN_ColorVector (QLA_Nc,(x))
#define QLA_ColorMatrix(x)  QLA_FN_ColorMatrix (QLA_Nc,(x))
#define QLA_DiracFermion(x) QLA_FN_DiracFermion(QLA_Nc,(x))
#else
#define QLA_ColorVector(x)  QLA_DN_ColorVector (QLA_Nc,(x))
#define QLA_ColorMatrix(x)  QLA_DN_ColorMatrix (QLA_Nc,(x))
#define QLA_DiracFermion(x) QLA_DN_DiracFermion(QLA_Nc,(x))
#endif
#define QLA_F_ColorVector(x)  QLA_FN_ColorVector (QLA_Nc,(x))
#define QLA_F_ColorMatrix(x)  QLA_FN_ColorMatrix (QLA_Nc,(x))
#define QLA_F_DiracFermion(x) QLA_FN_DiracFermion(QLA_Nc,(x))
#define QLA_D_ColorVector(x)  QLA_DN_ColorVector (QLA_Nc,(x))
#define QLA_D_ColorMatrix(x)  QLA_DN_ColorMatrix (QLA_Nc,(x))
#define QLA_D_DiracFermion(x) QLA_DN_DiracFermion(QLA_Nc,(x))
#endif

#else

#define NCPROT
#define NCPROT1
#define NCPROTVOID void
#define NCARG
#define NCARG1
#define NCARGVOID
#ifndef QLA_ColorMatrix
#define QLA_ColorMatrix(x) QLA_ColorMatrix x
#define QLA_F_ColorMatrix(x) QLA_F_ColorMatrix x
#define QLA_D_ColorMatrix(x) QLA_D_ColorMatrix x
#endif

#if QOP_Colors == 1
#undef NCPROT1
#define NCPROT1 int NC,
#undef NCARG1
#define NCARG1 NC,
#endif

#endif

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
#define QOP_free(x) free(x)

#define QOP_printf0 if(QDP_this_node==0) printf
#define VERB(PRI, ...) if(QOP_common.verbosity>=QOP_VERB_##PRI) QOP_printf0(__VA_ARGS__)

#ifdef DO_TRACE
#define TRACE QOP_printf0("%s %s %i\n", __FILE__, __func__, __LINE__);
#else
#define TRACE
#endif

#define QOP_error(...) do { \
    printf("QOP error: %s %s %i\n", __FILE__, __func__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    QDP_abort(1); \
  } while(0)

#define QLATYPE_V QLA_ColorVector
#define QLATYPE_D QLA_DiracFermion
#define QLATYPE_M QLA_ColorMatrix

#define QOP_qdp_eq_raw(abbr, qdp, raw, evenodd)			\
  if(evenodd==QOP_EVEN) {					\
    QDP_insert_##abbr(qdp, (QLATYPE_##abbr *)raw, QDP_even);	\
  } else if(evenodd==QOP_ODD) {					\
    QDP_insert_##abbr(qdp, (QLATYPE_##abbr *)raw, QDP_odd);	\
  } else {							\
    QDP_insert_##abbr(qdp, (QLATYPE_##abbr *)raw, QDP_all);	\
  }

#define QOP_raw_eq_qdp(abbr, raw, qdp, evenodd)		      \
  if(evenodd==QOP_EVEN) {				      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_even); \
  } else if(evenodd==QOP_ODD) {				      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_odd);  \
  } else {						      \
    QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_all);  \
  }

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
  int use_fat7_lepage;
} QOP_hisq_links_t;
extern QOP_hisq_links_t QOP_hisq_links;

//double QOP_time(void);
#define QOP_time() QDP_time()
QDP_Subset *QOP_get_sub32(QDP_Lattice *lat);
QOP_status_t QOP_asqtad_invert_init(void);

#define APPLY(m,...) m(__VA_ARGS__)
#define CAT3(a,b,c) a ## b ## c
#define CAT4(a,b,c,d) a ## b ## c ## d

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QOPO(x) QOP_D_##x
#define QOPPO(x) QOP_FD_##x
#define QOPOP(x) QOP_DF_##x
#define QDPO(x) QDP_D_##x
#define QDPPO(x) QDP_FD_##x
#define QDPOP(x) QDP_DF_##x
#define REAL float
#if QOP_Colors == 'N'
#define QOPPC(x) QOP_FN_##x
#else
#define QOPPC(x) APPLY(CAT4,QOP_F,QOP_Colors,_,x)
#endif
#else
#define QOPP(x) QOP_D_##x
#define QOPO(x) QOP_F_##x
#define QOPPO(x) QOP_DF_##x
#define QOPOP(x) QOP_FD_##x
#define QDPO(x) QDP_F_##x
#define QDPPO(x) QDP_DF_##x
#define QDPOP(x) QDP_FD_##x
#define REAL double
#if QOP_Colors == 'N'
#define QOPPC(x) QOP_DN_##x
#else
#define QOPPC(x) APPLY(CAT4,QOP_D,QOP_Colors,_,x)
#endif
#endif

#if QOP_Colors == 'N'
#define QOPFC(x) QOP_FN_##x
#else
#define QOPFC(x) APPLY(CAT4,QOP_F,QOP_Colors,_,x)
#endif
#define QDPPC(x) QDP_##x

//AB HISQ derivatives, rank-4 tensor
#if QOP_Precision == 'F'
#define QLA_ColorTensor4 QLA_F3_ColorTensor4
#define QOP_u3_un_analytic QOP_F3_u3_un_analytic
#define QOP_u3_un_der_analytic QOP_F3_u3_un_der_analytic
#else
#define QLA_ColorTensor4 QLA_D3_ColorTensor4
#define QOP_u3_un_analytic QOP_D3_u3_un_analytic
#define QOP_u3_un_der_analytic QOP_D3_u3_un_der_analytic
#endif

//AB 4-tensor definition, should be moved to QLA
typedef struct { QLA_F_Complex t4[3][3][3][3]; } QLA_F3_ColorTensor4;
typedef struct { QLA_D_Complex t4[3][3][3][3]; } QLA_D3_ColorTensor4;

int
QOP_F3_u3_un_analytic( QOP_info_t *info,
                       QLA_F3_ColorMatrix *V, QLA_F3_ColorMatrix *W );
void
QOP_F3_u3_un_der_analytic( QOP_info_t *info, 
                           QLA_F3_ColorMatrix *V, QLA_F3_ColorTensor4 *dwdv, 
                           QLA_F3_ColorTensor4 *dwdagdv, int *svd_calls, int *ff_counter );
QLA_F_Complex QOP_F3_su3_mat_det( QLA_F3_ColorMatrix *U);
int
QOP_D3_u3_un_analytic( QOP_info_t *info,
                       QLA_D3_ColorMatrix *V, QLA_D3_ColorMatrix *W );
void
QOP_D3_u3_un_der_analytic( QOP_info_t *info,
			   QLA_D3_ColorMatrix *V, QLA_D3_ColorTensor4 *dwdv, 
                           QLA_D3_ColorTensor4 *dwdagdv, int *svd_calls, int *ff_counter );
QLA_D_Complex QOP_D3_su3_mat_det( QLA_D3_ColorMatrix *U);

#if QOP_Precision == 'F'
#  include <qop_f_internal.h>
#else
#  include <qop_d_internal.h>
#endif

#if QOP_Colors == 1
#include <qop_f1_internal.h>
#include <qop_d1_internal.h>
#elif QOP_Colors == 2
#include <qop_f2_internal.h>
#include <qop_d2_internal.h>
#elif QOP_Colors == 3
#include <qop_f3_internal.h>
#include <qop_d3_internal.h>
#elif QOP_Colors == 'N'
#include <qop_fn_internal.h>
#include <qop_dn_internal.h>
#endif

#include <qop_mg_internal.h>

#endif /* _QOP_INTERNAL_H */
