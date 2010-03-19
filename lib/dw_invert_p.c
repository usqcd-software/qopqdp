#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#undef LU

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_dw_initQ;
extern int QOP_dw_style;
extern int QOP_dw_nsvec;
extern int QOP_dw_nvec;
extern int dslash_initQ;

// -----------------------------------------------------------------
// Temporary variable management

static int invert_initQ = 0; // Is dslash set up?
// The temporaries are born here
static QDP_DiracFermion **tiv=NULL;

#define QDP_to_QOP(subset) \
  (subset==QDP_even?QOP_EVEN:(subset==QDP_odd?QOP_ODD:QOP_EVENODD))

#define check_inv_setup(ls) \
{ \
  if( (!invert_initQ) ) { \
    reset_inv_temps(ls); \
  } \
}

void
free_inv_temps(int ls)
{
  if (invert_initQ) {    
    for (int s=0; s<ls; s++) {
      QDP_destroy_D(tiv[s]);
    }
    free(tiv);
  }

  invert_initQ = 0;
}

void
reset_inv_temps(int ls)
{
  free_inv_temps(ls);
  
  tiv = (QDP_DiracFermion**) malloc(ls*sizeof(QDP_DiracFermion*));
  for (int s=0; s<ls; s++) {
    tiv[s] = QDP_create_D();
  }

  invert_initQ = 1;
}

// -----------------------------------------------------------------
// Domain-wall inverter

void
QOP_dw_invert_multi(QOP_info_t *info,
            QOP_FermionLinksDW *links,
              QOP_invert_arg_t *inv_arg,
               QOP_resid_arg_t **res_arg[],
                          REAL *m0[],
                          REAL *M[],
                           int nmass[],
              QOP_DiracFermion ***out_pt[],
              QOP_DiracFermion **in_pt[],
                           int ls,
                           int nsrc)
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_invert_multi unimplemented");
  DW_INVERT_END;
}


extern void QOPPC(Qeo)(QOP_FermionLinksDW *fldw,
                       QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                       REAL m0, REAL M, int ls);
extern void QOPPC(Qooinv)(QOP_FermionLinksDW *fldw,
                          QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                          REAL m0, REAL M, int ls);
extern void QOPPC(Qeeinv)(QOP_FermionLinksDW *fldw,
                          QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                          REAL m0, REAL M, int ls);
extern void QOPPC(QoeQeeinv)(QOP_FermionLinksDW *fldw,
                             QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                             REAL m0, REAL M, int ls);

static QLA_Real gl_mq, gl_M5;
static int gl_ls;
static QDP_Subset gl_osubset;
static QOP_FermionLinksDW *gl_fldw;

void
QOPPC(dw_dslash2_wrap)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                       QDP_Subset subset)
{
  QOPPC(dw_dslash2_qdp)(NULL, gl_fldw, gl_M5, gl_mq, out, in, gl_ls,
                        QDP_to_QOP(gl_osubset),QDP_to_QOP(subset));
}

void
QOPPC(dw_invert)( QOP_info_t *info,
          QOP_FermionLinksDW *fldw,
            QOP_invert_arg_t *inv_arg,
             QOP_resid_arg_t *res_arg,
                        REAL M5,
                        REAL mq,
            QOP_DiracFermion *out[],
            QOP_DiracFermion *in[],
                         int ls)
{
  double dtime;
  double nflop;
  QDP_DiracFermion *qdptmp[ls], *qdpin[ls], *qdpout[ls], *tiv2[ls];
  QDP_Subset subset, osubset;
  int s;

  DW_INVERT_BEGIN;

  /* cg has 5 * 48 *ls = 240*ls flops/site/it */
#ifdef LU
  /*              mf     D_w   1st Qxx1  1     */
  /* MdagM -> 2*(-48 + ((168*8-24)+4*24+12)*ls) = -96 + 2856*ls flops/site */
  nflop = -96 + 2976*ls; /* cg only on 1/2 the sites */
  subset = QDP_odd;
  osubset = QDP_even;
#else
  /*              mf       D_w   m0 .5 Ds     */
  /* MdagM -> 2*(2*24 + ((168*8)+24+24+24)*ls) = 96 + 2832*ls flops/site */
  nflop = 96 + 3072*ls;
  subset = QDP_all;
  osubset = QDP_all;
#endif

  gl_osubset = osubset;
  gl_mq = mq;
  gl_M5 = M5;
  gl_ls = ls;
  gl_fldw = fldw;
  
  REAL M0 = M5-5;

  check_inv_setup(ls);

  for(s=0; s<ls; s++) {
    tiv2[s] = QDP_create_D();
    qdpin[s] = QDP_create_D();
    qdptmp[s] = in[s]->df;
    qdpout[s] = out[s]->df;
  }

  {
#ifdef LU
    QOPPC(QoeQeeinv)(fldw, qdpin, qdptmp, M0, mq, ls);
    for (s=0; s<ls; s++)
      QDP_D_eq_D_minus_D(qdpin[s], qdptmp[s], qdpin[s], QDP_odd);
    QOPPC(Qooinv)(fldw, tiv, qdpin, M0, mq, ls);
    QOPPC(dw_dslash_qdp)(info, fldw, M5, mq, -1, qdpin, tiv, ls,
                         QDP_to_QOP(osubset), QDP_to_QOP(subset));
#else
    QOPPC(dw_dslash_qdp)(info, fldw, M5, mq, -1, qdpin, qdptmp, ls,
                         QDP_to_QOP(osubset), QDP_to_QOP(subset));
#endif
  }

#ifndef LU
 res_arg->rsqmin /= M0*M0*M0*M0; // ??? Why?
#endif
  printf0("begin cgv\n");
  dtime = -QOP_time();
  QOP_common.verbosity = QOP_VERB_DEBUG;
  printf("max iterations: %d, relmin: %g\n",inv_arg->max_iter,res_arg->relmin);
  QOPPC(invert_cg_vD)(QOPPC(dw_dslash2_wrap), inv_arg, res_arg,
                      qdpout, qdpin, tiv2, subset, ls);

  dtime += QOP_time();
  printf0("end cgv\n");
#ifndef LU
  res_arg->rsqmin *= M0*M0*M0*M0;
#endif

#ifdef LU
  {
    QOPPC(Qeo)(fldw, tiv, qdpout, M0, mq, ls);
    for (s=0; s<ls; s++)
      QDP_D_eq_D_minus_D(tiv[s], qdptmp[s], tiv[s], QDP_even);
    QOPPC(Qeeinv)(fldw, qdpout, tiv, M0, mq, ls);
  }
#endif

  for (s=0; s<ls; s++) {
    QDP_destroy_D(tiv2[s]);
    QDP_destroy_D(qdpin[s]);
  }
  free_inv_temps(ls);

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  DW_INVERT_END;
}

