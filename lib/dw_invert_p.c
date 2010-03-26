#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU 1
#define SDC_DEBUG 0

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

#define QDP_to_QOP(subset) \
  (subset==QDP_even?QOP_EVEN:(subset==QDP_odd?QOP_ODD:QOP_EVENODD))

// External variables and functions
extern int QOP_dw_initQ;
extern int QOP_dw_style;
extern int dslash_initQ;

extern void
QOPPC(dw_schur2_qdp)(QOP_info_t *info, QOP_FermionLinksDW *fldw,
                     QLA_Real M5, QLA_Real mq,
                     QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                     int ls, QOP_evenodd_t eo);
extern void
QOPPC(dw_schur_qdp)(QOP_info_t *info, QOP_FermionLinksDW *fldw,
                    QLA_Real M5, QLA_Real mq, int sign,
                    QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                    int ls, QOP_evenodd_t eo);
extern void
QOPPC(dw_EO_project)(QOP_FermionLinksDW *fldw,
                     QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                     QLA_Real M5, QLA_Real mq, int ls, QOP_evenodd_t eo);
extern void
QOPPC(dw_EO_reconstruct)(QOP_FermionLinksDW *fldw,
                         QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                         QLA_Real M5, QLA_Real mq, int ls, QOP_evenodd_t eo);

// Local globals (for easy passing into the D2 wrapper)
static QLA_Real gl_mq, gl_M5;
static int gl_ls;
static QDP_Subset gl_osubset;
static QOP_FermionLinksDW *gl_fldw;


// -----------------------------------------------------------------
// Domain-wall inverter

void
QOPPC(dw_dslash2_wrap)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                       QDP_Subset subset)
{
  QOPPC(dw_dslash2_qdp)(NULL, gl_fldw, gl_M5, gl_mq, out, in, gl_ls,
                        QDP_to_QOP(gl_osubset),QDP_to_QOP(subset));
}

void
QOPPC(dw_schur2_wrap)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                       QDP_Subset subset)
{
  QOPPC(dw_schur2_qdp)(NULL, gl_fldw, gl_M5, gl_mq, out, in, gl_ls,
                       QDP_to_QOP(subset));
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
  QDP_DiracFermion *qdptmp[ls], *qdpin[ls], *qdpout[ls], *tin[ls], *tcg[ls];
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

  // Set globals that will be used by the wrapper functions
  gl_osubset = osubset;
  gl_mq = mq;
  gl_M5 = M5;
  gl_ls = ls;
  gl_fldw = fldw;

  for (s=0; s<ls; s++) {
    tin[s]  = QDP_create_D();
    tcg[s]  = QDP_create_D();
    qdpin[s] = QDP_create_D();
    qdptmp[s] = in[s]->df;
    qdpout[s] = out[s]->df;
  }

// Even-odd precondition the input
#ifdef LU
  QOPPC(dw_EO_project)(fldw, tin, qdptmp, M5, mq, ls, QDP_to_QOP(subset));
#endif

// Prepare the source for inversion using the normal equations
// (if necessary, which it always is, since only CG is currently implemented)
#ifdef LU
  QOPPC(dw_schur_qdp)(info, fldw, M5, mq, -1, qdpin, tin, ls,
                      QDP_to_QOP(subset));
#else
  QOPPC(dw_dslash_qdp)(info, fldw, M5, mq, -1, qdpin, qdptmp, ls,
                       QDP_to_QOP(osubset), QDP_to_QOP(subset));
#endif

// Invoke the CG inverter on the normal equation
  printf0("begin cgv\n");
  dtime = -QOP_time();

#ifdef LU
  printf0("max iterations: %d, relmin: %g\n",inv_arg->max_iter,res_arg->relmin);
  QOPPC(invert_cg_vD)(QOPPC(dw_schur2_wrap), inv_arg, res_arg,
                      qdpout, qdpin, tcg, subset, ls);  
#else
  res_arg->rsqmin /= (5-M5)*(5-M5)*(5-M5)*(5-M5);
  printf0("max iterations: %d, relmin: %g\n",inv_arg->max_iter,res_arg->relmin);
  QOPPC(invert_cg_vD)(QOPPC(dw_dslash2_wrap), inv_arg, res_arg,
                      qdpout, qdpin, tcg, subset, ls);
  res_arg->rsqmin *= (5-M5)*(5-M5)*(5-M5)*(5-M5);
#endif
  dtime += QOP_time();
  printf0("end cgv\n");

// Reconstruct the EO solution
#ifdef LU
  QOPPC(dw_EO_reconstruct)(fldw, qdpout, qdpout, M5, mq, ls,
                           QDP_to_QOP(subset));
#endif

  for (s=0; s<ls; s++) {
    QDP_destroy_D(tin[s]);
    QDP_destroy_D(tcg[s]);
    QDP_destroy_D(qdpin[s]);
  }

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  DW_INVERT_END;
}

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
