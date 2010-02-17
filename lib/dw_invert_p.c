#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU

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

extern void
(*wilson_dslash)(QOP_FermionLinksDW *fldw,
                 QDP_DiracFermion *dest, QDP_DiracFermion *src,
                 int sign, QDP_Subset subset, int ntmp);

extern void
wilson_dslash0(QOP_FermionLinksDW *fldw,
               QDP_DiracFermion *dest, QDP_DiracFermion *src,
               int sign, QDP_Subset subset, int ntmp);

extern void
wilson_dslash1(QOP_FermionLinksDW *fldw,
               QDP_DiracFermion *dest, QDP_DiracFermion *src,
               int sign, QDP_Subset subset, int ntmp);


extern void
dw_dslash1(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           int sign, QDP_Subset subset, QDP_Subset othersubset,
           QLA_Real m0, QLA_Real M, int ls);

extern void
dw_dslash2(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           QDP_Subset subset, QDP_Subset othersubset,
           QLA_Real m0, QLA_Real M, int ls);

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


extern void Qeo(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
                QDP_DiracFermion *in[], REAL m0, REAL M, int ls);
extern void Qooinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                   REAL m0, REAL M, int ls);
extern void Qeeinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                   REAL m0, REAL M, int ls);
extern void QoeQeeinv(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
                      QDP_DiracFermion *in[], REAL m0, REAL M, int ls);

static QLA_Real gl_m0, gl_M;
static int gl_ls;
static QDP_Subset gl_osubset;
static QOP_FermionLinksDW *gl_fldw;

void
QOPPC(dw_dslash2)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                  QDP_Subset subset)
{
  dw_dslash2(gl_fldw, out, in, subset, gl_osubset, gl_m0, gl_M, gl_ls);
}

void
QOPPC(dw_invert)( QOP_info_t *info,
          QOP_FermionLinksDW *fldw,
            QOP_invert_arg_t *inv_arg,
             QOP_resid_arg_t *res_arg,
                        REAL m0,
                        REAL M,
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
  gl_m0 = m0;
  gl_M = M;
  gl_ls = ls;
  gl_fldw = fldw;

  check_inv_setup(ls);

  for(s=0; s<ls; s++) {
    tiv2[s] = QDP_create_D();
    qdpin[s] = QDP_create_D();
    qdptmp[s] = in[s]->df;
    qdpout[s] = out[s]->df;
  }

  {
#ifdef LU
    QoeQeeinv(fldw, qdpin, qdptmp, m0, M, ls);
    for (s=0; s<ls; s++)
      QDP_D_eq_D_minus_D(qdpin[s], qdptmp[s], qdpin[s], QDP_odd);
    Qooinv(tiv, qdpin, m0, M, ls);
    dw_dslash1(fldw, qdpin, tiv, -1, subset, osubset, m0, M, ls);
#else
    dw_dslash1(fldw, qdpin, qdptmp, -1, subset, osubset, m0, M, ls);
#endif
  }

#ifndef LU
  res_arg->rsqmin /= m0*m0*m0*m0;
#endif
  printf0("begin cgv\n");
  dtime = -QOP_time();

  QOPPC(invert_cg_vD)(QOPPC(dw_dslash2), inv_arg, res_arg,
                      qdpout, qdpin, tiv2, subset, ls);

  dtime += QOP_time();
  printf0("end cgv\n");
#ifndef LU
  res_arg->rsqmin *= m0*m0*m0*m0;
#endif

#ifdef LU
  {
    Qeo(fldw, tiv, qdpout, m0, M, ls);
    for (s=0; s<ls; s++)
      QDP_D_eq_D_minus_D(tiv[s], qdptmp[s], tiv[s], QDP_even);
    Qeeinv(qdpout, tiv, m0, M, ls);
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

