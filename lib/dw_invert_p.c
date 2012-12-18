//#define DO_TRACE
#include <qop_internal.h>

#define LU 1

#define QDP_to_QOP(subset) \
  (subset==QDP_even?QOP_EVEN:(subset==QDP_odd?QOP_ODD:QOP_EVENODD))

// Local globals (for easy passing into the D2 wrapper)
static QLA_Real gl_mq, gl_M5;
static int gl_ls;
static QDP_Subset gl_osubset;
static QOP_FermionLinksDW *gl_fldw;

#if 0
static void
printnrm2(QDP_DiracFermion *x[], int ls, QDP_Subset sub)
{
  QLA_Real nrm2 = 0;
  for(int s=0; s<ls; s++) {
    QLA_Real t;
    QDP_r_eq_norm2_D(&t, x[s], sub);
    nrm2 += t;
  }
  if(QDP_this_node==0) printf("nrm2 = %g\n", nrm2);
}
#else
#define printnrm2(...) (void)0
#endif

// -----------------------------------------------------------------
// Domain-wall inverter

void
QOP_dw_dslash2_wrap(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		    QDP_Subset subset)
{
  QOP_dw_dslash2_qdp(NULL, gl_fldw, gl_M5, gl_mq, out, in, gl_ls,
		     QDP_to_QOP(gl_osubset),QDP_to_QOP(subset));
}

void
QOP_dw_schur2_wrap(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		   QDP_Subset subset)
{
  QOP_dw_schur2_qdp(NULL, gl_fldw, gl_M5, gl_mq, out, in, gl_ls,
		    QDP_to_QOP(subset));
}

void
QOP_dw_invert(QOP_info_t *info,
	      QOP_FermionLinksDW *fldw,
	      QOP_invert_arg_t *inv_arg,
	      QOP_resid_arg_t *res_arg,
	      REAL M5,
	      REAL mq,
	      QOP_DiracFermion *out_pt[],
	      QOP_DiracFermion *in_pt[],
	      int ls)
{
  QDP_DiracFermion *out[ls], *in[ls];
  for(int i=0; i<ls; i++) {
    out[i] = out_pt[i]->df;
    in[i] = in_pt[i]->df;
  }
  QOP_dw_invert_qdp(info, fldw, inv_arg, res_arg, M5, mq, out, in, ls);
}

void
QOP_dw_invert_qdp(QOP_info_t *info,
		  QOP_FermionLinksDW *fldw,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  REAL M5,
		  REAL mq,
		  QDP_DiracFermion *out[],
		  QDP_DiracFermion *in[],
		  int ls)
{
#define NC QDP_get_nc(out[0])
  double dtime;
  double nflop;
  QDP_DiracFermion *qdpin[ls], *tin[ls], *tcg[ls];
  QDP_Subset subset, osubset;

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

  for (int s=0; s<ls; s++) {
    tin[s]  = QDP_create_D();
    tcg[s]  = QDP_create_D();
    qdpin[s] = QDP_create_D();
  }

  printnrm2(in, ls, QDP_all);

// Even-odd precondition the input
#ifdef LU
  QOP_dw_EO_project(fldw, tin, in, M5, mq, ls, QDP_to_QOP(subset));
  printnrm2(tin, ls, subset);
#endif

// Prepare the source for inversion using the normal equations
// (if necessary, which it always is, since only CG is currently implemented)
#ifdef LU
  QOP_dw_schur_qdp(info, fldw, M5, mq, -1, qdpin, tin, ls,
		   QDP_to_QOP(subset));
  // The opposite subset goes through trivially; we'll set it now
  //for (s=0; s<ls; s++) QDP_D_eq_D(out[s], tin[s], osubset);
#else
  QOP_dw_dslash_qdp(info, fldw, M5, mq, -1, qdpin, in, ls,
		    QDP_to_QOP(osubset), QDP_to_QOP(subset));
#endif
  printnrm2(qdpin, ls, subset);

// Invoke the CG inverter on the normal equation
  //printf0("begin cgv\n");
  dtime = -QOP_time();
#ifdef LU
  //printf0("max iterations: %d, relmin: %g\n",inv_arg->max_iter,res_arg->relmin);
  QOP_invert_cg_vD(QOP_dw_schur2_wrap, inv_arg, res_arg,
		   out, qdpin, tcg, subset, ls);
#else
  //res_arg->rsqmin /= (5-M5)*(5-M5)*(5-M5)*(5-M5);
  //printf0("max iterations: %d, relmin: %g\n",inv_arg->max_iter,res_arg->relmin);
  QOP_invert_cg_vD(QOP_dw_dslash2_wrap, inv_arg, res_arg,
		   out, qdpin, tcg, subset, ls);
  //res_arg->rsqmin *= (5-M5)*(5-M5)*(5-M5)*(5-M5);
#endif
  dtime += QOP_time();
  //printf0("end cgv\n");

// Reconstruct the EO solution
#ifdef LU
  QOP_dw_EO_reconstruct(fldw, out, in, M5, mq, ls,
			QDP_to_QOP(subset));
#endif

#if 0
  // calculate final residual
  QOP_dw_dslash_qdp(info, fldw, M5, mq, 1, tcg, out, ls,
		    QOP_EVENODD, QOP_EVENODD);
  //#ifdef LU
  //#else
  //QOPPC(dw_dslash2_qdp)(info, fldw, M5, mq, tcg, out, ls,
  //		QOP_EVENODD, QOP_EVENODD);
  //#endif
  QLA_Real rsq=0, rsqin=0;
  for (int s=0; s<ls; s++) {
    QLA_Real t;
    QDP_D_eq_D_minus_D(tcg[s], in[s], tcg[s], QDP_all);
    QDP_r_eq_norm2_D(&t, tcg[s], QDP_all);
    rsq += t;
    QDP_r_eq_norm2_D(&t, in[s], QDP_all);
    rsqin += t;
  }
  res_arg->final_rsq = rsq/rsqin;
  printf0("final rsq = %g/%g = %g  (%g)\n", rsq, rsqin, res_arg->final_rsq, res_arg->rsqmin);
#endif

  for (int s=0; s<ls; s++) {
    QDP_destroy_D(tin[s]);
    QDP_destroy_D(tcg[s]);
    QDP_destroy_D(qdpin[s]);
  }

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  DW_INVERT_END;
#undef NC
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
		    int nsrc,
		    int ls)
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_invert_multi unimplemented");
  DW_INVERT_END;
}

void
QOP_dw_invert_multi_qdp(QOP_info_t *info,
			QOP_FermionLinksDW *links,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t **res_arg[],
			REAL *m0[],
			REAL *M[],
			int nmass[],
			QDP_DiracFermion ***out_pt[],
			QDP_DiracFermion **in_pt[],
			int nsrc,
			int ls)
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_invert_multi_qdp unimplemented");
  DW_INVERT_END;
}
