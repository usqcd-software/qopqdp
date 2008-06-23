/* adapted from MILC version 6 */

/**********************
** original comments **
**********************/
/******* d_congrad2.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 6 */
/* Wilson fermions */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
   the fermion spinors live on even sites only.  In other words, if
   Dslash_oe is the dslash operator with its source on even sites and
   its result on odd sites, etc.:

   without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
   with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
*/
/**************************
** end original comments **
**************************/

/**************************************************************************
Wilson inverter conventions:

in even-odd blocks (where (a,b) is either (e,o) or (o,e))

D = ( Ca    D_ab )
    ( D_ba  Cb   )

D^-1 = ( A^-1 Ca^-1              -A^-1 Ca^-1 D_ab Cb^-1                 )   
       ( -Cb^-1 D_ba A^-1 Ca^-1  Cb^-1 [1 + D_ba A^-1 Ca^-1 D_ab Cb^-1] )

with A = (1 - Ca^-1 D_ab Cb^-1 D_ba)

with even-odd preconditioning we can write the solution as

D^-1 (x) = ( A^-1 z                )
     (y)   ( Cb^-1 [ y - D_ba A^-1 z ] )

with z = Ca^-1 [ x - D_ab Cb^-1 y ].
***************************************************************************/


#include <string.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU

extern int QOP_wilson_cgtype;
extern int QOP_wilson_eigcg_numax;
extern int QOP_wilson_eigcg_m;
extern int QOP_wilson_eigcg_nev;
extern QDP_DiracFermion *QOPPC(wilson_dslash_get_tmp)(QOP_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);

/* inverter stuff */

static QOP_FermionLinksWilson *gl_flw;
static REAL gl_kappa;
static QOP_evenodd_t gl_eo;
static QDP_DiracFermion *gl_tmp, *gl_tmp2, *ctmp, *gtmp;

#define project(flw, kappa, sign, out, in, eo) \
{ \
  QOP_wilson_diaginv_qdp(NULL, flw, kappa, gl_tmp2, in, oppsub(eo)); \
  QOP_wilson_dslash_qdp(NULL,flw,kappa,sign, out, gl_tmp2, eo, oppsub(eo)); \
  QDP_D_eq_D_minus_D(out, in, out, qdpsub(eo)); \
}

#define reconstruct(flw, kappa, out, soln, src, eo) \
{ \
  QDP_D_eq_D(gl_tmp, soln, qdpsub(oppsub(eo))); \
  QOP_wilson_dslash_qdp(NULL, flw, kappa, 1, out, gl_tmp, eo, oppsub(eo)); \
  QDP_D_eq_D_minus_D(gl_tmp, src, out, qdpsub(eo)); \
  QOP_wilson_diaginv_qdp(NULL, flw, kappa, out, gl_tmp, eo); \
}

#define dslash2(flw, kappa, sign, out, tmp, in, eo) \
{ \
  if(flw->clov==NULL) { \
    QLA_Real mk2 = -4*kappa*kappa; \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, tmp, in, oppsub(eo), eo); \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, out, tmp, eo, oppsub(eo)); \
    QDP_D_eq_r_times_D_plus_D(out, &mk2, out, in, qdpsub(eo)); \
  } else { \
    if(sign==1) { \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, ctmp, in, oppsub(eo), eo); \
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, tmp, ctmp, oppsub(eo)); \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, ctmp, tmp, eo, oppsub(eo)); \
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, out, ctmp, eo); \
    QDP_D_eq_D_minus_D(out, in, out, qdpsub(eo)); \
    } else { \
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, gl_tmp, in, eo); \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, ctmp,gl_tmp,oppsub(eo),eo); \
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, tmp, ctmp, oppsub(eo)); \
    QOP_wilson_dslash_qdp(NULL, flw, kappa, sign, out, tmp, eo, oppsub(eo)); \
    QDP_D_eq_D_minus_D(out, in, out, qdpsub(eo)); \
    } \
  } \
}

void
QOPPC(wilson_invert_d)(QDP_DiracFermion *out, QDP_DiracFermion *in,
		       QDP_Subset subset)
{
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, 1, out, in, gl_eo, gl_eo);
}

void
QOPPC(wilson_invert_dne)(QDP_DiracFermion *out, QDP_DiracFermion *in,
			 QDP_Subset subset)
{
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, 1, gl_tmp, in, gl_eo, gl_eo);
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, -1, out, gl_tmp, gl_eo, gl_eo);
}

void
QOPPC(wilson_invert_d2)(QDP_DiracFermion *out, QDP_DiracFermion *in,
			QDP_Subset subset)
{
#if 0
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, 1,
			gl_tmp3, in, QOP_EVENODD, gl_eo);
  project(gl_flw, gl_kappa, 1, out, gl_tmp3, gl_eo); 
#else
  dslash2(gl_flw, gl_kappa, 1, out, gl_tmp2, in, gl_eo);
#endif
}

void
QOPPC(wilson_invert_d2ne)(QDP_DiracFermion *out, QDP_DiracFermion *in,
			  QDP_Subset subset)
{
#if 0
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, 1,
			gl_tmp3, in, QOP_EVENODD, gl_eo);
  project(gl_flw, gl_kappa, 1, gl_tmp, gl_tmp3, gl_eo); 
  QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, -1,
			gl_tmp3, gl_tmp, QOP_EVENODD, gl_eo);
  project(gl_flw, gl_kappa, -1, out, gl_tmp3, gl_eo); 
#else
  if(gl_flw->clov==NULL) {
    dslash2(gl_flw, gl_kappa, 1, gl_tmp, gl_tmp2, in, gl_eo);
    dslash2(gl_flw, gl_kappa, -1, out, gl_tmp2, gl_tmp, gl_eo);
  } else {
    dslash2(gl_flw, gl_kappa, 1, gtmp, gl_tmp2, in, gl_eo);
    dslash2(gl_flw, gl_kappa, -1, out, gl_tmp2, gtmp, gl_eo);
  }
#endif
}

void
QOP_wilson_invert(QOP_info_t *info,
		  QOP_FermionLinksWilson *flw,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  REAL kappa,
		  QOP_DiracFermion *out,
		  QOP_DiracFermion *in)
{
  double dtime=0;
  double nflop;
  double rsqminold;
  QLA_Real rsq, rsqstop, insq;
  QDP_DiracFermion *qdpin, *qdpout;
  QDP_DiracFermion *cgp, *cgr;
  QDP_Subset insub, cgsub;
  QOP_evenodd_t ineo, cgeo;
  int iter = 0;
  int max_iter_old = inv_arg->max_iter;
  int max_restarts_old = inv_arg->max_restarts;
  int nrestart = -1, max_restarts = inv_arg->max_restarts;
  if(max_restarts<=0) max_restarts = 5;

  WILSON_INVERT_BEGIN;

  if(QOP_wilson_cgtype==2) {
    flw->eigcg.numax = QOP_wilson_eigcg_numax;
    flw->eigcg.m = QOP_wilson_eigcg_m;
    flw->eigcg.nev = QOP_wilson_eigcg_nev;
  }

  ineo = inv_arg->evenodd;
  insub = qdpsub(ineo);

  qdpin = QDP_create_D();
  qdpout = QDP_create_D();
  if(flw->clov!=NULL) ctmp = QDP_create_D();
  if(QOP_wilson_cgtype==0) gtmp = QDP_create_D();

  gl_flw = flw;
  gl_kappa = kappa;

  /* cg has 5 * 48 = 240 flops/site/it */
  /* bicg has 9*4*24 = 864 flops/site/it */
#ifdef LU
  /* MdagM(no clov) -> 2*(2*(144+168*7)+48) = 5376 flops/site */
  /* MdagM(clov)    -> 2*(2*(144+168*7+12*42)+24) = 7344 flops/site */
  if(flw->clov==NULL) {
    nflop = 5376;
  } else {
    nflop = 7344;
  }
  if(QOP_wilson_cgtype==1) {
    nflop += 864;
  } else {
    nflop += 240;
  }
  nflop *= 0.5; /* half for half of the sites */
  cgeo = ineo;
  if(ineo==QOP_EVENODD) cgeo = QOP_EVEN;
#else
  /* MdagM -> 2*((144+168*7)+48) = 2736 flops/site */
  /* clov -> 2*(12*(42-2)) = 960 flops/site */
  nflop = 2736 + 240;
  if(QOP_wilson_cgtype==1) nflop += (864-240);
  if(flw->clov!=NULL) nflop += 960;
  cgeo = QOP_EVENODD;
#endif
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  cgp = QOPPC(wilson_dslash_get_tmp)(flw, oppsub(cgeo), 1);
  cgr = QOPPC(wilson_dslash_get_tmp)(flw, oppsub(cgeo), 2);
  gl_tmp = cgr;

  //printf("test1\n");
  QDP_D_eq_zero(qdpin, QDP_all);
#ifdef LU
  gl_tmp2 = QOPPC(wilson_dslash_get_tmp)(flw, cgeo, 1);
  if(ineo==cgeo) {
    QDP_D_eq_D(qdpin, in->df, insub);
  } else {
    QDP_D_eq_zero(cgp, QDP_all);
    QDP_D_eq_D(cgp, in->df, insub);
    project(flw, kappa, 1, qdpin, cgp, cgeo);
  }
  if(QOP_wilson_cgtype==1) {
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, cgp, qdpin, cgeo);
    QDP_D_eq_D(qdpin, cgp, cgsub);
  } else {
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, cgp, qdpin, cgeo);
    dslash2(flw, kappa, -1, qdpin, gl_tmp2, cgp, cgeo);
  }
#else
  gl_tmp2 = cgr;
  if(QOP_wilson_cgtype==1) {
    QDP_D_eq_D(qdpin, in->df, insub);
  } else {
    QOP_wilson_dslash_qdp(NULL, flw, kappa, -1, qdpin, in->df, cgeo, ineo);
  }
#endif
  if(ineo!=QOP_EVENODD && ineo!=cgeo) {
    QDP_D_eq_zero(qdpout, cgsub);
    reconstruct(flw, kappa, qdpout, out->df, qdpout, oppsub(ineo));
  }
  QDP_D_eq_D(qdpout, out->df, insub);

#if 0
  {
    QDP_r_eq_norm2_D(&rsq, qdpin, cgsub);
    printf("nrm = %g\n", rsq);
    QOP_wilson_diaginv_qdp(NULL, flw, kappa, cgp, qdpin, cgeo);
    QDP_r_eq_norm2_D(&rsq, cgp, cgsub);
    printf("nrm = %g\n", rsq);
    QOP_wilson_dslash_qdp(NULL, flw, kappa, 1, qdpin, cgp, cgeo, cgeo);
    QDP_r_eq_norm2_D(&rsq, qdpin, cgsub);
    printf("nrm = %g\n", rsq);
  }
#endif

  QDP_r_eq_norm2_D(&insq, in->df, insub);
  rsqstop = insq * res_arg->rsqmin;
  VERB(LOW, "WILSON_INVERT: rsqstop = %g\n", rsqstop);
  rsqminold = res_arg->rsqmin;
  res_arg->rsqmin *= 0.9;
  inv_arg->max_restarts = 0;
  do {
    inv_arg->max_iter = max_iter_old - iter;

    dtime -= QOP_time();

#ifdef LU
    if(QOP_wilson_cgtype==2) {
      QOPPC(invert_eigcg_D)(QOPPC(wilson_invert_d2ne), inv_arg, res_arg,
			    qdpout, qdpin, cgp, cgsub, &flw->eigcg);
    } else if(QOP_wilson_cgtype==1) {
      QOPPC(invert_bicgstab_D)(QOPPC(wilson_invert_d2), inv_arg, res_arg,
			       qdpout, qdpin, cgp, cgr, cgsub);
    } else {
      QOPPC(invert_cg_D)(QOPPC(wilson_invert_d2ne), inv_arg, res_arg,
			 qdpout, qdpin, cgp, cgsub);
    }
#else
    if(QOP_wilson_cgtype==2) {
      QOPPC(invert_eigcg_D)(QOPPC(wilson_invert_dne), inv_arg, res_arg,
			    qdpout, qdpin, cgp, cgsub, &flw->eigcg);
    } else if(QOP_wilson_cgtype==1) {
      QOPPC(invert_bicgstab_D)(QOPPC(wilson_invert_d), inv_arg, res_arg,
			       qdpout, qdpin, cgp, cgr, cgsub);
    } else {
      QOPPC(invert_cg_D)(QOPPC(wilson_invert_dne), inv_arg, res_arg,
			 qdpout, qdpin, cgp, cgsub);
    }
#endif

    dtime += QOP_time();
    //printf("finished cg\n");

    //QDP_r_eq_norm2_D(&rsq, qdpout, cgsub);
    //printf("nrm = %g\n", rsq);
#ifdef LU
    if(ineo==QOP_EVENODD) {
      reconstruct(flw, kappa, qdpout, qdpout, in->df, oppsub(cgeo));
      //reconstruct(flw, kappa, qdpout, qdpout, qdpin, oppsub(cgeo));
    } else {
      reconstruct(flw, kappa, qdpout, qdpout, qdpin, oppsub(cgeo));
    }
#endif
    // get final residual
    //QDP_r_eq_norm2_D(&rsq, qdpout, QDP_all);
    //printf("nrm = %g\n", rsq);
    QOP_wilson_dslash_qdp(NULL, flw, kappa, 1, cgr, qdpout, ineo, QOP_EVENODD);
    QDP_D_meq_D(cgr, in->df, insub);
    QDP_r_eq_norm2_D(&rsq, cgr, insub);
    //printf("%i %i rsq = %g\tprec rsq = %g\trsqstop = %g\n", nrestart,
    //res_arg->final_iter, rsq, res_arg->final_rsq, rsqstop);
    res_arg->rsqmin = 0.9*res_arg->final_rsq*rsqstop/rsq;
    iter += res_arg->final_iter;
    VERB(LOW, "WILSON_INVERT: iter %i rsq = %g\n", iter, rsq);
  } while((rsq>rsqstop)&&(nrestart++<max_restarts));

  QDP_D_eq_D(out->df, qdpout, insub);

  QDP_destroy_D(qdpin);
  QDP_destroy_D(qdpout);
  if(flw->clov!=NULL) QDP_destroy_D(ctmp);
  if(QOP_wilson_cgtype==0) QDP_destroy_D(gtmp);

  inv_arg->max_iter = max_iter_old;
  inv_arg->max_restarts = max_restarts_old;
  res_arg->rsqmin = rsqminold;
  res_arg->final_iter = iter;
  res_arg->final_rsq = rsq/insq;
  res_arg->final_restart = nrestart;

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  WILSON_INVERT_END;
}

void
QOP_wilson_invert_multi(QOP_info_t *info,
			QOP_FermionLinksWilson *links,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t **res_arg[],
			REAL *kappas[],
			int nkappa[],
			QOP_DiracFermion **out_pt[],
			QOP_DiracFermion *in_pt[],
			int nsrc)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_invert_multi unimplemented");
  WILSON_INVERT_END;
}
