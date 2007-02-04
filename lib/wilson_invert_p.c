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

//#define printf0 QOP_printf0
#define printf0(...)

#define LU
//#define CLOV_FUNC

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "phi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include <string.h>
//#define QDP_PROFILE
#include <qop_internal.h>

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_wilson_inited;
extern int QOP_wilson_style;
extern int QOP_wilson_nsvec;
extern int QOP_wilson_nvec;
extern int QOP_wilson_cgtype;

static int old_style=-1;
static int old_nsvec=-1;
static int old_nvec=-1;

static int congrad_setup = 0;
static QDP_HalfFermion *htemp[4][20];
static QDP_DiracFermion *dtemp[4][12];
static QDP_DiracFermion *ttt, *tt1, *tt2;

static void
free_temps(QOP_FermionLinksWilson *flw)
{
  if(congrad_setup) {
    int i, j;

    QDP_destroy_D(ttt);
    QDP_destroy_D(tt1);
    QDP_destroy_D(tt2);

    if(shiftd_style(old_style)) {
      for(i=0; i<4; i++) {
	for(j=0; j<12; j++) {
	  QDP_destroy_D(dtemp[i][j]);
	}
      }
    } else {
      for(i=0; i<4; i++) {
	for(j=0; j<20; j++) {
	  QDP_destroy_H(htemp[i][j]);
	}
      }
    }
  }
  congrad_setup = 0;
}

static void
double_store(QOP_FermionLinksWilson *flw)
{
  if( dblstore_style(QOP_wilson_style) && (!flw->dblstored) ) {
    int i;
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<4; i++) {
      QDP_M_eq_sM(m, flw->links[i], QDP_neighbor[i], QDP_backward, QDP_all);
      QDP_M_eq_Ma(flw->bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    flw->dblstored = 1;
  }
}

static void
reset_temps(QOP_FermionLinksWilson *flw)
{
  int i, j;

  if(QOP_wilson_style!=old_style) {
    if(!dblstore_style(QOP_wilson_style)) {
      if(congrad_setup) {
	if(flw->dblstored) {
	  for(i=0; i<4; i++) {
	    QDP_destroy_M(flw->bcklinks[i]);
	  }
	}
      }
    } else {
      for(i=0; i<4; i++) {
        flw->bcklinks[i] = QDP_create_M();
      }
      for(i=0; i<4; i++) {
        flw->dbllinks[2*i] = flw->links[i];
        flw->dbllinks[2*i+1] = flw->bcklinks[i];
      }
    }
    flw->dblstored = 0;
  }
  double_store(flw);

  free_temps(flw);

  ttt = QDP_create_D();
  tt1 = QDP_create_D();
  tt2 = QDP_create_D();

  if(shiftd_style(QOP_wilson_style)) {
    for(i=0; i<4; i++) {
      for(j=0; j<12; j++) {
	dtemp[i][j] = QDP_create_D();
      }
    }
  } else {
    for(i=0; i<4; i++) {
      for(j=0; j<20; j++) {
	htemp[i][j] = QDP_create_H();
      }
    }
  }
  congrad_setup = 1;
}


/* link routines */

static void
get_clov(QDP_DiracPropagator *clov, QDP_ColorMatrix *links[], double csw)
{
  QLA_DiracPropagator p;
  int i, j;
  QLA_P_eq_zero(&p);
  for(i=0; i<QLA_Nc; i++) for(j=0; j<QLA_Ns; j++)
    QLA_c_eq_r(QLA_elem_P(p, i, j, i, j), 1);
  QDP_P_eq_p(clov, &p, QDP_all);
}

static void
get_clovinv(REAL *clovinv, REAL *clov)
{
  memcpy(clovinv, clov, QDP_sites_on_node*2*6*6*sizeof(REAL));
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
{
  QOP_FermionLinksWilson *flw;
  QOP_GaugeField *gf;

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_wilson_convert_L_from_qdp(gf->links, NULL);

  QOP_malloc(flw->clov, REAL, QDP_sites_on_node*2*6*6);
  QOP_malloc(flw->clovinv, REAL, QDP_sites_on_node*2*6*6);
  memcpy(flw->clov, clov, QDP_sites_on_node*2*6*6*sizeof(REAL));
  get_clovinv(flw->clovinv, flw->clov);
  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qdpclov = NULL;
  flw->qopgf = gf;

  return flw;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			   QOP_GaugeField *gauge)
{
  QOP_FermionLinksWilson *flw;
  QDP_DiracPropagator *clov;
  if(coeffs->clov_c==0.) {
    clov = NULL;
  } else {
    clov = QDP_create_P();
    get_clov(clov, gauge->links, coeffs->clov_c);
  }
  flw = QOP_wilson_create_L_from_qdp(gauge->links, clov);
  if(clov) QDP_destroy_P(clov);
  return flw;
}

void
QOP_wilson_extract_L_to_raw(REAL *links[], REAL *clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_extract_L_to_raw unimplemented.");
}

void
QOP_wilson_destroy_L(QOP_FermionLinksWilson *flw)
{
  int i;

  if(flw->qopgf) {
    QOP_destroy_G(flw->qopgf);
  } else {
    for(i=0; i<4; i++) QDP_destroy_M(flw->links[i]);
    free(flw->links);
  }
  if(flw->dblstored) {
    for(i=0; i<4; i++) QDP_destroy_M(flw->bcklinks[i]);
  }
  QDP_destroy_D(flw->cgp);
  if(flw->qdpclov) {
    QDP_destroy_P(flw->qdpclov);
  }
  if(flw->clov) {
    free(flw->clov);
    free(flw->clovinv);
  }
  free(flw->bcklinks);
  free(flw->dbllinks);
  free(flw);
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], REAL *clov,
			      QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_convert_L_from_raw unimplemented");
  return NULL;
}

void
QOP_wilson_convert_L_to_raw(REAL ***links, REAL **clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_convert_L_to_raw unimplemented");
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			    QOP_GaugeField *gauge)
{
  QOP_error("QOP_wilson_convert_L_from_G unimplemented");
  return NULL;
}

QOP_GaugeField *
QOP_wilson_convert_L_to_G(QOP_FermionLinksWilson *links)
{
  QOP_error("QOP_wilson_convert_L_to_G unimplemented");
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_qdp(QDP_ColorMatrix *links[],
			     QDP_DiracPropagator *clov)
{
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_wilson_convert_L_from_qdp(newlinks, clov);
  flw->qdpclov = NULL;

  return flw;
}

void
QOP_wilson_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QDP_DiracPropagator *clov,
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_extract_L_to_qdp unimplemented");
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_qdp(QDP_ColorMatrix *links[],
			      QDP_DiracPropagator *clov)
{
  QOP_FermionLinksWilson *flw;
  int i;

  QOP_malloc(flw, QOPPC(FermionLinksWilson), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);
  if(clov!=NULL) {
    int size = QDP_sites_on_node*2*6*6;
    QOP_malloc(flw->clov, REAL, size);
    QOP_malloc(flw->clovinv, REAL, size);
    {
      QLA_DiracPropagator *dp;
      int x, b, i, ic, j, jc, is, js, k=0;
      dp = QDP_expose_P(clov);
      for(x=0; x<QDP_sites_on_node; x++) {
	for(b=0; b<2; b++) { // two chiral blocks
	  // first the diagonal
	  for(i=0; i<6; i++) {
	    ic = i/2;
	    is = 2*b + i%2;
	    flw->clov[k++] = QLA_real(QLA_elem_P(dp[x], ic, is, ic, is));
	  }
	  // now the offdiagonal
	  for(i=0; i<6; i++) {
	    ic = i/2;
	    is = 2*b + i%2;
	    for(j=i+1; j<6; j++) {
	      QLA_Complex z1, z2;
	      jc = j/2;
	      js = 2*b + j%2;
	      //QLA_c_eq_c_plus_ca(z1, QLA_elem_P(dp[x], ic, is, jc, js),
	      //                   QLA_elem_P(dp[x], jc, js, ic, is));
	      QLA_c_eq_c(z1, QLA_elem_P(dp[x], ic, is, jc, js));
	      QLA_c_peq_ca(z1, QLA_elem_P(dp[x], jc, js, ic, is));
	      QLA_c_eq_r_times_c(z2, 0.5, z1);
	      flw->clov[k++] = QLA_real(z2);
	      flw->clov[k++] = QLA_imag(z2);
	    }
	  }
	}
      }
      QDP_reset_P(clov);
    }
    get_clovinv(flw->clovinv, flw->clov);
  } else {
    flw->clov = NULL;
    flw->clovinv = NULL;
  }

  flw->dblstored = 0;
  for(i=0; i<4; i++) {
    flw->links[i] = links[i];
  }
  flw->cgp = QDP_create_D();

  if( (!congrad_setup) ||
      (QOP_wilson_style != old_style) ||
      (QOP_wilson_nsvec != old_nsvec) ||
      (QOP_wilson_nvec != old_nvec) ) {
    reset_temps(flw);
    old_style = QOP_wilson_style;
    old_nsvec = QOP_wilson_nsvec;
    old_nvec = QOP_wilson_nvec;
  }

  double_store(flw);

  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qopgf = NULL;
  flw->qdpclov = clov;
  return flw;
}

void
QOP_wilson_convert_L_to_qdp(QDP_ColorMatrix ***links,
			    QDP_DiracPropagator **clov,
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_convert_L_to_qdp unimplemented");
}


/* inverter stuff */

static void
clov(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out, QLA_Real *mkappa,
     QDP_DiracFermion *dsl, QDP_DiracFermion *in, QDP_Subset subset);

static void (*dslash_special_qdp)(QOP_FermionLinksWilson *flw,
				  QDP_DiracFermion *dest,
				  QDP_DiracFermion *src,
				  int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash0(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash1(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp);

static void wilson_mdslash1(QOP_FermionLinksWilson *flw,
			    QDP_DiracFermion *out, QDP_DiracFermion *in,
			    int sign, QDP_Subset subset,
			    QDP_Subset othersubset, QLA_Real mkappa);

static void wilson_mdslash2(QOP_FermionLinksWilson *flw,
			    QDP_DiracFermion *out, QDP_DiracFermion *in,
			    QDP_Subset subset, QDP_Subset othersubset,
			    QLA_Real mkappa);

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
  QOP_error("QOP_wilson_invert_multi unimplemented");
}



static QLA_Real gl_mkappa;
static QDP_Subset gl_osubset;
static QOP_FermionLinksWilson *gl_flw;

void
QOPPC(wilson_dslash2)(QDP_DiracFermion *out, QDP_DiracFermion *in,
		      QDP_Subset subset)
{
  wilson_mdslash2(gl_flw, out, in, subset, gl_osubset, gl_mkappa);
}

void
QOPPC(wilson_invert)(QOP_info_t *info,
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
  QLA_Real mkappa, rsq, rsqstop, insq;
  QDP_DiracFermion *qdpin;
  QDP_Subset subset, osubset;
  int iter = 0;
  int nrestart=0, max_restarts=inv_arg->max_restarts;
  if(max_restarts<=0) max_restarts = 5;

  /* cg has 5 * 48 = 240 flops/site/it */
  /* bicg has 9*4*24 = 864 flops/site/it */
#ifdef LU
  /* MdagM -> 2*(2*(144+168*7)+48) = 5376 flops/site */
  if(QOP_wilson_cgtype==1) {
    nflop = 0.5 * 6240;  /* half for half of the sites */
  } else {
    nflop = 0.5 * 5616;  /* half for half of the sites */
  }
  mkappa = -kappa*kappa;
  subset = QDP_even;
  osubset = QDP_odd;
#else
  /* MdagM -> 2*((144+168*7)+48) = 2736 flops/site */
  /* clov -> 2*(12*44) = 1056 flops/site */
  nflop = 2736 + 240;
  if(QOP_wilson_cgtype==1) nflop += (864-240);
  if(flw->clov!=NULL) nflop += 1056;
  mkappa = -kappa;
  subset = QDP_all;
  osubset = QDP_all;
#endif

  gl_osubset = osubset;
  gl_mkappa = mkappa;
  gl_flw = flw;

  if( (QOP_wilson_style != old_style) ||
      (QOP_wilson_nsvec != old_nsvec) ||
      (QOP_wilson_nvec != old_nvec) ) {
    reset_temps(flw);
    old_style = QOP_wilson_style;
    old_nsvec = QOP_wilson_nsvec;
    old_nvec = QOP_wilson_nvec;
  }

  if(dblstore_style(QOP_wilson_style)) {
    dslash_special_qdp = wilson_dslash1;
  } else {
    dslash_special_qdp = wilson_dslash0;
  }

  qdpin = QDP_create_D();

  QDP_r_eq_norm2_D(&insq, in->df, QDP_all);
  //printf("insq = %g\n", insq);
  rsqstop = insq * res_arg->rsqmin;
  rsqminold = res_arg->rsqmin;
  res_arg->rsqmin *= 0.5;
  do {

    {
      QLA_Real kappa2 = 2.0*kappa;
#ifdef LU
      QDP_D_eq_D(tt1, in->df, osubset);
      dslash_special_qdp(flw, ttt, tt1, 1, subset, 1);
      //kappa = -kappa;
      QDP_D_eq_r_times_D_plus_D(ttt, &kappa, ttt, in->df, subset);
      //kappa = -kappa;
      QDP_D_eq_r_times_D(ttt, &kappa2, ttt, subset);
      QDP_D_eq_r_times_D(qdpin, &kappa2, in->df, osubset);
      if(QOP_wilson_cgtype==1) {
	QDP_D_eq_D(qdpin, ttt, subset);
      } else {
	wilson_mdslash1(flw, qdpin, ttt, -1, subset, osubset, mkappa);
      }
#else
      if(QOP_wilson_cgtype==1) {
	QDP_D_eq_r_times_D(qdpin, &kappa2, in->df, subset);
      } else {
	QDP_D_eq_r_times_D(ttt, &kappa2, in->df, subset);
	wilson_mdslash1(flw, qdpin, ttt, -1, subset, osubset, mkappa);
      }
#endif
    }

    //printf("starting cg\n");
    dtime -= QOP_time();

    if(QOP_wilson_cgtype==1) {
      QOPPC(invert_bicgstab_D)(QOPPC(wilson_dslash2), inv_arg, res_arg,
			       out->df, qdpin, flw->cgp, ttt, subset);
    } else {
      QOPPC(invert_cg_D)(QOPPC(wilson_dslash2), inv_arg, res_arg,
			 out->df, qdpin, flw->cgp, subset);
    }

    dtime += QOP_time();
    //printf("finished cg\n");

#ifdef LU
    {
      QDP_D_eq_D(ttt, out->df, subset);
      dslash_special_qdp(flw, tt2, ttt, 1, osubset, 2);
      QDP_D_eq_r_times_D_plus_D(ttt, &kappa, tt2, qdpin, osubset);
      QDP_D_eq_D(out->df, ttt, QDP_all);
    }
#endif

    // get final residual
    {
      QLA_Real mk = -kappa;
      QLA_Real m2ki = -0.5/kappa;
#ifdef LU
      QDP_D_eq_D(ttt, out->df, QDP_even);
      dslash_special_qdp(flw, tt2, ttt, 1, QDP_odd, 2);
      QDP_D_eq_D(tt1, out->df, QDP_odd);
      dslash_special_qdp(flw, tt2, tt1, 1, QDP_even, 1);
#else
      QDP_D_eq_D(ttt, out->df, QDP_all);
      dslash_special_qdp(flw, tt2, ttt, 1, QDP_all, 1);
#endif
      //QDP_D_eq_r_times_D_plus_D(ttt, &mk, tt2, out->df, QDP_all);
      clov(flw, ttt, &mk, tt2, out->df, QDP_all);
      QDP_D_eq_r_times_D_plus_D(ttt, &m2ki, ttt, in->df, QDP_all);
      QDP_r_eq_norm2_D(&rsq, ttt, QDP_all);
      //printf("rsq = %g\tprec rsq = %g\trsqstop = %g\n",
      //     rsq, res_arg->final_rsq, rsqstop);
    }
    res_arg->rsqmin *= 0.5*rsqstop/rsq;
    iter += res_arg->final_iter;
  } while((rsq>rsqstop)&&(nrestart++<max_restarts));
  QDP_destroy_D(qdpin);

  res_arg->rsqmin = rsqminold;
  res_arg->final_iter = iter;
  res_arg->final_rsq = rsq;

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}

QLA_Real clov_mkappa, *clov_clov;
QLA_DiracFermion *clov_dsl, *clov_in;

#define cmplx(x) (*((QLA_Complex *)(&(x))))

#ifdef CLOV_FUNC
static void
clov_func(QLA_DiracFermion *out, int coords[])
{
  int b, i, j, ic, jc, is, js, k, x, xb;
  QLA_DiracFermion *dsl;

  x = QDP_index(coords);  // site offset
  if(clov_dsl!=NULL) dsl = &clov_dsl[x]; else dsl = out;

  for(b=0; b<2; b++) {
    xb = 36*(2*x+b);  // chiral block offset (in REALs)
    for(i=0; i<6; i++) {
      QLA_Complex z;
      ic = i/2;
      is = 2*b + i%2;

      QLA_c_eq_r(z, 0.);

      //#if 0
      // lower triangular part comes from adjoint of upper
      k = xb + 6 + 2*(i-1); // block + skip diag + cmplx elem
      for(j=0; j<i; j++) {
	jc = j/2;
	js = 2*b + j%2;
	QLA_c_peq_ca_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	k += 2*(4 - j);
      }
      // diagonal part
      QLA_c_peq_r_times_c(z, clov_clov[xb], QLA_elem_D(clov_in[x],ic,is));
      // upper triangular part
      for(j=i+1; j<6; j++) {
	jc = j/2;
	js = 2*b + j%2;
	k += 2;
	QLA_c_peq_c_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
      }
      //#endif

      QLA_c_peq_r_times_c(z, clov_mkappa, QLA_elem_D(*dsl, ic, is));
      QLA_c_eq_c_plus_c(QLA_elem_D(*out, ic, is), z, QLA_elem_D(clov_in[x],ic,is));
    }
  }
}
#endif

static void
apply_clov(REAL *clov, QDP_DiracFermion *out, QLA_Real *mkappa,
	   QDP_DiracFermion *dsl, QDP_DiracFermion *in, QDP_Subset subset)
{
#ifdef CLOV_FUNC
  clov_in = QDP_expose_D(in);
  clov_mkappa = *mkappa;
  clov_clov = clov;
  if(dsl==out) clov_dsl = NULL; else clov_dsl = QDP_expose_D(dsl);
  QDP_D_eq_func(out, clov_func, subset);
  if(dsl!=out) QDP_reset_D(dsl);
  QDP_reset_D(in);
#else
  QLA_DiracFermion *clov_out;
  clov_out = QDP_expose_D(out);
  clov_in = QDP_expose_D(in);
  clov_mkappa = *mkappa;
  clov_clov = clov;
  if(dsl==out) clov_dsl = clov_out;
  else clov_dsl = QDP_expose_D(dsl);
  {
    int x, start, end;
    if(subset==QDP_odd) start = QDP_subset_len(QDP_even);
    else start = 0;
    end = start + QDP_subset_len(subset);
    for(x=start; x<end; x++) {
      int b;
      QLA_DiracFermion *dsl;
      dsl = &clov_dsl[x];
      for(b=0; b<2; b++) {
	int xb;
	xb = 36*(2*x+b);  // chiral block offset (in REALs)
#if 0 // loop version
	int i, j, ic, jc, is, js, k;
	for(i=0; i<6; i++) {
	  QLA_Complex z;
	  ic = i/2;
	  is = 2*b + i%2;

	  QLA_c_eq_r(z, 0.);
	  //#if 0
	  // lower triangular part comes from adjoint of upper
	  k = xb + 6 + 2*(i-1); // block + skip diag + cmplx elem
	  for(j=0; j<i; j++) {
	    jc = j/2;
	    js = 2*b + j%2;
	    QLA_c_peq_ca_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	    k += 2*(4 - j);
	  }
	  // diagonal part
	  QLA_c_peq_r_times_c(z, clov_clov[xb], QLA_elem_D(clov_in[x],ic,is));
	  // upper triangular part
	  for(j=i+1; j<6; j++) {
	    jc = j/2;
	    js = 2*b + j%2;
	    k += 2;
	    QLA_c_peq_c_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	  }
	  //#endif

	  QLA_c_peq_r_times_c(z, clov_mkappa, QLA_elem_D(*dsl, ic, is));
	  QLA_c_eq_c_plus_c(QLA_elem_D(clov_out[x], ic, is), z, QLA_elem_D(clov_in[x],ic,is));
	}
#else  // unrolled version
#define clov_diag(i) clov_clov[xb+i]
#define clov_offd(i) cmplx(clov_clov[xb+6+i])
#define src(i) QLA_elem_D(clov_in[x],i/2,2*b+i%2)
#define dsrc(i) QLA_elem_D(*dsl,i/2,2*b+i%2)
#define dest(i) QLA_elem_D(clov_out[x],i/2,2*b+i%2)
	{
	  QLA_Complex z;

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_r_times_c(z, clov_diag(0), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(0), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(1), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(2), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(3), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(4), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(0));
	  QLA_c_eq_c_plus_c(dest(0), z, src(0));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(0), src(0));
	  QLA_c_peq_r_times_c(z, clov_diag(1), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(5), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(6), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(7), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(8), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(1));
	  QLA_c_eq_c_plus_c(dest(1), z, src(1));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(1), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(5), src(1));
	  QLA_c_peq_r_times_c(z, clov_diag(2), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(9), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(10), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(11), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(2));
	  QLA_c_eq_c_plus_c(dest(2), z, src(2));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(2), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(6), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(9), src(2));
	  QLA_c_peq_r_times_c(z, clov_diag(3), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(12), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(13), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(3));
	  QLA_c_eq_c_plus_c(dest(3), z, src(3));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(3), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(7), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(10), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(12), src(3));
	  QLA_c_peq_r_times_c(z, clov_diag(4), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(14), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(4));
	  QLA_c_eq_c_plus_c(dest(4), z, src(4));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(4), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(8), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(11), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(13), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(14), src(4));
	  QLA_c_peq_r_times_c(z, clov_diag(5), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(5));
	  QLA_c_eq_c_plus_c(dest(5), z, src(5));
	}
#endif
      }
    }
  }
  QDP_reset_D(in);
  QDP_reset_D(out);
  if(dsl!=out) QDP_reset_D(dsl);
#endif
}

static void
clov(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out, QLA_Real *mkappa,
     QDP_DiracFermion *dsl, QDP_DiracFermion *in, QDP_Subset subset)
{
  if(flw->clov==NULL) {
    QDP_D_eq_r_times_D_plus_D(out, mkappa, dsl, in, subset);
  } else {
    apply_clov(flw->clov, out, mkappa, dsl, in, subset);
  }
}

static void
clovinv(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out, QLA_Real *mkappa,
	QDP_DiracFermion *dsl, QDP_DiracFermion *in, QDP_Subset subset)
{
  if(flw->clov==NULL) {
    QDP_D_eq_r_times_D_plus_D(out, mkappa, dsl, in, subset);
  } else {
    apply_clov(flw->clovinv, out, mkappa, dsl, in, subset);
  }
}


/************ dslash *************/

/* Special dslash for use by congrad.  Uses restart_gather() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather, otherwise use restart_gather.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
static void
wilson_dslash0(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp)
{
  int mu;
  QDP_DiracFermion *vsrc[4];
  QDP_DiracFermion *vdest[4];
  QDP_ShiftDir fwd[4], bck[4];
  int dir[4], sgn[4], msgn[4];
  QDP_Subset othersubset;

  sign = -sign;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vdest[mu] = dest;
    fwd[mu] = QDP_forward;
    bck[mu] = QDP_backward;
    dir[mu] = mu;
    sgn[mu] = sign;
    msgn[mu] = -sign;
  }
  //sgn[1] = -sign;
  //msgn[1] = sign;

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  printf0("dslash0\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_D_veq_sD\n");
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, QDP_neighbor+mu, fwd+mu, subset,
		   QOP_wilson_nsvec);
      printf0("end QDP_D_veq_sD\n");
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_D\n");
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, QDP_neighbor+mu, fwd+mu,
		   subset, QOP_wilson_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_wilson_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
#if 0
      QDP_H_veq_spproj_Ma_times_D(hf+mu, fwdlinks+mu, vsrc+mu,
                               dir+mu, msgn+mu, othersubset, QOP_wilson_nsvec);
      QDP_D_veq_sprecon_H(dtemp[ntmp]+8+mu, hf+mu,
                               dir+mu, msgn+mu, othersubset, QOP_wilson_nsvec);
#endif
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_Ma_times_D\n");
      QDP_H_veq_spproj_Ma_times_D(htemp[ntmp]+12+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+12+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  printf0("dslash0 - fwd\n");
  QDP_D_eq_zero(dest, subset);

  if(shiftd_style(QOP_wilson_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->links+mu, dtemp[ntmp]+mu,
				  dir+mu, sgn+mu, subset, QOP_wilson_nvec);
#if 0
      QDP_H_veq_spproj_M_times_D(hf+mu, flw->links+mu, dtemp[ntmp]+mu,
				 dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      QDP_D_vpeq_sprecon_H(vdest+mu, hf+mu,
			   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
#endif
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->links+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_D(vdest+mu, dtemp[ntmp]+4+mu, subset, QOP_wilson_nvec);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_H(vdest+mu, htemp[ntmp]+4+mu, dir+mu, msgn+mu, subset,
			   QOP_wilson_nvec);
    }
  }

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
} /* end of dslash_special_qdp() */

/* Special dslash for use by congrad.  Uses restart_gather() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather, otherwise use restart_gather.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
static void
wilson_dslash1(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp)
{
  int mu;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];
  QDP_Subset othersubset=QDP_all;

  if(!shiftd_style(QOP_wilson_style)) {
    if(subset==QDP_even) othersubset = QDP_odd;
    else if(subset==QDP_odd) othersubset = QDP_even;
    else othersubset = QDP_all;
  }

  sign = -sign;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vsrc[mu+4] = src;
    vdest[mu] = dest;
    vdest[mu+4] = dest;
    dir[2*mu] = mu;
    dir[2*mu+1] = mu;
    sgn[2*mu] = sign;
    sgn[2*mu+1] = -sign;
    sh[2*mu] = QDP_neighbor[mu];
    sh[2*mu+1] = QDP_neighbor[mu];
    sd[2*mu] = QDP_forward;
    sd[2*mu+1] = QDP_backward;
  }
  //sgn[2] = -sign;
  //sgn[3] = sign;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  //printf0("ds1 1\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, sh+mu, sd+mu, subset,
		   QOP_wilson_nsvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, sh+mu, sd+mu, subset,
		   QOP_wilson_nsvec);
    }
  }
  //printf0("ds1 2\n");

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add to dest */

  QDP_D_eq_zero(dest, subset);
  //printf0("ds1 3\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->dbllinks+mu, dtemp[ntmp]+mu,
			          dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->dbllinks+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  }
  //printf0("ds1 4\n");

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
} /* end of dslash_special_qdp() */

static void
wilson_mdslash1(QOP_FermionLinksWilson *flw,
		QDP_DiracFermion *out, QDP_DiracFermion *in,
		int sign, QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
#ifdef LU

  //printf0("here3\n");
  dslash_special_qdp(flw, tt2, in, sign, QDP_odd, 2);
  dslash_special_qdp(flw, out, tt2, sign, QDP_even, 3);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);

#else

  //printf0("here6\n");
  dslash_special_qdp(flw, out, in, sign, QDP_all, 1);
  //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
  clov(flw, out, &mkappa, out, in, QDP_all);

#endif
}

static void
wilson_mdslash2(QOP_FermionLinksWilson *flw,
		QDP_DiracFermion *out, QDP_DiracFermion *in,
		QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
#ifdef LU

  if(QOP_wilson_cgtype==1) {
    if(in==flw->cgp) {
      dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 0);
      dslash_special_qdp(flw, out, tt1, 1, QDP_even, 1);
      QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);
    } else {
      dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 2);
      dslash_special_qdp(flw, out, tt1, 1, QDP_even, 1);
      QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);
    }
  } else {
    dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 0);
    dslash_special_qdp(flw, ttt, tt1, 1, QDP_even, 1);
    QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_even);

    dslash_special_qdp(flw, tt2, ttt, -1, QDP_odd, 2);
    dslash_special_qdp(flw, out, tt2, -1, QDP_even, 3);
    QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_even);
  }

#else

  if(QOP_wilson_cgtype==1) {
    if(in==flw->cgp) {
      dslash_special_qdp(flw, out, in, 1, QDP_all, 0);
      //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
      clov(flw, out, &mkappa, out, in, QDP_all);
    } else {
      dslash_special_qdp(flw, out, in, 1, QDP_all, 1);
      //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
      clov(flw, out, &mkappa, out, in, QDP_all);
    }
  } else {
    //printf0("here6\n");
    dslash_special_qdp(flw, ttt, in, 1, QDP_all, 0);
    //QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_all);
    clov(flw, ttt, &mkappa, ttt, in, QDP_all);
    //printf0("here7\n");
    dslash_special_qdp(flw, out, ttt, -1, QDP_all, 1);
    //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_all);
    clov(flw, out, &mkappa, out, ttt, QDP_all);
  }

#endif
}
