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

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "phi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

//#define QDP_PROFILE
#include <qop_internal.h>

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_wilson_inited;
extern int QOP_wilson_style;
extern int QOP_wilson_nsvec;
extern int QOP_wilson_nvec;

static int old_style=-1;
static int old_nsvec=-1;
static int old_nvec=-1;

static int congrad_setup = 0;
static QDP_HalfFermion *htemp[4][20];
static QDP_DiracFermion *dtemp[4][12];
static QDP_DiracFermion *psi, *chi, *cgp, *cgr, *mp, *ttt, *tt1, *tt2, *t1, *t2, *t3;

static void
free_temps(QOP_FermionLinksWilson *flw)
{
  if(congrad_setup) {
    int i, j;

    QDP_destroy_D(psi);
    QDP_destroy_D(chi);
    QDP_destroy_D(cgp);
    QDP_destroy_D(cgr);
    QDP_destroy_D(mp);
    QDP_destroy_D(ttt);
    QDP_destroy_D(tt1);
    QDP_destroy_D(tt2);
    QDP_destroy_D(t1);
    QDP_destroy_D(t2);
    QDP_destroy_D(t3);

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
        for(i=0; i<4; i++) {
          QDP_destroy_M(flw->bcklinks[i]);
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

  psi = QDP_create_D();
  chi = QDP_create_D();
  cgp = QDP_create_D();
  cgr = QDP_create_D();
  mp = QDP_create_D();
  ttt = QDP_create_D();
  tt1 = QDP_create_D();
  tt2 = QDP_create_D();
  t1 = QDP_create_D();
  t2 = QDP_create_D();
  t3 = QDP_create_D();

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

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(REAL *links[], QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_GaugeField *gauge)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

void
QOP_wilson_extract_L_to_raw(REAL *links[], QOP_FermionLinksWilson *src,
			       QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
}

void
QOP_wilson_destroy_L(QOP_FermionLinksWilson *field)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

REAL **
QOP_wilson_convert_L_to_raw(QOP_FermionLinksWilson *src,
			    QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_G(QOP_GaugeField *gauge)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

QOP_GaugeField *
QOP_wilson_convert_L_to_G(QOP_FermionLinksWilson *links)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_FermionLinksWilson *flw;
  int i;

  QOP_malloc(flw, QOPPC(FermionLinksWilson), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);

  flw->dblstored = 0;
  for(i=0; i<4; i++) {
    flw->links[i] = QDP_create_M();
    QDP_M_eq_M(flw->links[i], links[i], QDP_all);
  }

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

  return flw;
}

void
QOP_wilson_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QOP_FermionLinksWilson *src)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_qdp(QDP_ColorMatrix *links[])
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

QDP_ColorMatrix **
QOP_wilson_convert_L_to_qdp(QOP_FermionLinksWilson *src)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}

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

QOP_status_t
QOP_wilson_invert_multi(QOP_FermionLinksWilson *links,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t **res_arg[],
			REAL *kappas[], int nkappa[],
			QOP_DiracFermion **out_pt[],
			QOP_DiracFermion *in_pt[],
			int nsrc)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return QOP_SUCCESS;
}

QOP_status_t
QOPPC(wilson_invert)(QOP_FermionLinksWilson *flw,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     REAL kappa,
		     QOP_DiracFermion *out_pt,
		     QOP_DiracFermion *in_pt)
{
  double source_norm;
  double rsqstop;
  QLA_Real a, b;
  double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
  QLA_Real mkappa;
  QLA_Real sum;
  double dtime;
  double nflop = 0.5 * 6096;
  int iteration;	/* counter for iterations */
  QDP_DiracFermion *in, *out;
  QDP_Subset subset, osubset;

  in = in_pt->df;
  out = out_pt->df;
#ifdef LU
  mkappa = -kappa*kappa;
  subset = QDP_even;
  osubset = QDP_odd;
#else
  mkappa = -kappa;
  subset = QDP_all;
  osubset = QDP_all;
#endif

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

  //printf0("here1\n");
  {
#ifdef LU
    QDP_D_eq_D(tt1, in, osubset);
    dslash_special_qdp(flw, ttt, tt1, 1, subset, 1);
    QDP_D_eq_r_times_D_plus_D(cgp, &kappa, ttt, in, subset);
    wilson_mdslash1(flw, chi, cgp, -1, subset, osubset, mkappa);
    QDP_D_eq_D(psi, out, subset);
#else
    QLA_Real kappa2 = 2.0*kappa;
    QDP_D_eq_r_times_D(cgp, &kappa2, in, subset);
    wilson_mdslash1(flw, chi, cgp, -1, subset, osubset, mkappa);
    QDP_D_eq_D(psi, out, subset);
#endif
  }
  //printf0("here2\n");

  dtime = -QOP_time();

  iteration=0;
 start:
  /* mp <-  M_adjoint*M*psi
     r,p <- chi - mp
     rsq = |r|^2
     source_norm = |chi|^2
  */
  rsq = source_norm = 0.0;

  QDP_D_eq_D(cgp, psi, subset);

  printf0("test0\n");
  wilson_mdslash2(flw, mp, cgp, subset, osubset, mkappa);
  printf0("test1\n");

  QDP_D_eq_D_minus_D(cgr, chi, mp, subset);
  QDP_D_eq_D(cgp, cgr, subset);

  QDP_r_eq_norm2_D(&sum, chi, subset);
  source_norm = sum;
  QDP_r_eq_norm2_D(&sum, cgr, subset);
  rsq = sum;

  iteration++ ;	/* iteration counts number of multiplications
		   by M_adjoint*M */
  /**if(this_node==0)printf("congrad2: source_norm = %e\n",source_norm);
     if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
     iteration,(double)rsq,(double)pkp,(double)a );**/
  rsqstop = res_arg->rsqmin * source_norm;
  if( rsq <= rsqstop ){
    res_arg->final_rsq = (float)rsq;
    return (iteration);
  }

  /* main loop - do until convergence or time to restart */
  /* 
     oldrsq <- rsq
     mp <- M_adjoint*M*p
     pkp <- p.M_adjoint*M.p
     a <- rsq/pkp
     psi <- psi + a*p
     r <- r - a*mp
     rsq <- |r|^2
     b <- rsq/oldrsq
     p <- r + b*p
  */
  do {
    //if(iteration%10==0) printf0("%i %g\n", iteration, rsq);
    oldrsq = rsq;

    wilson_mdslash2(flw, mp, cgp, subset, osubset, mkappa);

    QDP_r_eq_re_D_dot_D(&sum, cgp, mp, subset);
    pkp = sum;
    iteration++;

    a = rsq / pkp;
    QDP_D_peq_r_times_D(psi, &a, cgp, subset);
    QDP_D_meq_r_times_D(cgr, &a, mp, subset);
    QDP_r_eq_norm2_D(&sum, cgr, subset);
    rsq = sum;

    /**if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
       iteration,(double)rsq,(double)pkp,(double)a );**/
    if( rsq <= rsqstop ){
      res_arg->final_rsq = (float)rsq;
      dtime += QOP_time();
#if 0
      printf0("CONGRAD2: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	      dtime, rsq, iteration,
	      (double)6480*iteration*QDP_sites_on_node/(2*dtime*1e6));
      //(double)5616*iteration*even_sites_on_node/(dtime*1e6));
#endif
#ifdef LU
      {
	QLA_Real kappa2 = 2.0*kappa;
	QDP_D_eq_D(cgp, psi, subset);
	dslash_special_qdp(flw, tt1, cgp, 1, osubset, 0);
	QDP_D_eq_r_times_D_plus_D(psi, &kappa, tt1, in, osubset);
	QDP_D_eq_r_times_D(out, &kappa2, psi, QDP_all);
      }
#else
      QDP_D_eq_D(out, psi, subset);
#endif

      inv_arg->final_flop = nflop*iteration*QDP_sites_on_node;
      inv_arg->final_sec = dtime;
      res_arg->final_iter = iteration;
      return (iteration);
    }

    b = rsq / oldrsq;
    QDP_D_eq_r_times_D_plus_D(cgp, &b, cgp, cgr, subset);

  } while( iteration%inv_arg->restart != 0);

#ifdef LU
  {
    QLA_Real kappa2 = 2.0*kappa;
    QDP_D_eq_D(cgp, psi, subset);
    dslash_special_qdp(flw, tt1, cgp, 1, osubset, 0);
    QDP_D_eq_r_times_D_plus_D(psi, &kappa, tt1, in, osubset);
    QDP_D_eq_r_times_D(out, &kappa2, psi, QDP_all);
  }
#else
  QDP_D_eq_D(out, psi, subset);
#endif

  if( iteration < inv_arg->max_iter ) goto start;
  dtime += QOP_time();
  res_arg->final_rsq = (float)rsq;
  inv_arg->final_flop = nflop*iteration*QDP_sites_on_node;
  inv_arg->final_sec = dtime;
  res_arg->final_iter = iteration;

  return QOP_SUCCESS;
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
  sgn[1] = -sign;
  msgn[1] = sign;

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
  sgn[2] = -sign;
  sgn[3] = sign;

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
  dslash_special_qdp(flw, tt1, in, sign, QDP_odd, 0);
  dslash_special_qdp(flw, out, tt1, sign, QDP_even, 1);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);

#else

  //printf0("here6\n");
  dslash_special_qdp(flw, out, in, sign, QDP_all, 0);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);

#endif
}

static void
wilson_mdslash2(QOP_FermionLinksWilson *flw,
		QDP_DiracFermion *out, QDP_DiracFermion *in,
		QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
#ifdef LU

  //printf0("here3\n");
  dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 0);
  dslash_special_qdp(flw, ttt, tt1, 1, QDP_even, 1);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_even);
  //printf0("here4\n");
  dslash_special_qdp(flw, tt2, ttt, -1, QDP_odd, 2);
  dslash_special_qdp(flw, out, tt2, -1, QDP_even, 3);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_even);

#else

  //printf0("here6\n");
  dslash_special_qdp(flw, ttt, in, 1, QDP_all, 0);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_all);
  //printf0("here7\n");
  dslash_special_qdp(flw, out, ttt, -1, QDP_all, 1);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_all);

#endif
}
