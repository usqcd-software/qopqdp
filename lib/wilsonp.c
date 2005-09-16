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

#define LU

#ifdef LU
#define MYSUBSET QDP_even
#else
#define MYSUBSET QDP_all
#endif

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

static int old_style=0;
static int old_nsvec=0;
static int old_nvec=0;

static int congrad_setup = 0;
static int dblstored = 0;
static QDP_ColorMatrix *fwdlinks[4];
static QDP_ColorMatrix *bcklinks[4], *dbllinks[8];
static QDP_HalfFermion *htemp[4][20];
static QDP_DiracFermion *dtemp[4][12];
static QDP_DiracFermion *psi, *chi, *cgp, *cgr, *mp, *ttt, *tt1, *tt2, *t1, *t2, *t3;

static void
free_temps(void)
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
double_store(void)
{
  if( dblstore_style(QOP_wilson_style) && (!dblstored) ) {
    int i;
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<4; i++) {
      QDP_M_eq_sM(m, fwdlinks[i], QDP_neighbor[i], QDP_backward, QDP_all);
      QDP_M_eq_Ma(bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    dblstored = 1;
  }
}

static void
reset_temps(void)
{
  int i, j;

  if(QOP_wilson_style!=old_style) {
    if(!dblstore_style(QOP_wilson_style)) {
      if(congrad_setup) {
        for(i=0; i<4; i++) {
          QDP_destroy_M(bcklinks[i]);
        }
      }
    } else {
      for(i=0; i<4; i++) {
        bcklinks[i] = QDP_create_M();
      }
      for(i=0; i<4; i++) {
        dbllinks[2*i] = fwdlinks[i];
        dbllinks[2*i+1] = bcklinks[i];
      }
    }
    dblstored = 0;
  }
  double_store();

  free_temps();

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

QOP_status_t
PREC(wilson_invert_load_links_qdp)(QDP_ColorMatrix *lnk[])
{
  int i;

  dblstored = 0;
  for(i=0; i<4; i++) {
    fwdlinks[i] = lnk[i];
  }

  if(!congrad_setup) {
    reset_temps();
  }

  if( (QOP_wilson_style != old_style) ||
      (QOP_wilson_nsvec != old_nsvec) ||
      (QOP_wilson_nvec != old_nvec) ) {
    reset_temps();
    old_style = QOP_wilson_style;
    old_nsvec = QOP_wilson_nsvec;
    old_nvec = QOP_wilson_nvec;
  }

  double_store();

  return QOP_SUCCESS;
}

QOP_status_t
PREC(wilson_invert_unload_links)(void)
{
  dblstored = 0;
  return QOP_SUCCESS;
}

static void wilson_mdslash2(QDP_DiracFermion *out, QDP_DiracFermion *in,
			    QDP_Subset subset, QDP_Subset othersubset,
			    QLA_Real mkappa);

int
PREC(wilson_inv_qdp)(QOP_invert_arg *inv_arg,
		     QDP_DiracFermion *out, QDP_DiracFermion *in)
{
  double source_norm;
  double rsqstop;
  QLA_Real a, b;
  double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
  QLA_Real mkappa;
  QLA_Real sum;
  double dtime;
  int iteration;	/* counter for iterations */
#ifdef LU
  mkappa = -inv_arg->mass*inv_arg->mass;
#else
  mkappa = -inv_arg->mass;
#endif

  if( (QOP_wilson_style != old_style) ||
      (QOP_wilson_nsvec != old_nsvec) ||
      (QOP_wilson_nvec != old_nvec) ) {
    reset_temps();
    old_style = QOP_wilson_style;
    old_nsvec = QOP_wilson_nsvec;
    old_nvec = QOP_wilson_nvec;
  }

  //printf0("here1\n");
  QDP_D_eq_D(chi, in, QDP_all);
  QDP_D_eq_D(psi, out, QDP_all);
  //printf0("here2\n");

  dtime = -dclock();

  iteration=0;
 start:
  /* mp <-  M_adjoint*M*psi
     r,p <- chi - mp
     rsq = |r|^2
     source_norm = |chi|^2
  */
  rsq = source_norm = 0.0;

  QDP_D_eq_D(cgp, psi, QDP_even);

#ifdef LU

  wilson_mdslash2(mp, cgp, QDP_even, QDP_odd, mkappa);

  QDP_D_eq_D_minus_D(cgr, chi, mp, QDP_even);
  QDP_D_eq_D(cgp, cgr, QDP_even);

  //printf0("here5\n");
  QDP_r_eq_norm2_D(&sum, chi, QDP_even);
  source_norm = sum;
  QDP_r_eq_norm2_D(&sum, cgr, QDP_even);
  rsq = sum;

#else

  wilson_mdslash2(mp, cgp, QDP_all, QDP_all, mkappa);

  //printf0("here8\n");
  QDP_D_eq_D_minus_D(cgr, chi, mp, QDP_all);
  QDP_D_eq_D(cgp, cgr, QDP_all);

  //printf0("here9\n");
  QDP_r_eq_norm2_D(&sum, chi, QDP_all);
  source_norm = sum;
  QDP_r_eq_norm2_D(&sum, cgr, QDP_all);
  rsq = sum;

#endif
  //printf0("here10\n");

  iteration++ ;	/* iteration counts number of multiplications
		   by M_adjoint*M */
  /**if(this_node==0)printf("congrad2: source_norm = %e\n",source_norm);
     if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
     iteration,(double)rsq,(double)pkp,(double)a );**/
  rsqstop = inv_arg->rsqmin * source_norm;
  if( rsq <= rsqstop ){
    inv_arg->final_rsq = (float)rsq;
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
#ifdef LU
    wilson_mdslash2(mp, cgp, QDP_even, QDP_odd, mkappa);

    QDP_r_eq_re_D_dot_D(&sum, cgp, mp, QDP_even);
    pkp = sum;
#else
    wilson_mdslash2(mp, cgp, QDP_all, QDP_all, mkappa);

    QDP_r_eq_re_D_dot_D(&sum, cgp, mp, QDP_all);
    pkp = sum;
#endif
    iteration++;

    a = rsq / pkp;
    QDP_D_peq_r_times_D(psi, &a, cgp, MYSUBSET);
    QDP_D_meq_r_times_D(cgr, &a, mp, MYSUBSET);
    QDP_r_eq_norm2_D(&sum, cgr, MYSUBSET);
    rsq = sum;

    /**if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
       iteration,(double)rsq,(double)pkp,(double)a );**/
    if( rsq <= rsqstop ){
      inv_arg->final_rsq = (float)rsq;
      dtime += dclock();
#if 0
      printf0("CONGRAD2: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	      dtime, rsq, iteration,
	      (double)6480*iteration*QDP_sites_on_node/(2*dtime*1e6));
      //(double)5616*iteration*even_sites_on_node/(dtime*1e6));
#endif
      QDP_D_eq_D(out, psi, QDP_all);
      inv_arg->final_flop = 0.5*6480.0*iteration*QDP_sites_on_node;
      inv_arg->final_sec = dtime;
      inv_arg->final_iter = iteration;
      return (iteration);
    }

    b = rsq / oldrsq;
    QDP_D_eq_r_times_D_plus_D(cgp, &b, cgp, cgr, MYSUBSET);

  } while( iteration%inv_arg->restart != 0);

  QDP_D_eq_D(out, psi, QDP_all);

  if( iteration < inv_arg->max_iter ) goto start;
  dtime += dclock();
  inv_arg->final_rsq = (float)rsq;
  inv_arg->final_flop = 0.5*6480.0*iteration*QDP_sites_on_node;
  inv_arg->final_sec = dtime;
  inv_arg->final_iter = iteration;
  return(iteration);
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
wilson_dslash0(QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp)
{
  int mu;
  QDP_DiracFermion *vsrc[4];
  QDP_DiracFermion *vdest[4];
  QDP_ShiftDir fwd[4], bck[4];
  int dir[4], sgn[4], msgn[4];
  QDP_Subset othersubset;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vdest[mu] = dest;
    fwd[mu] = QDP_forward;
    bck[mu] = QDP_backward;
    dir[mu] = mu;
    sgn[mu] = sign;
    msgn[mu] = -sign;
  }

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, QDP_neighbor+mu, fwd+mu, subset,
		   QOP_wilson_nsvec);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, QDP_neighbor+mu, fwd+mu,
		   subset, QOP_wilson_nsvec);
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  if(shiftd_style(QOP_wilson_style)) {
    QDP_HalfFermion *hf1, *hf2;
    hf1 = QDP_create_H();
    hf2 = QDP_create_H();
    for(mu=0; mu<4; mu++) {
      QDP_H_eq_spproj_D(hf1, src, mu, -sign, othersubset);
      QDP_H_eq_Ma_times_H(hf2, fwdlinks[mu], hf1, othersubset);
      QDP_D_eq_sprecon_H(dtemp[ntmp][8+mu], hf2, mu, -sign, othersubset);
      QDP_D_eq_sD(dtemp[ntmp][4+mu], dtemp[ntmp][8+mu], QDP_neighbor[mu],
		  QDP_backward, subset);
    }
    QDP_destroy_H(hf1);
    QDP_destroy_H(hf2);
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+12+mu, vsrc+mu, dir+mu, msgn+mu,
			othersubset, QOP_wilson_nsvec);
      QDP_H_veq_Ma_times_H(htemp[ntmp]+16+mu, fwdlinks+mu, htemp[ntmp]+12+mu,
			   othersubset, QOP_wilson_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+16+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
    }
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  QDP_D_eq_zero(dest, subset);

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_wilsonspin_M_times_D(vdest+mu, fwdlinks+mu, dtemp[ntmp]+mu,
				      dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, fwdlinks+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_D(vdest+mu, dtemp[ntmp]+4+mu, subset, QOP_wilson_nvec);
      QDP_discard_D(dtemp[ntmp][4+mu]);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_H(vdest+mu, htemp[ntmp]+4+mu, dir+mu, msgn+mu, subset,
			   QOP_wilson_nvec);
      QDP_discard_H(htemp[ntmp][4+mu]);
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
wilson_dslash1(QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp)
{
  int mu;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];
  QDP_Subset othersubset;

  if(!shiftd_style(QOP_wilson_style)) {
    if(subset==QDP_even) othersubset = QDP_odd;
    else if(subset==QDP_odd) othersubset = QDP_even;
    else othersubset = QDP_all;
  }

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
      QDP_D_vpeq_wilsonspin_M_times_D(vdest+mu, dbllinks+mu, dtemp[ntmp]+mu,
				      dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, dbllinks+mu, htemp[ntmp]+mu,
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

static void (*dslash_special_qdp)(QDP_DiracFermion *dest,
				  QDP_DiracFermion *src,
				  int sign, QDP_Subset subset, int ntmp);

static void
wilson_mdslash2(QDP_DiracFermion *out, QDP_DiracFermion *in,
		QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
  if(dblstore_style(QOP_wilson_style)) {
    dslash_special_qdp = wilson_dslash1;
  } else {
    dslash_special_qdp = wilson_dslash0;
  }

#ifdef LU

  //printf0("here3\n");
  dslash_special_qdp(tt1, cgp, 1, QDP_odd, 0);
  dslash_special_qdp(ttt, tt1, 1, QDP_even, 1);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_even);
  //printf0("here4\n");
  dslash_special_qdp(tt2, ttt, -1, QDP_odd, 2);
  dslash_special_qdp(mp, tt2, -1, QDP_even, 3);
  QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_even);

#else

  //printf0("here6\n");
  dslash_special_qdp(ttt, cgp, 1, QDP_all, 0);
  QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, cgp, QDP_all);
  //printf0("here7\n");
  dslash_special_qdp(mp, ttt, -1, QDP_all, 1);
  QDP_D_eq_r_times_D_plus_D(mp, &mkappa, mp, ttt, QDP_all);

#endif
}
