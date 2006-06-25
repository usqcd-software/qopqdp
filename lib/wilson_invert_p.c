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

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(REAL *links[], REAL *clov[],
			     QOP_evenodd_t evenodd)
{
  QOP_FermionLinksWilson *flw;
  QOP_GaugeField *gf;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_wilson_convert_L_from_qdp(gf->links, NULL);

  flw->raw = NULL;
  flw->qopgf = gf;

  return flw;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			   QOP_GaugeField *gauge)
{
  QOP_error("QOP_wilson_create_L_from_G unimplemented.");
  return NULL;
}

void
QOP_wilson_extract_L_to_raw(REAL *links[], REAL *clov[],
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
  free(flw->bcklinks);
  free(flw->dbllinks);
  free(flw);
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], REAL *clov[],
			      QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_convert_L_from_raw unimplemented");
  return NULL;
}

void
QOP_wilson_convert_L_to_raw(REAL ***links, REAL ***clov,
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
			     QDP_DiracPropagator *clov[])
{
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_wilson_convert_L_from_qdp(newlinks, NULL);

  return flw;
}

void
QOP_wilson_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QDP_DiracPropagator *clov[],
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_extract_L_to_qdp unimplemented");
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_qdp(QDP_ColorMatrix *links[],
			      QDP_DiracPropagator *clov[])
{
  QOP_FermionLinksWilson *flw;
  int i;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

  QOP_malloc(flw, QOPPC(FermionLinksWilson), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);

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

  flw->raw = NULL;
  return flw;
}

void
QOP_wilson_convert_L_to_qdp(QDP_ColorMatrix ***links,
			    QDP_DiracPropagator ***clov,
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_convert_L_to_qdp unimplemented");
}


/* inverter stuff */

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
  double dtime;
  double nflop;
  QLA_Real mkappa;
  QDP_DiracFermion *qdpin;
  QDP_Subset subset, osubset;

  /* cg has 5 * 48 = 240 flops/site/it */
#ifdef LU
  /* MdagM -> 2*(2*(144+168*7)+48) = 5376 flops/site */
  nflop = 0.5 * 5616;  /* half for half of the sites */
  mkappa = -kappa*kappa;
  subset = QDP_even;
  osubset = QDP_odd;
#else
  /* MdagM -> 2*((144+168*7)+48) = 2736 flops/site */
  nflop = 2976;
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
    wilson_mdslash1(flw, qdpin, ttt, -1, subset, osubset, mkappa);
#else
    QDP_D_eq_r_times_D(ttt, &kappa2, in->df, subset);
    wilson_mdslash1(flw, qdpin, ttt, -1, subset, osubset, mkappa);
#endif
  }

  dtime = -QOP_time();

  QOPPC(invert_cg_D)(QOPPC(wilson_dslash2), inv_arg, res_arg,
		     out->df, qdpin, flw->cgp, subset);

  dtime += QOP_time();

#ifdef LU
      {
	QDP_D_eq_D(ttt, out->df, subset);
	dslash_special_qdp(flw, tt2, ttt, 1, osubset, 2);
	QDP_D_eq_r_times_D_plus_D(ttt, &kappa, tt2, qdpin, osubset);
	QDP_D_eq_D(out->df, ttt, QDP_all);
      }
#endif
  QDP_destroy_D(qdpin);

  //res_arg->final_rsq = rsq;
  //res_arg->final_iter = iteration;
  //inv_arg->final_iter = iteration;
  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
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
