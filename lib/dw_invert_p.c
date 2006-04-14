#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

//#define LU

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_dw_inited;
extern int QOP_dw_style;
extern int QOP_dw_nsvec;
extern int QOP_dw_nvec;

static int old_style=-1;
static int old_nsvec=-1;
static int old_nvec=-1;

static int congrad_setup = 0;
static QDP_HalfFermion *htemp[4][20];
static QDP_DiracFermion *dtemp[4][12];
static QDP_DiracFermion *tt1, *tt2, **ttv;

static void
free_temps(QOP_FermionLinksDW *flw)
{
  if(congrad_setup) {
    int i, j;

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
double_store(QOP_FermionLinksDW *flw)
{
  if( dblstore_style(QOP_dw_style) && (!flw->dblstored) ) {
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
reset_temps(QOP_FermionLinksDW *flw)
{
  int i, j;

  if(QOP_dw_style!=old_style) {
    if(!dblstore_style(QOP_dw_style)) {
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

  tt1 = QDP_create_D();
  tt2 = QDP_create_D();

  if(shiftd_style(QOP_dw_style)) {
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

QOP_FermionLinksDW *
QOP_dw_create_L_from_raw(REAL *links[], QOP_evenodd_t evenodd)
{
  QOP_FermionLinksDW *flw;
  QOP_GaugeField *gf;

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_dw_convert_L_from_G(gf);

  flw->raw = NULL;
  flw->qopgf = gf;

  return flw;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_G(QOP_GaugeField *gauge)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return NULL;
}

void
QOP_dw_extract_L_to_raw(REAL *links[], QOP_FermionLinksDW *src,
			       QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
}

void
QOP_dw_destroy_L(QOP_FermionLinksDW *flw)
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
  free(flw->bcklinks);
  free(flw->dbllinks);
  free(flw);
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_raw(REAL *links[], QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return NULL;
}

REAL **
QOP_dw_convert_L_to_raw(QOP_FermionLinksDW *src,
			    QOP_evenodd_t evenodd)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return NULL;
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_G(QOP_GaugeField *gauge)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return NULL;
}

QOP_GaugeField *
QOP_dw_convert_L_to_G(QOP_FermionLinksDW *links)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return NULL;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_FermionLinksDW *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_dw_convert_L_from_qdp(newlinks);

  return flw;
}

void
QOP_dw_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QOP_FermionLinksDW *src)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_FermionLinksDW *flw;
  int i;

  QOP_malloc(flw, QOPPC(FermionLinksDW), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);

  flw->dblstored = 0;
  for(i=0; i<4; i++) {
    flw->links[i] = links[i];
  }

  if( (!congrad_setup) ||
      (QOP_dw_style != old_style) ||
      (QOP_dw_nsvec != old_nsvec) ||
      (QOP_dw_nvec != old_nvec) ) {
    reset_temps(flw);
    old_style = QOP_dw_style;
    old_nsvec = QOP_dw_nsvec;
    old_nvec = QOP_dw_nvec;
  }

  double_store(flw);

  flw->raw = NULL;
  return flw;
}

QDP_ColorMatrix **
QOP_dw_convert_L_to_qdp(QOP_FermionLinksDW *src)
{
  fprintf(stderr, "unimpleented\n");
  QDP_abort();
  return NULL;
}


/******************/
/* inverter stuff */
/******************/

static void
(*wilson_dslash)(QOP_FermionLinksDW *flw,
		 QDP_DiracFermion *dest, QDP_DiracFermion *src,
		 int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash0(QOP_FermionLinksDW *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash1(QOP_FermionLinksDW *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QDP_Subset subset, int ntmp);


static void dw_dslash1(QOP_FermionLinksDW *flw,
		       QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		       int sign, QDP_Subset subset, QDP_Subset othersubset,
		       QLA_Real m0, QLA_Real M, int Ls);

static void dw_dslash2(QOP_FermionLinksDW *flw,
		       QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		       QDP_Subset subset, QDP_Subset othersubset,
		       QLA_Real m0, QLA_Real M, int Ls);

QOP_status_t
QOP_dw_invert_multi(QOP_FermionLinksDW *links,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t **res_arg[],
		    REAL *m0[],
		    REAL *M[],
		    int nmass[],
		    QOP_DiracFermion ***out_pt[],
		    QOP_DiracFermion **in_pt[],
		    int Ls,
		    int nsrc)
{
  fprintf(stderr, "unimplemented\n");
  QDP_abort();
  return QOP_SUCCESS;
}

static QLA_Real gl_m0, gl_M;
static int gl_Ls;
static QDP_Subset gl_osubset;
static QOP_FermionLinksDW *gl_flw;

void
QOPPC(dw_dslash2)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		  QDP_Subset subset)
{
  dw_dslash2(gl_flw, out, in, subset, gl_osubset, gl_m0, gl_M, gl_Ls);
}

QOP_status_t
QOPPC(dw_invert)(QOP_FermionLinksDW *flw,
		 QOP_invert_arg_t *inv_arg,
		 QOP_resid_arg_t *res_arg,
		 REAL m0,
		 REAL M,
		 QOP_DiracFermion *out[],
		 QOP_DiracFermion *in[],
		 int Ls)
{
  double dtime;
  double nflop;
  QLA_Real kappa, mkappa;
  QDP_DiracFermion *qdptmp[Ls], *qdpin[Ls], *qdpout[Ls];
  QDP_Subset subset, osubset;
  int i;

  kappa = m0;
  /* cg has 5 * 48 *Ls = 240*Ls flops/site/it */
#ifdef LU
  /* ???? MdagM -> 2*(2*(144+168*7)+48) = 5376 flops/site */
  nflop = 0.5 * 6096;
  mkappa = -kappa*kappa;
  subset = QDP_even;
  osubset = QDP_odd;
#else
  /* MdagM -> 2*(48 + ((144+168*7)+96)*Ls) = 96 + 2832*Ls flops/site */
  nflop = 96 + 3072*Ls;
  mkappa = -kappa;
  subset = QDP_all;
  osubset = QDP_all;
#endif

  gl_osubset = osubset;
  gl_m0 = m0;
  gl_M = M;
  gl_Ls = Ls;
  gl_flw = flw;

  if( (QOP_dw_style != old_style) ||
      (QOP_dw_nsvec != old_nsvec) ||
      (QOP_dw_nvec != old_nvec) ) {
    reset_temps(flw);
    old_style = QOP_dw_style;
    old_nsvec = QOP_dw_nsvec;
    old_nvec = QOP_dw_nvec;
  }

  if(dblstore_style(QOP_dw_style)) {
    wilson_dslash = wilson_dslash1;
  } else {
    wilson_dslash = wilson_dslash0;
  }

  ttv = (QDP_DiracFermion **) malloc(Ls*sizeof(QDP_DiracFermion *));
  for(i=0; i<Ls; i++) {
    ttv[i] = QDP_create_D();
    qdpin[i] = QDP_create_D();
    qdptmp[i] = in[i]->df;
    qdpout[i] = out[i]->df;
  }

  {
#ifdef LU
    QDP_D_eq_D(tt1, in->df, osubset);
    dslash_special_qdp(flw, ttt, tt1, 1, subset, 1);
    QDP_D_eq_r_times_D_plus_D(ttt, &kappa, ttt, in->df, subset);
    dw_mdslash1(flw, qdpin, ttt, -1, subset, osubset, mkappa);
#else
    dw_dslash1(flw, qdpin, qdptmp, -1, subset, osubset, m0, M, Ls);
#endif
  }

  printf0("begin cgv\n");
  dtime = -QOP_time();

  QOPPC(invert_cgv_D)(QOPPC(dw_dslash2), inv_arg, res_arg,
		      qdpout, qdpin, subset, Ls);

  dtime += QOP_time();
  printf0("end cgv\n");

#ifdef LU
  {
    QLA_Real kappa2 = 2.0*kappa;
    QDP_D_eq_D(ttt, out->df, subset);
    dslash_special_qdp(flw, tt2, ttt, 1, osubset, 2);
    QDP_D_eq_r_times_D_plus_D(ttt, &kappa, tt2, in->df, osubset);
    QDP_D_eq_r_times_D(out->df, &kappa2, ttt, QDP_all);
  }
#endif

  //res_arg->final_rsq = rsq;
  //res_arg->final_iter = iteration;
  //inv_arg->final_iter = iteration;
  inv_arg->final_sec = dtime;
  inv_arg->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;

  for(i=0; i<Ls; i++) {
    QDP_destroy_D(ttv[i]);
    QDP_destroy_D(qdpin[i]);
  }
  free(ttv);

  return QOP_SUCCESS;
}







/************ dslash *************/
static void set_out(void);

/* Special dslash for use by congrad.  Uses restart_gather() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather, otherwise use restart_gather.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
static void
wilson_dslash0(QOP_FermionLinksDW *flw,
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

  /* Take DW projection for src displaced in up direction, gather
     it to "our site" */

  printf0("dslash0\n");
  if(shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_D_veq_sD\n");
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, QDP_neighbor+mu, fwd+mu, subset,
		   QOP_dw_nsvec);
      printf0("end QDP_D_veq_sD\n");
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_H_veq_spproj_D\n");
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_dw_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, QDP_neighbor+mu, fwd+mu,
		   subset, QOP_dw_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Take DW projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_dw_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_dw_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_dw_nsvec);
#if 0
      QDP_H_veq_spproj_Ma_times_D(hf+mu, fwdlinks+mu, vsrc+mu,
                               dir+mu, msgn+mu, othersubset, QOP_dw_nsvec);
      QDP_D_veq_sprecon_H(dtemp[ntmp]+8+mu, hf+mu,
                               dir+mu, msgn+mu, othersubset, QOP_dw_nsvec);
#endif
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_dw_nsvec);
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_H_veq_spproj_Ma_times_D\n");
      QDP_H_veq_spproj_Ma_times_D(htemp[ntmp]+12+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_dw_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+12+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_dw_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Set dest to zero */
  /* Take DW projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  printf0("dslash0 - fwd\n");
  //QDP_D_eq_zero(dest, subset);
  set_out();

  if(shiftd_style(QOP_dw_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->links+mu, dtemp[ntmp]+mu,
				  dir+mu, sgn+mu, subset, QOP_dw_nvec);
#if 0
      QDP_H_veq_spproj_M_times_D(hf+mu, flw->links+mu, dtemp[ntmp]+mu,
				 dir+mu, sgn+mu, subset, QOP_dw_nvec);
      QDP_D_vpeq_sprecon_H(vdest+mu, hf+mu,
			   dir+mu, sgn+mu, subset, QOP_dw_nvec);
#endif
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->links+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_dw_nvec);
    }
  }

  /* Take DW projection for src displaced in down direction,
     expand it, and add to dest */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_D(vdest+mu, dtemp[ntmp]+4+mu, subset, QOP_dw_nvec);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_H(vdest+mu, htemp[ntmp]+4+mu, dir+mu, msgn+mu, subset,
			   QOP_dw_nvec);
    }
  }

  if(shiftd_style(QOP_dw_style)) {
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
wilson_dslash1(QOP_FermionLinksDW *flw,
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

  if(!shiftd_style(QOP_dw_style)) {
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

  /* Take DW projection for src displaced in up direction, gather
     it to "our site" */

  //printf0("ds1 1\n");
  if(shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<8; mu+=QOP_dw_nsvec) {
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, sh+mu, sd+mu, subset,
		   QOP_dw_nsvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_dw_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_dw_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, sh+mu, sd+mu, subset,
		   QOP_dw_nsvec);
    }
  }
  //printf0("ds1 2\n");

  /* Set dest to zero */
  /* Take DW projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add to dest */

  //QDP_D_eq_zero(dest, subset);
  set_out();

  //printf0("ds1 3\n");
  if(shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<8; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->dbllinks+mu, dtemp[ntmp]+mu,
			          dir+mu, sgn+mu, subset, QOP_dw_nvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->dbllinks+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_dw_nvec);
    }
  }
  //printf0("ds1 4\n");

  if(shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
} /* end of dslash_special_qdp() */

static QLA_Real mw;
static int dirs[2], signs[2], lupass=0;
static QDP_DiracFermion *vout[2], *vin[2], *out0, *in0;
static QDP_Subset lusubset;

static void
set_out(void)
{
#ifdef LU
  if(lupass==0) {
    QDP_D_eq_zero(out0, lusubset);
  } else {
    QDP_D_eq_r_times_D(out0, &mw, in0, lusubset);
    QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, lusubset, 2);
  }
#else
  //QDP_D_eq_r_times_D(out0, &mw, in0, QDP_all);
  //QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_all, 2);
  QDP_D_eq_r_times_D(out0, &mw, in0, QDP_even);
  QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_even, 2);
  QDP_D_eq_r_times_D(out0, &mw, in0, QDP_odd);
  QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_odd, 2);
#endif
}

static void
dw_dslash1(QOP_FermionLinksDW *flw,
	   QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	   int sign, QDP_Subset subset, QDP_Subset othersubset,
	   QLA_Real m0, QLA_Real M, int Ls)
{
#ifdef LU

  QLA_Real half = -0.5;
  QLA_Real mf = -M;
  int i;

  lupass = 0;
  lusubset = osubset;
  for(i=0; i<Ls; i++) {
    out0 = out[i];
    QDP_D_eq_r_times_D(tt1, &half, in[i], subset);
    wilson_dslash(flw, out[i], tt1, sign, osubset, 0);
  }

  //Qoo-1(out);

  lupass = 1;
  lusubset = subset;
  mw = m0;
  dirs[0] = dirs[1] = 4;
  signs[0] = -sign;
  signs[1] = sign;
  for(i=0; i<Ls; i++) {
    vout[0] = vout[1] = out[i];
    vin[0] = in[(i+1)%Ls];
    vin[1] = in[(i-1+Ls)%Ls];
    out0 = out[i];
    in0 = in[i];
    if(i==0) {
      QDP_D_eq_r_times_D(tt2, &mf, in[Ls-1], subset);
      vin[1] = tt2;
    } else if(i==Ls-1) {
      QDP_D_eq_r_times_D(tt2, &mf, in[0], subset);
      vin[0] = tt2;
    }

    QDP_D_eq_r_times_D(tt1, &half, out[i], osubset);
    wilson_dslash(flw, out[i], tt1, sign, subset, 1);
  }

#else

#if 0
  QLA_Real half = -0.5;
  QLA_Real mw[Ls];
  QLA_Real mf = -M;
  int dirs[2*Ls], signs[2*Ls];
  QDP_DiracFermion *vout[2*Ls], *vin[2*Ls];
  int i;
  for(i=0; i<Ls; i++) {
    QDP_D_eq_r_times_D(tt1, &half, in[i], QDP_all);
    wilson_dslash(flw, out[i], tt1, sign, QDP_all, 0);
    //QDP_D_eqm_D(out[i], out[i], QDP_all);
    mw[i] = m0;
    vout[2*i] = out[i];
    vin[2*i] = in[(i+1)%Ls];
    dirs[2*i] = 4;
    signs[2*i] = -sign;
    vout[2*i+1] = out[i];
    vin[2*i+1] = in[(i-1+Ls)%Ls];
    dirs[2*i+1] = 4;
    signs[2*i+1] = sign;
  }
  QDP_D_eq_r_times_D(tt1, &mf, in[0], QDP_all);
  QDP_D_eq_r_times_D(tt2, &mf, in[Ls-1], QDP_all);
  vin[2*(Ls-1)] = tt1;
  vin[1] = tt2;
  QDP_D_vpeq_r_times_D(out, mw, in, QDP_all, Ls);
  QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_all, 2*Ls);

#else

  QLA_Real half = -0.5;
  QLA_Real mf = -M;
  int i;

  mw = m0;
  dirs[0] = dirs[1] = 4;
  signs[0] = -sign;
  signs[1] = sign;
  for(i=0; i<Ls; i++) {
    vout[0] = vout[1] = out[i];
    vin[0] = in[(i+1)%Ls];
    vin[1] = in[(i-1+Ls)%Ls];
    out0 = out[i];
    in0 = in[i];
    if(i==0) {
      QDP_D_eq_r_times_D(tt2, &mf, in[Ls-1], QDP_all);
      vin[1] = tt2;
    } else if(i==Ls-1) {
      QDP_D_eq_r_times_D(tt2, &mf, in[0], QDP_all);
      vin[0] = tt2;
    }

    QDP_D_eq_r_times_D(tt1, &half, in[i], QDP_all);
    //QDP_D_eq_r_times_D(out[i], &mw, in[i], QDP_all);
    wilson_dslash(flw, out[i], tt1, sign, QDP_all, 0);
    //QDP_D_eqm_D(out[i], out[i], QDP_all);

    //QDP_D_peq_r_times_D(out[i], &mw, in[i], QDP_all);
    //QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_all, 2);
  }

#endif

#endif
}

static void
dw_dslash2(QOP_FermionLinksDW *flw,
	   QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	   QDP_Subset subset, QDP_Subset othersubset,
	   QLA_Real m0, QLA_Real M, int Ls)
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

  dw_dslash1(flw, ttv, in, 1, subset, othersubset, m0, M, Ls);
  dw_dslash1(flw, out, ttv, -1, subset, othersubset, m0, M, Ls);

#endif
}
