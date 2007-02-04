#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU

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
static QDP_HalfFermion *htemp[2][16];
static QDP_DiracFermion *dtemp[2][12];
static QDP_HalfFermion **hf0, **hf1, *tth;
static QDP_DiracFermion *tt1, *tt2, **ttv;

static void
free_temps(QOP_FermionLinksDW *flw)
{
  if(congrad_setup) {
    int i, j;

    QDP_destroy_H(tth);
    QDP_destroy_D(tt1);
    QDP_destroy_D(tt2);

    if(shiftd_style(old_style)) {
      for(i=0; i<2; i++) {
	for(j=0; j<12; j++) {
	  QDP_destroy_D(dtemp[i][j]);
	}
      }
    } else {
      for(i=0; i<2; i++) {
	for(j=0; j<16; j++) {
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

  tth = QDP_create_H();
  tt1 = QDP_create_D();
  tt2 = QDP_create_D();

  if(shiftd_style(QOP_dw_style)) {
    for(i=0; i<2; i++) {
      for(j=0; j<12; j++) {
	dtemp[i][j] = QDP_create_D();
      }
    }
  } else {
    for(i=0; i<2; i++) {
      for(j=0; j<16; j++) {
	htemp[i][j] = QDP_create_H();
      }
    }
  }
  congrad_setup = 1;
}


/* link routines */

QOP_FermionLinksDW *
QOP_dw_create_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
{
  QOP_FermionLinksDW *flw;
  QOP_GaugeField *gf;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_dw_convert_L_from_qdp(gf->links, NULL);

  flw->raw = NULL;
  flw->qopgf = gf;

  return flw;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_G(QOP_info_t *info, QOP_dw_coeffs_t *coeffs,
		       QOP_GaugeField *gauge)
{
  QOP_error("QOP_dw_create_L_from_G unimplemented");
  return NULL;
}

void
QOP_dw_extract_L_to_raw(REAL *links[], REAL *clov, QOP_FermionLinksDW *src,
			QOP_evenodd_t evenodd)
{
  QOP_error("QOP_dw_extract_L_to_raw unimplemented");
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
QOP_dw_convert_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
{
  QOP_error("QOP_dw_convert_L_from_raw unimplemented");
  return NULL;
}

void
QOP_dw_convert_L_to_raw(REAL ***links, REAL **clov, QOP_FermionLinksDW *src,
			QOP_evenodd_t evenodd)
{
  QOP_error("QOP_dw_convert_L_to_raw unimplemented");
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_G(QOP_info_t *info, QOP_dw_coeffs_t *coeffs,
			QOP_GaugeField *gauge)
{
  QOP_error("QOP_dw_convert_L_from_G unimplemented");
  return NULL;
}

QOP_GaugeField *
QOP_dw_convert_L_to_G(QOP_FermionLinksDW *links)
{
  QOP_error("QOP_dw_convert_L_to_G unimplemented");
  return NULL;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_qdp(QDP_ColorMatrix *links[], QDP_DiracPropagator *clov)
{
  QOP_FermionLinksDW *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_dw_convert_L_from_qdp(newlinks, clov);

  return flw;
}

void
QOP_dw_extract_L_to_qdp(QDP_ColorMatrix *links[],
			QDP_DiracPropagator *clov,
			QOP_FermionLinksDW *src)
{
  QOP_error("QOP_dw_extract_L_to_qdp unimpleented");
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_qdp(QDP_ColorMatrix *links[],
			  QDP_DiracPropagator *clov)
{
  QOP_FermionLinksDW *flw;
  int i;

  if(clov!=NULL) {
    QOP_error("clover term is not implemented yet.");
  }

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

void
QOP_dw_convert_L_to_qdp(QDP_ColorMatrix ***link, QDP_DiracPropagator **clov,
			QOP_FermionLinksDW *src)
{
  QOP_error("QOP_dw_convert_L_to_qdp unimpleented");
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
		    int Ls,
		    int nsrc)
{
  QOP_error("QOP_dw_invert_multi unimplemented");
}


static void Qeo(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
		QDP_DiracFermion *in[], REAL m0, REAL M, int Ls);
static void Qooinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		   REAL m0, REAL M, int Ls);
static void Qeeinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
		   REAL m0, REAL M, int Ls);
static void QoeQeeinv(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
		      QDP_DiracFermion *in[], REAL m0, REAL M, int Ls);

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

void
QOPPC(dw_invert)(QOP_info_t *info,
		 QOP_FermionLinksDW *flw,
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
  QDP_DiracFermion *qdptmp[Ls], *qdpin[Ls], *qdpout[Ls], *ttv2[Ls];
  QDP_Subset subset, osubset;
  int i;

  /* cg has 5 * 48 *Ls = 240*Ls flops/site/it */
#ifdef LU
  /*              mf     D_w   1st Qxx1  1     */
  /* MdagM -> 2*(-48 + ((168*8-24)+4*24+12)*Ls) = -96 + 2856*Ls flops/site */
  nflop = -96 + 2976*Ls; /* cg only on 1/2 the sites */
  subset = QDP_odd;
  osubset = QDP_even;
#else
  /*              mf       D_w   m0 .5 Ds     */
  /* MdagM -> 2*(2*24 + ((168*8)+24+24+24)*Ls) = 96 + 2832*Ls flops/site */
  nflop = 96 + 3072*Ls;
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

  hf0 = (QDP_HalfFermion **) malloc(Ls*sizeof(QDP_HalfFermion *));
  hf1 = (QDP_HalfFermion **) malloc(Ls*sizeof(QDP_HalfFermion *));
  ttv = (QDP_DiracFermion **) malloc(Ls*sizeof(QDP_DiracFermion *));
  for(i=0; i<Ls; i++) {
    hf0[i] = QDP_create_H();
    hf1[i] = QDP_create_H();
    ttv[i] = QDP_create_D();
    ttv2[i] = QDP_create_D();
    qdpin[i] = QDP_create_D();
    qdptmp[i] = in[i]->df;
    qdpout[i] = out[i]->df;
  }

  {
#ifdef LU
    QoeQeeinv(flw, qdpin, qdptmp, m0, M, Ls);
    for(i=0; i<Ls; i++) QDP_D_eq_D_minus_D(qdpin[i], qdptmp[i], qdpin[i], QDP_odd);
    Qooinv(ttv, qdpin, m0, M, Ls);
    dw_dslash1(flw, qdpin, ttv, -1, subset, osubset, m0, M, Ls);
#else
    dw_dslash1(flw, qdpin, qdptmp, -1, subset, osubset, m0, M, Ls);
#endif
  }

#ifndef LU
  res_arg->rsqmin /= m0*m0*m0*m0;
#endif
  printf0("begin cgv\n");
  dtime = -QOP_time();

  QOPPC(invert_cg_vD)(QOPPC(dw_dslash2), inv_arg, res_arg,
		      qdpout, qdpin, ttv2, subset, Ls);

  dtime += QOP_time();
  printf0("end cgv\n");
#ifndef LU
  res_arg->rsqmin *= m0*m0*m0*m0;
#endif

#ifdef LU
  {
    Qeo(flw, ttv, qdpout, m0, M, Ls);
    for(i=0; i<Ls; i++) QDP_D_eq_D_minus_D(ttv[i], qdptmp[i], ttv[i], QDP_even);
    Qeeinv(qdpout, ttv, m0, M, Ls);
  }
#endif

  for(i=0; i<Ls; i++) {
    QDP_destroy_H(hf0[i]);
    QDP_destroy_H(hf1[i]);
    QDP_destroy_D(ttv[i]);
    QDP_destroy_D(ttv2[i]);
    QDP_destroy_D(qdpin[i]);
  }
  free(hf0);
  free(hf1);
  free(ttv);

  //res_arg->final_rsq = rsq;
  //res_arg->final_iter = iteration;
  //inv_arg->final_iter = iteration;
  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}




/************ dslashes *************/

static void set_out(void);

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
  //sgn[1] = -sign;
  //msgn[1] = sign;

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
}

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
  //sgn[2] = -sign;
  //sgn[3] = sign;

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
}


/* domain wall operators */

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
  QDP_D_eq_r_times_D(out0, &mw, in0, QDP_all);
  QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_all, 2);
  //QDP_D_eq_r_times_D(out0, &mw, in0, QDP_even);
  //QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_even, 2);
  //QDP_D_eq_r_times_D(out0, &mw, in0, QDP_odd);
  //QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_odd, 2);
#endif
}

static void
Qxy(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[], QDP_DiracFermion *in[],
    REAL m0, REAL M, int Ls, int sign, QDP_Subset subset, QDP_Subset osubset)
{
  QLA_Real half = -0.5;
  int i, ntmp;
  if(subset==QDP_odd) ntmp = 0; else ntmp = 1;
  lupass = 0;
  lusubset = subset;
  for(i=0; i<Ls; i++) {
    out0 = out[i];
    QDP_D_eq_r_times_D(tt1, &half, in[i], osubset);
    wilson_dslash(flw, out[i], tt1, sign, subset, ntmp);
  }
}

static void
Qeo(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
    QDP_DiracFermion *in[], REAL m0, REAL M, int Ls)
{
  Qxy(flw, out, in, m0, M, Ls, 1, QDP_even, QDP_odd);
}

static void
Qoe(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
    QDP_DiracFermion *in[], REAL m0, REAL M, int Ls)
{
  Qxy(flw, out, in, m0, M, Ls, 1, QDP_odd, QDP_even);
}

static void
Seo(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
    QDP_DiracFermion *in[], REAL m0, REAL M, int Ls)
{
  Qxy(flw, out, in, m0, M, Ls, -1, QDP_even, QDP_odd);
}

static void
Soe(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
    QDP_DiracFermion *in[], REAL m0, REAL M, int Ls)
{
  Qxy(flw, out, in, m0, M, Ls, -1, QDP_odd, QDP_even);
}

static void
Qxxinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int Ls, int sign, QDP_Subset subset)
{
  QLA_Real fac, den, mba;
  QDP_DiracFermion *tin[2];
  QDP_HalfFermion *thf[2];
  int dirs[2]={4,4}, sgns[2]={1,-1};
  int i;

  sign = -sign;
  sgns[0] *= sign;
  sgns[1] *= sign;

  mba = 1.0/m0;
  den = 1.0/(1.0 + (M/m0)*pow(mba, Ls-1));
  fac = -(M/m0)*den;
  QDP_H_eq_zero(hf0[Ls-1], subset);
  for(i=0; i<Ls; i++) {
    thf[0] = hf0[i];
    if(i==Ls-1) thf[0] = tth;
    thf[1] = hf1[i];
    tin[0] = in[i];
    tin[1] = in[i];
    QDP_H_veq_spproj_D(thf, tin, dirs, sgns, subset, 2);

    // L_A^-1
    if(i<Ls-1) {
      QDP_H_peq_r_times_H(hf0[Ls-1], &fac, hf0[i], subset);
      fac *= mba;
    } else {
      QDP_H_peq_r_times_H(hf0[Ls-1], &den, tth, subset);
    }

    // L_B^-1
    if(i==0) {
      QDP_H_eq_r_times_H(hf1[i], &mba, hf1[i], subset);
    } else {
      QDP_H_peq_H(hf1[i], hf1[i-1], subset);
      QDP_H_eq_r_times_H(hf1[i], &mba, hf1[i], subset);
    }
  }

  fac = -(M/m0)*pow(mba, Ls-2);
  for(i=Ls-1; i>=0; i--) {
    // R_A^-1
    if(i==Ls-1) {
      QDP_H_eq_r_times_H(hf0[i], &mba, hf0[i], subset);
    } else {
      QDP_H_peq_H(hf0[i], hf0[i+1], subset);
      QDP_H_eq_r_times_H(hf0[i], &mba, hf0[i], subset);
    }

    // R_B^-1
    if(i==Ls-1) {
      QDP_H_eq_r_times_H(hf1[i], &den, hf1[i], subset);
    } else {
      QDP_H_peq_r_times_H(hf1[i], &fac, hf1[Ls-1], subset);
      fac *= m0;
    }

    QDP_D_eq_sprecon_H(out[i], hf0[i], dirs[0], sgns[0], subset);
    QDP_D_peq_sprecon_H(out[i], hf1[i], dirs[1], sgns[1], subset);
  }
}

static void
Qooinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int Ls)
{
  Qxxinv(out, in, m0, M, Ls, 1, QDP_odd);
}

static void
Qeeinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int Ls)
{
  Qxxinv(out, in, m0, M, Ls, 1, QDP_even);
}

static void
Sooinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int Ls)
{
  Qxxinv(out, in, m0, M, Ls, -1, QDP_odd);
}

static void
Seeinv(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int Ls)
{
  Qxxinv(out, in, m0, M, Ls, -1, QDP_even);
}

static void
QoeQeeinv(QOP_FermionLinksDW *flw, QDP_DiracFermion *out[],
	  QDP_DiracFermion *in[], REAL m0, REAL M, int Ls)
{
  Qeeinv(out, in, m0, M, Ls);
  Qoe(flw, out, out, m0, M, Ls);
}

#if 0
static void
Qxx(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
    QLA_Real m0, QLA_Real M, int Ls, int sign, QDP_Subset subset)
{
  int i;
  sign = -sign;
  for(i=0; i<Ls; i++) {
    QDP_D_eq_r_times_D(out[i], &m0, in[i], subset);

    if(i<Ls-1) {
      QDP_D_meq_spproj_D(out[i], in[i+1], 4, sign, subset);
    } else {
      QDP_D_eq_r_times_D(tt1, &M, in[0], subset);
      QDP_D_peq_spproj_D(out[i], tt1, 4, sign, subset);
    }
    if(i==0) {
      QDP_D_eq_r_times_D(tt1, &M, in[Ls-1], subset);
      QDP_D_peq_spproj_D(out[i], tt1, 4, -sign, subset);
    } else {
      QDP_D_meq_spproj_D(out[i], in[i-1], 4, -sign, subset);
    }
  }
}

static void
Qoo(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
    QLA_Real m0, QLA_Real M, int Ls)
{
  Qxx(out, in, m0, M, Ls, 1, QDP_odd);
}
#endif

static void
dw_dslash1(QOP_FermionLinksDW *flw,
	   QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	   int sign, QDP_Subset subset, QDP_Subset othersubset,
	   QLA_Real m0, QLA_Real M, int Ls)
{
#ifdef LU
  int i;

#if 0
  Qoo(out, in, m0, M, Ls);
  Qooinv(out, out, m0, M, Ls);
  for(i=0; i<Ls; i++) {
    QLA_Real t;
    QDP_D_meq_D(out[i], in[i], QDP_odd);
    QDP_r_eq_norm2_D(&t, out[i], QDP_odd);
    printf("%g\n", t);
  }
#endif

  if(sign>0) {
    Qeo(flw, out, in, m0, M, Ls);
    Qeeinv(out, out, m0, M, Ls);
    Qoe(flw, out, out, m0, M, Ls);
    Qooinv(out, out, m0, M, Ls);
  } else {
    Sooinv(out, in, m0, M, Ls);
    Seo(flw, out, out, m0, M, Ls);
    Seeinv(out, out, m0, M, Ls);
    Soe(flw, out, out, m0, M, Ls);
  }
  for(i=0; i<Ls; i++) {
    QDP_D_eq_D_minus_D(out[i], in[i], out[i], subset);
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
    //wilson_dslash(flw, out[i], in[i], sign, QDP_all, 0);
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
  dw_dslash1(flw, ttv, in, 1, subset, othersubset, m0, M, Ls);
  dw_dslash1(flw, out, ttv, -1, subset, othersubset, m0, M, Ls);
}
