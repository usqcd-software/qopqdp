#include <qop_internal.h>

static int old_style=0;
static int old_nsvec=0;
static int old_nvec=0;

static int congrad_setup = 0;
static QDP_ColorVector *tttt;
static QDP_ColorVector *temp1[24], *temp2[24];

static void QOPPC(asqtad_mdslash2)(QOPPC(FermionLinksAsqtad) *,
				   QDP_ColorVector *out, QDP_ColorVector *in,
				   QDP_Subset subset, QDP_Subset othersubset,
				   QLA_Real m2x4);

static void
free_temps(void)
{
  if(congrad_setup) {
    int i, n;

    QDP_destroy_V(tttt);

    if(old_style==0) n = 24;
    else n = 16;
    for(i=0; i<n; i++) {
      QDP_destroy_V(temp1[i]);
      QDP_destroy_V(temp2[i]);
    }
  }
  congrad_setup = 0;
}

static void
double_store(QOPPC(FermionLinksAsqtad) *fla)
{
  if( (QOP_asqtad.style==1) && (!fla->dblstored) ) {
    int i;
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<8; i++) {
      QDP_M_eq_sM(m, fla->fwdlinks[i], QOP_asqtad.shifts[i],
		  QDP_backward, QDP_all);
      QDP_M_eqm_Ma(fla->bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    fla->dblstored = 1;
  }
}

static void
reset_temps(QOPPC(FermionLinksAsqtad) *fla)
{
  int i, n;

  if(QOP_asqtad.style!=old_style) {
    if(QOP_asqtad.style==0) {
      if(congrad_setup) {
	for(i=0; i<8; i++) {
	  QDP_destroy_M(fla->bcklinks[i]);
	}
      }
    } else {
      for(i=0; i<8; i++) {
	fla->bcklinks[i] = QDP_create_M();
      }
      for(i=0; i<4; i++) {
	fla->dbllinks[4*i] = fla->fwdlinks[2*i];
	fla->dbllinks[4*i+1] = fla->fwdlinks[2*i+1];
	fla->dbllinks[4*i+2] = fla->bcklinks[2*i];
	fla->dbllinks[4*i+3] = fla->bcklinks[2*i+1];
      }
    }
    fla->dblstored = 0;
  }
  double_store(fla);

  free_temps();

  tttt = QDP_create_V();

  if(QOP_asqtad.style==0) n = 24;
  else n = 16;
  for(i=0; i<n; i++) {
    temp1[i] = QDP_create_V();
    temp2[i] = QDP_create_V();
  }
  congrad_setup = 1;
}


/********************/
/*  link functions  */
/********************/

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_convert_L_from_qdp)(QDP_ColorMatrix *fatlinks[],
				 QDP_ColorMatrix *longlinks[])
{
  int i;
  QOPPC(FermionLinksAsqtad) *fla;
  if(!QOP_asqtad.inited) QOP_asqtad_invert_init();
  QOP_malloc(fla, QOPPC(FermionLinksAsqtad), 1);
  QOP_malloc(fla->fatlinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(fla->longlinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(fla->fwdlinks, QDPPC(ColorMatrix) *, 8);
  QOP_malloc(fla->bcklinks, QDPPC(ColorMatrix) *, 8);
  QOP_malloc(fla->dbllinks, QDPPC(ColorMatrix) *, 16);
  for(i=0; i<4; i++) {
    fla->fatlinks[i] = fatlinks[i];
    fla->longlinks[i] = longlinks[i];
  }

  fla->dblstored = 0;
  for(i=0; i<4; i++) {
    fla->fwdlinks[2*i] = fatlinks[i];
    fla->fwdlinks[2*i+1] = longlinks[i];
  }
  congrad_setup = 0;
  old_style = 0;

  if( (!congrad_setup) ||
      (QOP_asqtad.style != old_style) ||
      (QOP_asqtad.nsvec != old_nsvec) ||
      (QOP_asqtad.nvec != old_nvec) ) {
    reset_temps(fla);
    old_style = QOP_asqtad.style;
    old_nsvec = QOP_asqtad.nsvec;
    old_nvec = QOP_asqtad.nvec;
  }

  double_store(fla);

  return fla;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_qdp)(QDP_ColorMatrix *fatlinks[],
				QDP_ColorMatrix *longlinks[])
{
  QOPPC(FermionLinksAsqtad) *fla;
  QDP_ColorMatrix *fl[4], *ll[4];
  int i;

  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    QDP_M_eq_M(fl[i], fatlinks[i], QDP_all);
    ll[i] = QDP_create_M();
    QDP_M_eq_M(ll[i], longlinks[i], QDP_all);
  }

  fla = QOPPC(asqtad_convert_L_from_qdp)(fl, ll);

  return fla;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_raw)(REAL *fatlinks[], REAL *longlinks[],
				QOP_evenodd_t evenodd)
{
  QDP_ColorMatrix *fl[4], *ll[4];
  int i;
  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    QOP_qdp_eq_raw(M, fl[i], fatlinks[i], evenodd);
    ll[i] = QDP_create_M();
    QOP_qdp_eq_raw(M, ll[i], longlinks[i], evenodd);
  }
  return QOPPC(asqtad_create_L_from_qdp)(fl, ll);
}

/* Computes the staple :
                 mu
              +-------+
        nu    |       |
              |       |
              X       X
  Where the mu link can be any su3_matrix. The result is saved in staple.
  if staple==NULL then the result is not saved.
  It also adds the computed staple to the fatlink[mu] with weight coef.
*/
static void
compute_gen_staple(QDP_ColorMatrix *staple, int mu, int nu,
		   QDP_ColorMatrix *link, double dcoef,
		   QDP_ColorMatrix **gauge, QDP_ColorMatrix **fl)
{
  QLA_Real coef = dcoef;
  QDP_ColorMatrix *ts0, *ts1;
  QDP_ColorMatrix *tmat1, *tmat2;
  QDP_ColorMatrix *tempmat;

  ts0 = QDP_create_M();
  ts1 = QDP_create_M();
  tmat1 = QDP_create_M();
  tmat2 = QDP_create_M();
  tempmat = QDP_create_M();

  /* Upper staple */
  QDP_M_eq_sM(ts0, link, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);

  if(staple!=NULL) {  /* Save the staple */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(staple, gauge[nu], tmat1, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(tmat2, gauge[nu], tmat1, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, tmat2, QDP_all);
  }

  /* lower staple */
  QDP_M_eq_sM(ts0, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(tmat1, gauge[nu], link, QDP_all);
  QDP_M_eq_M_times_M(tempmat, tmat1, ts0, QDP_all);
  QDP_M_eq_sM(ts0, tempmat, QDP_neighbor[nu], QDP_backward, QDP_all);

  if(staple!=NULL) { /* Save the staple */
    QDP_M_peq_M(staple, ts0, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, staple, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_peq_r_times_M(fl[mu], &coef, ts0, QDP_all);
  }

  QDP_destroy_M(ts0);
  QDP_destroy_M(ts1);
  QDP_destroy_M(tmat1);
  QDP_destroy_M(tmat2);
  QDP_destroy_M(tempmat);
} /* compute_gen_staple */

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_G)(QOP_asqtad_coeffs_t *coeffs,
			      QOPPC(GaugeField) *gauge)
{
  QOPPC(FermionLinksAsqtad) *fla;
  QDP_ColorMatrix *fl[4], *ll[4], **gf;
  int i, dir;
  QDP_ColorMatrix *staple, *tempmat1;
  int  nu,rho,sig ;
  QLA_Real one_link;
  int nflop = 61632;
  double dtime;

  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    ll[i] = QDP_create_M();
  }
  staple = QDP_create_M();
  tempmat1 = QDP_create_M();
  gf = gauge->links;

  dtime = -QOP_time();

  /* to fix up the Lepage term, included by a trick below */
  one_link = coeffs->one_link - 6.0*coeffs->lepage;

  for(dir=0; dir<4; dir++) {
    QDP_M_eq_r_times_M(fl[dir], &one_link, gf[dir], QDP_all);
    for(nu=0; nu<4; nu++) if(nu!=dir) {
      compute_gen_staple(staple, dir, nu, gf[dir],
			 coeffs->three_staple, gf, fl);
      compute_gen_staple(NULL, dir, nu, staple, coeffs->lepage, gf, fl);
      for(rho=0; rho<4; rho++) if((rho!=dir)&&(rho!=nu)) {
	compute_gen_staple(tempmat1, dir, rho, staple,
			   coeffs->five_staple, gf, fl);
	for(sig=0; sig<4; sig++) {
	  if((sig!=dir)&&(sig!=nu)&&(sig!=rho)) {
	    compute_gen_staple(NULL, dir, sig, tempmat1,
			       coeffs->seven_staple, gf, fl);
	  }
	} /* sig */
      } /* rho */
    } /* nu */
  } /* dir */

  /* long links */
  for(dir=0; dir<4; dir++) {
    QLA_Real naik = coeffs->naik;
    QDP_M_eq_sM(staple, gf[dir], QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(tempmat1, gf[dir], staple, QDP_all);
    QDP_discard_M(staple);
    QDP_M_eq_sM(staple, tempmat1, QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(ll[dir], gf[dir], staple, QDP_all);
    QDP_M_eq_r_times_M(ll[dir], &naik, ll[dir], QDP_all);
  }

  QDP_destroy_M(staple);
  QDP_destroy_M(tempmat1);

  dtime += QOP_time();
  //node0_printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,
  //       (Real)nflop*volume/(1e6*dtime*numnodes()) );

  fla = QOPPC(asqtad_convert_L_from_qdp)(fl, ll);

  return fla;
}

void
QOPPC(asqtad_extract_L_to_raw)(REAL *fatlinks[], REAL *longlinks[],
			       QOPPC(FermionLinksAsqtad) *src,
			       QOP_evenodd_t evenodd)
{
  int i;
  for(i=0; i<4; i++) {
    QOP_raw_eq_qdp(M, fatlinks[i], src->fatlinks[i], evenodd);
    QOP_raw_eq_qdp(M, longlinks[i], src->longlinks[i], evenodd);
  }
}

void
QOPPC(asqtad_destroy_L)(QOPPC(FermionLinksAsqtad) *field)
{
  int i;
  for(i=0; i<4; i++) {
    QDP_destroy_M(field->fatlinks[i]);
    QDP_destroy_M(field->longlinks[i]);
  }
}


/* inverter */

QOPPC(linop_t_V) QOP_asqtad_dslash2;
static QOPPC(FermionLinksAsqtad) *gl_fla;
static QDP_Subset gl_osubset;
static QLA_Real gl_m2x4;

void
QOPPC(asqtad_dslash2)(QDP_ColorVector *out, QDP_ColorVector *in,
		      QDP_Subset subset)
{
  QOPPC(asqtad_mdslash2)(gl_fla, out, in, subset, gl_osubset, gl_m2x4);
}

QOP_status_t
QOPPC(asqtad_invert)(QOPPC(FermionLinksAsqtad) *fla,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     REAL mass,
		     QOPPC(ColorVector) *out,
		     QOPPC(ColorVector) *in)
{
  /* cg has 5 * 12 = 60 flops/site/it */
  /* MdagM -> 2*(66+72*15)+12 = 2304 flops/site */
  double dtimec;
  double nflop = 0.5 * 2364;
  QDP_Subset subset, othersubset;

  if(inv_arg->evenodd==QOP_EVEN) {
    subset = QDP_even;
    othersubset = QDP_odd;
  } else if(inv_arg->evenodd==QOP_ODD) {
    subset = QDP_odd;
    othersubset = QDP_even;
  } else {
    subset = QDP_all;
    othersubset = QDP_all;
    nflop *= 2;
  }
  gl_osubset = othersubset;

  gl_m2x4 = 4*mass*mass;
  gl_fla = fla;

  if( (QOP_asqtad.style != old_style) ||
      (QOP_asqtad.nsvec != old_nsvec) ||
      (QOP_asqtad.nvec != old_nvec) ) {
    reset_temps(fla);
    old_style = QOP_asqtad.style;
    old_nsvec = QOP_asqtad.nsvec;
    old_nvec = QOP_asqtad.nvec;
  }
  //QOPPC(invert_cg_init_V)();

  dtimec = -QDP_time();

  QOPPC(invert_cg_V)(QOPPC(asqtad_dslash2), inv_arg, res_arg,
		     out->cv, in->cv, subset);

  dtimec += QDP_time();

  //res_arg->final_rsq = rsq;
  //res_arg->final_iter = iteration;
  //inv_arg->final_iter = iteration;
  inv_arg->final_sec = dtimec;
  inv_arg->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;

  return QOP_SUCCESS;
}

QOP_status_t
QOPPC(asqtad_invert_multi)(QOP_FermionLinksAsqtad *asqtad,
			   QOP_invert_arg_t *inv_arg,
			   QOP_resid_arg_t **res_arg[],
			   REAL *masses[], int nmass[],
			   QOP_ColorVector **out_pt[],
			   QOP_ColorVector *in_pt[],
			   int nsrc)
{
  QOP_invert_arg_t tinvarg;
  int i, j;

  //inv_arg->final_iter = 0;
  inv_arg->final_sec = 0;
  inv_arg->final_flop = 0;
  for(i=0; i<nsrc; i++) {
    for(j=0; j<nmass[i]; j++) {
      QOPPC(asqtad_invert)(asqtad, &tinvarg, res_arg[i][j], masses[i][j],
			   out_pt[i][j], in_pt[i]);
      //inv_arg->final_iter += tinvarg.final_iter;
      inv_arg->final_sec += tinvarg.final_sec;
      inv_arg->final_flop += tinvarg.final_flop;
    }
  }

  return QOP_SUCCESS;
}




/* internal functions */

static void
QOPPC(asqtad_dslash0)(QOPPC(FermionLinksAsqtad) *fla,
		      QDP_ColorVector *dest, QDP_Subset dsubset,
		      QDP_ColorVector *src, QDP_Subset ssubset,
		      QDP_ColorVector *temp[])
{
  QDP_ColorVector *vsrc[8], *vdest[8];
  int i;

  for(i=0; i<8; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  /* Start gathers from positive directions */
  for(i=0; i<8; i+=QOP_asqtad.nsvec) {
    QDP_V_veq_sV(temp+i, vsrc+i, QOP_asqtad.shifts+i, QOP_common.shiftfwd+i,
		 dsubset, QOP_asqtad.nsvec);
  }

  /* Multiply by adjoint matrix at other sites */
  /* Start gathers from negative directions */
  for(i=0; i<8; i+=QOP_asqtad.nsvec) {
    QDP_V_veq_Ma_times_V(temp+16+i, fla->fwdlinks+i, vsrc+i, ssubset,
			 QOP_asqtad.nsvec);
    QDP_V_veq_sV(temp+8+i, temp+16+i, QOP_asqtad.shifts+i,
		 QOP_common.shiftbck+i, dsubset, QOP_asqtad.nsvec);
  }

  /* Multiply by matrix for positive directions and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<8; i+=QOP_asqtad.nvec) {
    QDP_V_vpeq_M_times_V(vdest+i, fla->fwdlinks+i, temp+i, dsubset,
			 QOP_asqtad.nvec);
  }
  for(i=0; i<8; i++) QDP_discard_V(temp[i]);

  /* Wait gathers from negative directions, accumulate (negative) */
  for(i=0; i<8; i+=QOP_asqtad.nvec) {
    QDP_V_vmeq_V(vdest+i, temp+8+i, dsubset, QOP_asqtad.nvec);
  }
  for(i=0; i<8; ++i) QDP_discard_V(temp[i+8]);
}

static void
QOPPC(asqtad_dslash1)(QOPPC(FermionLinksAsqtad) *fla,
		      QDP_ColorVector *dest, QDP_Subset dsubset,
		      QDP_ColorVector *src, QDP_Subset ssubset,
		      QDP_ColorVector *temp[])
{
  QDP_ColorVector *vsrc[16], *vdest[16];
  int i;

  for(i=0; i<16; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  /* Start gathers from all directions */
  for(i=0; i<16; i+=QOP_asqtad.nsvec) {
    QDP_V_veq_sV(temp+i, vsrc+i, QOP_asqtad.shifts_dbl+i,
		 QOP_asqtad.shiftdirs_dbl+i, dsubset, QOP_asqtad.nsvec);
  }

  /* Wait gathers from all directions, multiply by matrix and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<16; i+=QOP_asqtad.nvec) {
    QDP_V_vpeq_M_times_V(vdest+i, fla->dbllinks+i, temp+i, dsubset,
			 QOP_asqtad.nvec);
  }
  for(i=0; i<16; ++i) QDP_discard_V(temp[i]);
}

static void
QOPPC(asqtad_mdslash2)(QOPPC(FermionLinksAsqtad) *fla,
		       QDP_ColorVector *out, QDP_ColorVector *in,
		       QDP_Subset subset, QDP_Subset othersubset,
		       QLA_Real m2x4)
{
  if(QOP_asqtad.style==0) {
    QOPPC(asqtad_dslash0)(fla, tttt, othersubset, in, subset, temp1);
    QOPPC(asqtad_dslash0)(fla, out, subset, tttt, othersubset, temp2);
  } else {
    QOPPC(asqtad_dslash1)(fla, tttt, othersubset, in, subset, temp1);
    QOPPC(asqtad_dslash1)(fla, out, subset, tttt, othersubset, temp2);
  }
  //QDP_V_meq_r_times_V(out, &m2x4, in, subset);
  QDP_V_eq_r_times_V_minus_V(out, &m2x4, in, out, subset);
}
