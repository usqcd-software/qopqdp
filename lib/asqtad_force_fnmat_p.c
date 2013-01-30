//#define DO_TRACE
#include <qop_internal.h>

static void
QOP_asqtad_force_multi_fnmat3_qdp(QOP_info_t *info, QOP_GaugeField *gauge,
				  QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
				  QDP_ColorMatrix *mid_fat[], QDP_ColorMatrix *mid_naik[]);

void 
QOP_asqtad_force_multi_fnmat_qdp(QOP_info_t *info, QOP_GaugeField *gauge,
				 QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
				 REAL eps[], QDP_ColorVector *x[], int nterms)
{
#define NC QDP_get_nc(x[0])
  ASQTAD_FORCE_BEGIN;
  if(!QOP_asqtad.inited) QOP_asqtad_invert_init();

  double dtime = QOP_time();
  double nflops = 0;
  QOP_info_t tinfo;

  QDP_ColorMatrix *mid[8];  // 4 fat and 4 naik
  for(int i=0; i<8; i++) mid[i] = QDP_create_M();

  QOP_get_mid(&tinfo, mid, QDP_neighbor, 4, eps, x, nterms);
  nflops += tinfo.final_flop;
  QOP_get_mid(&tinfo, mid+4, QOP_common.neighbor3, 4, eps, x, nterms);
  nflops += tinfo.final_flop;
  // compensate for -1 on odd sites here instead of at end
  for(int dir=0; dir<4; dir++) {
    QDP_M_eqm_M(mid[dir], mid[dir], QDP_odd);
    QDP_M_eqm_M(mid[4+dir], mid[4+dir], QDP_odd);
  }

  QOP_asqtad_force_multi_fnmat3_qdp(&tinfo, gauge, force, coef, mid, mid+4);
  nflops += tinfo.final_flop;

  for(int i=0; i<8; i++) QDP_destroy_M(mid[i]);

  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops;
  info->status = QOP_SUCCESS;

  ASQTAD_FORCE_END;
#undef NC
}

// gather vectors for mid link
void 
QOP_get_mid(QOP_info_t *info, QDP_ColorMatrix *mid[], QDP_Shift shifts[], int ns,
	    REAL eps[], QDP_ColorVector *x[], int nterms)
{
#define NC QDP_get_nc(x[0])
  double dtime = QOP_time();

#if 0
  {
    for(int i=0; i<nterms; i++) {
      QLA_Real nrm2;
      QDP_r_eq_norm2_V(&nrm2, x[i], QDP_all);
      QOP_printf0("norm2 x[%i] = %g\n", i, nrm2);
    }
    int i;
    QDP_loop_sites(i, QDP_all, {
	int c[4];
	QDP_get_coords(c, 0, i);
	printf(" %i %i %i %i", c[0], c[1], c[2], c[3]);
	for(int j=0; j<QLA_Nc; j++) {
	  QLA_ColorVector *v = QDP_site_ptr_readonly_V(x[0], i);
	  printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	  printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
	}
	printf("\n");
      });
    exit(0);
  }
#endif

  QDP_ColorVector *tsrc[2], *vtmp[2];
  QDP_ColorMatrix *tmat;
  for(int i=0; i<2; i++) tsrc[i] = QDP_create_V();
  for(int i=0; i<2; i++) vtmp[i] = QDP_create_V();
  tmat = QDP_create_M();

  for(int s=0; s<ns; s++) {
    QDP_M_eq_zero(mid[s], QDP_all);
    for(int i=-1; i<nterms; i++) {
      if(i+1<nterms) {
	int k = (i+1)%2;
	QDP_V_eq_V(tsrc[k], x[i+1], QDP_all);
	QDP_V_eq_sV(vtmp[k], tsrc[k], shifts[s], QDP_forward, QDP_all);
      }
      if(i>=0) {
	int k = i%2;
	QDP_M_eq_V_times_Va(tmat, tsrc[k], vtmp[k], QDP_all);
	QDP_discard_V(vtmp[k]);
	QDP_M_peq_r_times_M(mid[s], &eps[i], tmat, QDP_all);
      }
    }
  }

  for(int i=0; i<2; i++) QDP_destroy_V(tsrc[i]);
  for(int i=0; i<2; i++) QDP_destroy_V(vtmp[i]);
  QDP_destroy_M(tmat);

  info->final_sec = QOP_time() - dtime;
  info->final_flop = (54+36)*ns*nterms*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
#undef NC
}

#define NMTMP 4
#define NFTMP 3
#define NBTMP 3
static QDP_ColorMatrix *mtmp[NMTMP], *ftmp0[NFTMP], *ftmp[NFTMP][4],
  *btmp0[NBTMP], *btmp[NBTMP][4];
static int setcount=0;
#define set_temps() if(!setcount++) set_temps0()
#define free_temps() if(!--setcount) free_temps0()

#if QOP_Colors == 'N'
static int gnc;
#define NC gnc
#define SETNC(x) gnc = x
#else
#define SETNC(x) (void)0
#endif
#define SETNCF(x) SETNC(QDP_get_nc(x))

static void
set_temps0(void)
{
  for(int i=0; i<NMTMP; i++) mtmp[i] = QDP_create_M();
  for(int i=0; i<NFTMP; i++) {
    ftmp0[i] = QDP_create_M();
    for(int j=0; j<4; j++) ftmp[i][j] = QDP_create_M();
  }
  for(int i=0; i<NBTMP; i++) {
    btmp0[i] = QDP_create_M();
    for(int j=0; j<4; j++) btmp[i][j] = QDP_create_M();
  }
}

static void
free_temps0(void)
{
  for(int i=0; i<NMTMP; i++) QDP_destroy_M(mtmp[i]);
  for(int i=0; i<NFTMP; i++) {
    QDP_destroy_M(ftmp0[i]);
    for(int j=0; j<4; j++) QDP_destroy_M(ftmp[i][j]);
  }
  for(int i=0; i<NBTMP; i++) {
    QDP_destroy_M(btmp0[i]);
    for(int j=0; j<4; j++) QDP_destroy_M(btmp[i][j]);
  }
}

static void
staple(QDP_ColorMatrix *out, QDP_ColorMatrix *in0, QDP_ColorMatrix *link0, int mu, int nu)
{
#define link     ftmp0[0]
#define linkmu   ftmp[0][mu]
#define in       ftmp0[1]
#define innu     ftmp[1][nu]
#define linkin   mtmp[0]
#define back     btmp0[0]
#define backnu   btmp[0][nu]
#define linkinnu mtmp[1]

  QDP_M_eq_M(link, link0, QDP_all);
  QDP_M_eq_sM(linkmu, link, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(in, in0, QDP_all);
  QDP_M_eq_sM(innu, in, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(linkin, link, in, QDP_all);
  QDP_M_eq_M_times_M(back, linkin, linkmu, QDP_all);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(linkinnu, link, innu, QDP_all);
  QDP_discard_M(innu);
  QDP_M_peq_M_times_Ma(out, linkinnu, linkmu, QDP_all);
  QDP_discard_M(linkmu);
  QDP_M_peq_M(out, backnu, QDP_all);
  QDP_discard_M(backnu);
#define STAPLE_FLOPS (3*198+216+18)

#undef link
#undef linkmu
#undef in
#undef innu
#undef linkin
#undef back
#undef backnu
#undef linkinnu
}

static void
side_force(QDP_ColorMatrix *force, QDP_ColorMatrix *bot0, QDP_ColorMatrix *side0,
	   QDP_ColorMatrix *top0, int mu, int nu, QDP_ColorMatrix *stpl)
{
#define side     ftmp0[0]
#define sidemu   ftmp[0][mu]
#define top      ftmp0[1]
#define topnu    ftmp[1][nu]
#define bot      ftmp0[2]
#define botnu    ftmp[2][nu]
#define sidebot  mtmp[0]
#define sidetop  mtmp[1]
#define topnusidebot  btmp0[0]
#define fbmu          btmp[0][mu]
#define botnusidetop  btmp0[1]
#define fmbmu         btmp[1][mu]
#define sidebotsidemu btmp0[2]
#define stm           btmp[2][nu]
#define botnusidemu   mtmp[2]
#define botsidemu     mtmp[3]

  // force += bot * sidemu * topnu+
  // force -= bot-mu+ * side-mu * topnu-mu
  // -= top <-> bot
  // stpl += side * botnu * sidemu+
  // stpl += side-nu+ * bot-nu * sidemu-nu
  QDP_M_eq_M(side, side0, QDP_all);
  QDP_M_eq_sM(sidemu, side, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(top, top0, QDP_all);
  QDP_M_eq_sM(topnu, top, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_M(bot, bot0, QDP_all);
  QDP_M_eq_sM(botnu, bot, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(sidebot, side, bot, QDP_all);
  QDP_M_eq_Ma_times_M(sidetop, side, top, QDP_all);
  QDP_M_eq_Ma_times_M(topnusidebot, topnu, sidebot, QDP_all);
  QDP_M_eq_sM(fbmu, topnusidebot, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_M_eq_Ma_times_M(botnusidetop, botnu, sidetop, QDP_all);
  QDP_M_eq_sM(fmbmu, botnusidetop, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(sidebotsidemu, sidebot, sidemu, QDP_all);
  QDP_M_eq_sM(stm, sidebotsidemu, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_Ma(botnusidemu, botnu, sidemu, QDP_all);
  QDP_discard_M(botnu);
  QDP_M_peq_M_times_M(stpl, side, botnusidemu, QDP_all);
  //QDP_M_meq_M_times_Ma(force, top, botnusidemu, QDP_all);
  QDP_M_peq_M_times_Ma(force, top, botnusidemu, QDP_all);
  QDP_M_eq_M_times_M(botsidemu, bot, sidemu, QDP_all);
  QDP_discard_M(sidemu);
  QDP_M_peq_M_times_Ma(force, botsidemu, topnu, QDP_all);
  QDP_discard_M(topnu);
  //QDP_M_meq_Ma(force, fbmu, QDP_all);
  QDP_M_peq_Ma(force, fbmu, QDP_all);
  QDP_discard_M(fbmu);
  QDP_M_peq_Ma(force, fmbmu, QDP_all);
  QDP_discard_M(fmbmu);
  QDP_M_peq_M(stpl, stm, QDP_all);
  QDP_discard_M(stm);
#define SIDE_FORCE_FLOPS (7*198+3*216+3*18)

#undef side
#undef sidemu
#undef top
#undef topnu
#undef bot
#undef botnu
#undef sidebot
#undef sidetop
#undef topnusidebot
#undef fbmu
#undef botnusidetop
#undef fmbmu
#undef sidebotsidemu
#undef stm
#undef botnusidemu
#undef botsidemu
}

static void
QOP_fat_deriv(QOP_info_t *info, QDP_ColorMatrix *gauge[],
	      QDP_ColorMatrix *deriv[], QOP_asqtad_coeffs_t *coef,
	      QDP_ColorMatrix *mid[])
{
  double dtime = QOP_time();
  int nflops = 0;

  QLA_Real coef1 = coef->one_link;
  QLA_Real coef3 = coef->three_staple;
  QLA_Real coef5 = coef->five_staple;
  QLA_Real coef7 = coef->seven_staple;
  QLA_Real coefL = coef->lepage;
  coef1 -= 6*coefL;
  int have5 = coef5 || coef7 || coefL;
  int have3 = coef3 || have5;
  QDP_ColorMatrix *stpl3[4];
  QDP_ColorMatrix *mid5[4];
  QDP_ColorMatrix *stpl5=NULL, *mid3=NULL;
  if(have3) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
	stpl3[mu] = QDP_create_M();
	mid5[mu] = QDP_create_M();
      }
    }
    stpl5 = QDP_create_M();
    mid3 = QDP_create_M();
    set_temps();
  }

  for(int sig=0; sig<4; sig++) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
	if(mu==sig) continue;
	QDP_M_eq_zero(stpl3[mu], QDP_all);
	staple(stpl3[mu], gauge[sig], gauge[mu], sig, mu);
	QDP_M_eq_zero(mid5[mu], QDP_all);
	nflops += STAPLE_FLOPS;
      }
      for(int rho=0; rho<4; rho++) {
	if(rho==sig) continue;
	QDP_M_eq_zero(stpl5, QDP_all);

	if(coef7) {
	  for(int mu=0; mu<4; mu++) {
	    if(mu==sig||mu==rho) continue;
	    for(int nu=0; nu<4; nu++) {
	      if(nu==sig||nu==rho||nu==mu) continue;
	      staple(stpl5, stpl3[mu], gauge[nu], sig, nu);
	      nflops += STAPLE_FLOPS;
	    }
	  }
	  QDP_M_eq_r_times_M(stpl5, &coef7, stpl5, QDP_all);
	  nflops += 18;
	}

	QDP_M_eq_zero(mid3, QDP_all);
	if(coefL) {
	  QDP_M_peq_r_times_M(stpl5, &coefL, stpl3[rho], QDP_all);
	  nflops += 36;
	}
	if(coefL || coef7) {
	  side_force(deriv[rho], mid[sig], gauge[rho], stpl5, sig, rho, mid3);
	  nflops += SIDE_FORCE_FLOPS;
	}
	if(coefL) {
	  QDP_M_peq_r_times_M(mid5[rho], &coefL, mid3, QDP_all);
	  nflops += 36;
	}

	//QDP_M_eqm_r_times_M(mid3, &coef7, mid3, QDP_all);
	QDP_M_eq_r_times_M(mid3, &coef7, mid3, QDP_all);
	QDP_M_peq_r_times_M(mid3, &coef5, mid[sig], QDP_all);
	nflops += 18+36;
	for(int mu=0; mu<4; mu++) {
	  if(mu==sig||mu==rho) continue;
	  for(int nu=0; nu<4; nu++) {
	    if(nu==sig||nu==rho||nu==mu) continue;
	    side_force(deriv[mu], mid3, gauge[mu], stpl3[nu], sig, mu, mid5[nu]);
	    nflops += SIDE_FORCE_FLOPS;
	  }
	}
      }
    }

    if(have3) {
      QDP_M_eq_zero(mid3, QDP_all);
      for(int mu=0; mu<4; mu++) {
	if(mu==sig) continue;
	if(have5) {
	  //QDP_M_eq_r_times_M_minus_M(stpl5, &coef3, mid[sig], mid5[mu], QDP_all);
	  QDP_M_eq_r_times_M_plus_M(stpl5, &coef3, mid[sig], mid5[mu], QDP_all);
	  nflops += 36;
	} else {
	  QDP_M_eq_r_times_M(stpl5, &coef3, mid[sig], QDP_all);
	  nflops += 18;
	}
	side_force(deriv[mu], stpl5, gauge[mu], gauge[sig], sig, mu, mid3);
	nflops += SIDE_FORCE_FLOPS;
      }
      //QDP_M_meq_M(deriv[sig], mid3, QDP_all);
      QDP_M_peq_M(deriv[sig], mid3, QDP_all);
      nflops += 18;
    }
    if(coef1) {
      QDP_M_peq_r_times_M(deriv[sig], &coef1, mid[sig], QDP_all);
      nflops += 36;
    }
  }

#if 1
  if(coefL) {
    // fix up Lepage term
    QLA_Real fixL = -coefL;
    for(int mu=0; mu<4; mu++) {
      QDP_M_eq_zero(ftmp0[0], QDP_all);
      for(int nu=0; nu<4; nu++) {
	if(nu==mu) continue;
	QDP_M_eq_Ma_times_M(btmp0[0], mid[nu], gauge[nu], QDP_all);
	QDP_M_eq_sM(btmp[0][nu], btmp0[0], QDP_neighbor[nu], QDP_backward, QDP_all);
	QDP_M_eq_M_times_Ma(stpl5, mid[nu], gauge[nu], QDP_all);
	QDP_M_meq_M(stpl5, btmp[0][nu], QDP_all);
	QDP_discard_M(btmp[0][nu]);
	QDP_M_peq_M(ftmp0[0], stpl5, QDP_all);
	QDP_M_peq_Ma(ftmp0[0], stpl5, QDP_all);
      }
      QDP_M_eq_sM(ftmp[0][mu], ftmp0[0], QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_M_times_M(stpl5, ftmp0[0], gauge[mu], QDP_all);
      QDP_M_meq_M_times_M(stpl5, gauge[mu], ftmp[0][mu], QDP_all);
      QDP_discard_M(ftmp[0][mu]);
      QDP_M_peq_r_times_M(deriv[mu], &fixL, stpl5, QDP_all);
    }
    nflops += 4*(3*(2*198+3*18)+198+216+36);
  }
#endif

  if(have3) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
	QDP_destroy_M(stpl3[mu]);
	QDP_destroy_M(mid5[mu]);
      }
    }
    QDP_destroy_M(stpl5);
    QDP_destroy_M(mid3);
    free_temps();
  }

  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}

void
QOP_asqtad_deriv(QOP_info_t *info, QDP_ColorMatrix *gauge[],
		 QDP_ColorMatrix *deriv[], QOP_asqtad_coeffs_t *coef,
		 QDP_ColorMatrix *mid_fat[],
		 QDP_ColorMatrix *mid_naik[])
{
  SETNCF(deriv[0]);
  double dtime = QOP_time();
  int nflops = 0;
  QOP_info_t tinfo;
  set_temps();

#if 0
  {
    QLA_Real nrm2=0, t;
    if(mid_fat) {
      for(int i=0; i<4; i++) {
	QDP_r_eq_norm2_M(&t, mid_fat[i], QDP_all);
	nrm2 += t;
      }
      QOP_printf0("norm2 mid_fat: %g\n", nrm2);
    }
    if(mid_naik) {
      nrm2 = 0;
      for(int i=0; i<4; i++) {
	QDP_r_eq_norm2_M(&t, mid_naik[i], QDP_all);
	nrm2 += t;
      }
      QOP_printf0("norm2 mid_naik: %g\n", nrm2);
    }
  }
#endif

  // fat
  if(mid_fat) {
    QOP_fat_deriv(&tinfo, gauge, deriv, coef, mid_fat);
    nflops += tinfo.final_flop;
  }

  // Naik
  QLA_Real coefN = coef->naik;
  if(coefN && mid_naik) {
#define U       gauge[mu]
#define mid     mid_naik[mu]
#define Uf      ftmp0[0]
#define Umu     ftmp[0][mu]
#define Umid    btmp0[0]
#define Umidbmu btmp[0][mu]
#define UmuU    ftmp0[1]
#define UmuUs   ftmp[1][mu]
#define f3b     btmp0[1]
#define f3      btmp[1][mu]
#define f       mtmp[0]
    for(int mu=0; mu<4; mu++) {
      // force += mid * Umumu+ * Umu+
      // force -= U-mu+ * mid-mu * Umu+
      // force += U-mu+ * U-mu-mu+ * mid-mu-mu
      QDP_M_eq_M(Uf, U, QDP_all);
      QDP_M_eq_sM(Umu, Uf, QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_Ma_times_M(Umid, Uf, mid, QDP_all);
      QDP_M_eq_sM(Umidbmu, Umid, QDP_neighbor[mu], QDP_backward, QDP_all);
      QDP_M_eq_Ma_times_Ma(UmuU, Umu, Uf, QDP_all);
      QDP_M_eq_sM(UmuUs, UmuU, QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_Ma_times_M(f3b, Uf, Umidbmu, QDP_all);
      QDP_M_eq_sM(f3, f3b, QDP_neighbor[mu], QDP_backward, QDP_all);
      QDP_M_eq_M_times_M(f, mid, UmuUs, QDP_all);
      QDP_discard_M(UmuUs);
      //QDP_M_meq_M_times_Ma(f, Umidbmu, Umu, QDP_all);
      QDP_M_peq_M_times_Ma(f, Umidbmu, Umu, QDP_all);
      QDP_discard_M(Umidbmu);
      QDP_discard_M(Umu);
      QDP_M_peq_M(f, f3, QDP_all);
      QDP_discard_M(f3);
      QDP_M_peq_r_times_M(deriv[mu], &coefN, f, QDP_all);
    }
#undef U
#undef mid
#undef Uf
#undef Umu
#undef Umid
#undef Umidbmu
#undef UmuU
#undef UmuUs
#undef f3b
#undef f3
#undef f
    nflops += 4*(4*198+216+18+36)*QDP_sites_on_node;
  }
  free_temps();

#if 0
  {
    QLA_Real nrm2=0, t;
    for(int i=0; i<4; i++) {
      QDP_r_eq_norm2_M(&t, deriv[i], QDP_all);
      nrm2 += t;
    }
    QOP_printf0("norm2 deriv: %g\n", nrm2);
  }
#endif

  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops;
  info->status = QOP_SUCCESS;
}

void
QOP_asqtad_force_multi_fnmat3_qdp(QOP_info_t *info, QOP_GaugeField *gauge,
				  QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
				  QDP_ColorMatrix *mid_fat[],
				  QDP_ColorMatrix *mid_naik[])
{
  SETNCF(force[0]);
  double dtime = QOP_time();
  QOP_info_t tinfo;
  QDP_ColorMatrix *cm, *deriv[4], *gcm[4];
  cm = QDP_create_M();
  for(int mu=0; mu<4; mu++) {
    deriv[mu] = QDP_create_M();
    QDP_M_eq_zero(deriv[mu], QDP_all);
    gcm[mu] = gauge->links[mu];
  }
  QOP_asqtad_deriv(&tinfo, gcm, deriv, coef, mid_fat, mid_naik);
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M_times_Ma(cm, gauge->links[mu], deriv[mu], QDP_all);
    QDP_M_eq_antiherm_M(deriv[mu], cm, QDP_all);
    //QDP_M_peq_M(force[mu], deriv[mu], QDP_even);
    //QDP_M_meq_M(force[mu], deriv[mu], QDP_odd);
    QDP_M_peq_M(force[mu], deriv[mu], QDP_all);
    QDP_destroy_M(deriv[mu]);
  }
  QDP_destroy_M(cm);

  info->final_sec = QOP_time() - dtime;
  info->final_flop = tinfo.final_flop;
  info->final_flop += 4*(198+24+18)*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}
