//#define DO_TRACE
#include <qop_internal.h>
//#undef QOP_printf0
//#define QOP_printf0(...)
#define QOP_trace(...)

#define OPP_DIR(i) (7-i)
#define GOES_FORWARDS(dir) (dir<4)
#define GOES_BACKWARDS(dir) (dir>=4)

static void u_shift_mat(QDP_ColorMatrix *src, QDP_ColorMatrix *dest,
			int dir, QDP_ColorMatrix *tmpmat);

static void shift_mat(QDP_ColorMatrix *src, QDP_ColorMatrix *dest, int dir);
static void shift_mid(QDP_ColorMatrix *src[], QDP_ColorMatrix *dest[], int dir);

static void add_forces_to_mom(QDP_ColorMatrix *back, QDP_ColorMatrix *forw,
			      QDP_ColorMatrix *mid, int dir, REAL coeff);

static void side_link_forces(int mu, int nu, REAL coeff,
			     QDP_ColorMatrix *Path, QDP_ColorMatrix *Path_nu,
			     QDP_ColorMatrix *Path_mu, QDP_ColorMatrix *Path_numu,
			     QDP_ColorMatrix *mid, QDP_ColorMatrix *mid_mu);

//static su3_matrix *backwardlink[4];
static QDP_ColorMatrix *tempmom_qdp[4];

//static QDP_ColorMatrix *bcklink[4];
//static QDP_ColorMatrix *fwdlink[4];
//static QDP_ColorMatrix *hw_qdp[8][2];

static QDP_ColorMatrix *fblink[8];
static QDP_Shift fbshift[8];
static QDP_ShiftDir fbshiftdir[8];
static QDP_ColorMatrix *tm;

//extern QOP_asqtad_force_t QOP_asqtad_ff;

void
QOPPC(asqtad_force_multi_fnmat2)(QOP_info_t *info, QOP_GaugeField *gauge,
				 QOP_Force *force, QOP_asqtad_coeffs_t *coef,
				 QDP_ColorMatrix *mid_fat[], QDP_ColorMatrix *mid_naik[]);

void 
QOPPC(asqtad_force_multi_fnmat)(QOP_info_t *info, QOP_GaugeField *Gauge,
				QOP_Force *Force, QOP_asqtad_coeffs_t *asq_coeff,
				REAL eps[], QOP_ColorVector *x[], int nterms)
{
  ASQTAD_FORCE_BEGIN;

  double nflop1 = 253935;
  double nflop2 = 433968;
  double nflop = nflop1 + (nflop2-nflop1)*(3-1);
  double dtime;
  dtime = -QOP_time();

  QDP_ColorVector *vec_tmp[2];
  QDP_ColorMatrix *mid_fat[4], *mid_naik[4], *tmat;
  for(int i=0; i<2; i++) vec_tmp[i] = QDP_create_V();
  tmat = QDP_create_M();

  if(!QOP_asqtad.inited) QOP_asqtad_invert_init();

  // gather vectors for fat link term
  for(int mu=0; mu<4; mu++) {
    mid_fat[mu] = QDP_create_M();
    QDP_M_eq_zero(mid_fat[mu], QDP_all);
    for(int i=0; i<=nterms; i++) {
      if(i<nterms) {
	QDP_V_eq_sV(vec_tmp[i%2], x[i]->cv, QDP_neighbor[mu], QDP_forward, QDP_all);
      }
      if(i>0) {
	QDP_M_eq_V_times_Va(tmat, x[i-1]->cv, vec_tmp[(i-1)%2], QDP_all);
	QDP_discard_V(vec_tmp[(i-1)%2]);
	QDP_M_peq_r_times_M(mid_fat[mu], &eps[i-1], tmat, QDP_all);
      }
    }
  }
  // gather vectors for naik term
  for(int mu=0; mu<4; mu++) {
    mid_naik[mu] = QDP_create_M();
    QDP_M_eq_zero(mid_naik[mu], QDP_all);
    for(int i=0; i<=nterms; i++) {
      if(i<nterms) {
	QDP_V_eq_sV(vec_tmp[i%2], x[i]->cv, QOP_common.neighbor3[mu], QDP_forward, QDP_all);
      }
      if(i>0) {
	QDP_M_eq_V_times_Va(tmat, x[i-1]->cv, vec_tmp[(i-1)%2], QDP_all);
	QDP_discard_V(vec_tmp[(i-1)%2]);
	QDP_M_peq_r_times_M(mid_naik[mu], &eps[i-1], tmat, QDP_all);
      }
    }
  }

  for(int i=0; i<2; i++) QDP_destroy_V(vec_tmp[i]);
  QDP_destroy_M(tmat);

  QOPPC(asqtad_force_multi_fnmat2)(info,Gauge,Force,asq_coeff,mid_fat,mid_naik);

  for(int mu=0; mu<4; mu++) {
    QDP_destroy_M(mid_fat[mu]);
    QDP_destroy_M(mid_naik[mu]);
  }

  dtime += QOP_time();

  info->final_sec = dtime;
  info->final_flop = nflop*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  ASQTAD_FORCE_END;
}

void
QOPPC(asqtad_force_multi_fnmat2)(QOP_info_t *info, QOP_GaugeField *gauge,
				 QOP_Force *force, QOP_asqtad_coeffs_t *coef,
				 QDP_ColorMatrix *mid_fat[], QDP_ColorMatrix *mid_naik[])
{
  REAL OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt;
  REAL mNaik, mLepage, mFiveSt, mThreeSt, mSevenSt;

  QDP_ColorMatrix *P3[8];

  QDP_ColorMatrix *P5[8];
  QDP_ColorMatrix *P5tmp[8][8];
  QDP_ColorMatrix *P5s[4];
  QDP_ColorMatrix *P5tmps[4][8];

  QDP_ColorMatrix *Pmu;
  QDP_ColorMatrix *Pmutmp[8];
  QDP_ColorMatrix *Pnumu;
  QDP_ColorMatrix *Pnumutmp[8];
  QDP_ColorMatrix *Prhonumu;
  QDP_ColorMatrix *Prhonumutmp[8];
  QDP_ColorMatrix *P7;
  QDP_ColorMatrix *P7tmp[8];
  QDP_ColorMatrix *P7rho;

  QDP_ColorMatrix *mid[8];
  QDP_ColorMatrix *midmu[8];
  QDP_ColorMatrix *midnumu[8];
  QDP_ColorMatrix *midrhonumu[8];

  /* setup parallel transport */
  tm = QDP_create_M(); // global, also used later
  for(int i=0; i<QOP_common.ndim; i++) {
    fbshift[i] = QDP_neighbor[i];
    fbshiftdir[i] = QDP_forward;
    fblink[i] = gauge->links[i];
    fbshift[OPP_DIR(i)] = QDP_neighbor[i];
    fbshiftdir[OPP_DIR(i)] = QDP_backward;
    fblink[OPP_DIR(i)] = QDP_create_M();
    QDP_M_eq_sM(tm, fblink[i], QDP_neighbor[i], QDP_backward, QDP_all);
    QDP_M_eq_Ma(fblink[OPP_DIR(i)], tm, QDP_all);
  }

  /* Allocate temporary vectors */
  Pmu = QDP_create_M();
  Pnumu = QDP_create_M();
  Prhonumu = QDP_create_M();
  P7 = QDP_create_M();
  P7rho = QDP_create_M();
  for(int dir=0; dir<8; dir++) {
    Pmutmp[dir] = QDP_create_M();
    Pnumutmp[dir] = QDP_create_M();
    Prhonumutmp[dir] = QDP_create_M();
    P7tmp[dir] = QDP_create_M();
    midmu[dir] = QDP_create_M();
    midnumu[dir] = QDP_create_M();
    midrhonumu[dir] = QDP_create_M();
  }
#define midtmp Pmutmp
#define midmutmp Pnumutmp
#define midnumutmp Prhonumutmp
#if 1
  for(int mu=0; mu<4; mu++) {
    P5s[mu] = QDP_create_M();
    for(int dir=0; dir<8; dir++) {
      P5tmps[mu][dir] = QDP_create_M();
    }
  }
#else
  for(int mu=0; mu<8; mu++) {
    P5[mu] = QDP_create_M();
    for(int dir=0; dir<8; dir++) {
      P5tmp[mu][dir] = QDP_create_M();
      //printf("%p %p\n", P5tmp[mu][dir], &(P5tmp[mu][dir])); fflush(stdout);
      if(P5tmp[mu][dir]==NULL) {
	fprintf(stderr, "error: can't create V\n");
	QDP_abort();
      }
    }
  }
#endif
  //printf("%p\n", P5tmp[0][4][0]); fflush(stdout);

  for(int mu=0; mu<8; mu++) {
    P3[mu] = QDP_create_M();
    //P5[mu] = QDP_create_M();
  }

  for(int mu=0; mu<4; mu++) {
    tempmom_qdp[mu] = force->force[mu];
    QDP_M_eqm_M(tempmom_qdp[mu], tempmom_qdp[mu], QDP_odd);
  }

  for(int mu=0; mu<4; mu++) {
    mid[mu] = mid_fat[mu];
    mid[OPP_DIR(mu)] = QDP_create_M();
    QDP_M_eq_sM(tm, mid[mu], QDP_neighbor[mu], QDP_backward, QDP_all);
    QDP_M_eq_Ma(mid[OPP_DIR(mu)], tm, QDP_all);
  }

  /* Load path coefficients from table */
  /* epsilon is now absorbed in mid matrices */
  OneLink = coef->one_link     ;
  Naik    = coef->naik         ; mNaik    = -Naik;
  ThreeSt = coef->three_staple ; mThreeSt = -ThreeSt;
  FiveSt  = coef->five_staple  ; mFiveSt  = -FiveSt;
  SevenSt = coef->seven_staple ; mSevenSt = -SevenSt;
  Lepage  = coef->lepage       ; mLepage  = -Lepage;

  /* *************************************** */

  QOP_trace("start force loop\n");
  for(int mu=0; mu<8; mu++) {
    //u_shift_mat(unit, Pmu, OPP_DIR(mu), unittmp[OPP_DIR(mu)]);
    QDP_M_eq_M(Pmu, fblink[OPP_DIR(mu)], QDP_all);
    //shift_mid(mid, midmu, OPP_DIR(mu));

    for(int sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) ) {
	u_shift_mat(Pmu, P3[sig], sig, Pmutmp[sig]);
	//shift_mat(mid[sig], midmu[sig], OPP_DIR(mu));
	u_shift_mat(mid[sig], midmu[sig], OPP_DIR(mu), midtmp[sig]);
	if(GOES_FORWARDS(sig)) {
	  /* Add the force F_sig[x+mu]:         x--+             *
	   *                                   |   |             *
	   *                                   o   o             *
	   * the 1 link in the path: - (numbering starts form 0) */
	  //add_forces_to_mom(P3[sig], Pmu, midmu[sig], sig, mThreeSt);
	  add_forces_to_mom(P3[sig], midmu[sig], NULL, sig, mThreeSt);
	}
      }

    for(int nu=0; nu<8; nu++) if( (nu!=mu)&&(nu!=OPP_DIR(mu)) ) {
	int nP5 = 0;
	u_shift_mat(Pmu, Pnumu, OPP_DIR(nu), Pmutmp[OPP_DIR(nu)]);
	//shift_mid(midmu, midnumu, OPP_DIR(nu));

	for(int sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
					 (sig!=nu)&&(sig!=OPP_DIR(nu)) ) {
	    P5[sig] = P5s[nP5];
	    for(int dir=0; dir<8; dir++) P5tmp[sig][dir] = P5tmps[nP5][dir];
	    nP5++;
	    u_shift_mat(Pnumu, P5[sig], sig, Pnumutmp[sig]);
	    //shift_mat(midmu[sig], midnumu[sig], OPP_DIR(nu));
	    u_shift_mat(midmu[sig], midnumu[sig], OPP_DIR(nu), midmutmp[sig]);
	    if(GOES_FORWARDS(sig)) {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      //add_forces_to_mom(P5[sig], Pnumu, midnumu[sig], sig, FiveSt);
	      add_forces_to_mom(P5[sig], midnumu[sig], NULL, sig, FiveSt);
	    }
	  }

	for(int rho=0; rho<8; rho++) if( (rho!=mu)&&(rho!=OPP_DIR(mu)) &&
					 (rho!=nu)&&(rho!=OPP_DIR(nu)) ) {
	    u_shift_mat(Pnumu, Prhonumu, OPP_DIR(rho),
			      Pnumutmp[OPP_DIR(rho)]);
	    //shift_mid(midnumu, midrhonumu, OPP_DIR(rho));
	    for(int sig=0; sig<8; sig++) if( (sig!=mu )&&(sig!=OPP_DIR(mu )) &&
					     (sig!=nu )&&(sig!=OPP_DIR(nu )) &&
					     (sig!=rho)&&(sig!=OPP_DIR(rho)) ) {
		/* Length 7 paths */
		u_shift_mat(Prhonumu, P7, sig, Prhonumutmp[sig]);
		//shift_mat(midnumu[sig], midrhonumu[sig], OPP_DIR(rho));
		u_shift_mat(midnumu[sig], midrhonumu[sig], OPP_DIR(rho), midnumutmp[sig]);
		if(GOES_FORWARDS(sig)) {
		  /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 3 link in the path: - (numbering starts form 0) */
		  //add_forces_to_mom(P7, Prhonumu, midrhonumu[sig], sig, mSevenSt);
		  add_forces_to_mom(P7, midrhonumu[sig], NULL, sig, mSevenSt);
		}
		/* Add the force F_rho the 2(4) link in the path: +     */
		u_shift_mat(P7, P7rho, rho, P7tmp[rho]);
		//side_link_forces(rho, sig, SevenSt, Pnumu, P7, Prhonumu, P7rho,
		//midnumu[sig], midrhonumu[sig]);
		side_link_forces(rho, sig, SevenSt, midnumu[sig], P7, midrhonumu[sig], P7rho,
				 NULL, NULL);
		/* Add the P7rho vector to P5 */
		{
		  QLA_Real coeff = 0;
		  if(FiveSt!=0) coeff = SevenSt/FiveSt;
		  QDP_M_peq_r_times_M(P5[sig], &coeff, P7rho, QDP_all);
		}
	      } /* sig */
	  } /* rho */

#define P5nu P7
	for(int sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
					 (sig!=nu)&&(sig!=OPP_DIR(nu)) ) {
	    /* Length 5 paths */
	    /* Add the force F_nu the 1(3) link in the path: -     */
	    u_shift_mat(P5[sig], P5nu, nu, P5tmp[sig][nu]);
	    //side_link_forces(nu, sig, mFiveSt, Pmu, P5[sig], Pnumu, P5nu,
	    //midmu[sig], midnumu[sig]);
	    side_link_forces(nu, sig, mFiveSt, midmu[sig], P5[sig], midnumu[sig], P5nu,
			     NULL, NULL);
	    /* Add the P5nu vector to P3 */
	    {
	      QLA_Real coeff = 0;
	      if(ThreeSt!=0) coeff = FiveSt/ThreeSt; 
	      QDP_M_peq_r_times_M(P3[sig], &coeff, P5nu, QDP_all);
	    }
	  } /* sig */
      } /* nu */

#define Pmumu Pnumu
#define Pmumutmp Pnumutmp
#define P5sig Prhonumu
#define P5sigtmp Prhonumutmp
#define P3mu P7
#define Popmu P7
#define Pmumumu P7
    /* Now the Lepage term... It is the same as 5-link paths with
       nu=mu and FiveSt=Lepage. */
    u_shift_mat(Pmu, Pmumu, OPP_DIR(mu), Pmutmp[OPP_DIR(mu)]);
    //shift_mid(midmu, midnumu, OPP_DIR(mu));

    for(int sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) ) {
	u_shift_mat(Pmumu, P5sig, sig, Pmumutmp[sig]);
	//shift_mat(midmu[sig], midnumu[sig], OPP_DIR(mu));
	u_shift_mat(midmu[sig], midnumu[sig], OPP_DIR(mu), midmutmp[sig]);
	if(GOES_FORWARDS(sig)) {
	  /* Add the force F_sig[x+mu+nu]:      x--+             *
	   *                                   |   |             *
	   *                                   o   o             *
	   * the 2 link in the path: + (numbering starts form 0) */
	  //add_forces_to_mom(P5sig, Pmumu, midnumu[sig], sig, Lepage);
	  add_forces_to_mom(P5sig, midnumu[sig], NULL, sig, Lepage);
	}
	/* Add the force F_nu the 1(3) link in the path: -     */
	u_shift_mat(P5sig, P5nu, mu, P5sigtmp[mu]);
	//side_link_forces(mu, sig, mLepage, Pmu, P5sig, Pmumu, P5nu,
	//midmu[sig], midnumu[sig]);
	side_link_forces(mu, sig, mLepage, midmu[sig], P5sig, midnumu[sig], P5nu,
			 NULL, NULL);
	/* Add the P5nu vector to P3 */
	{
	  QLA_Real coeff = 0;
	  if(ThreeSt!=0) coeff = Lepage/ThreeSt;
	  QDP_M_peq_r_times_M(P3[sig], &coeff, P5nu, QDP_all);
	}
	/* Length 3 paths (Not the Naik term) */
	/* Add the force F_mu the 0(2) link in the path: +     */
	if(GOES_FORWARDS(mu)) {
	  QDP_M_eq_M(P5sig, P3[sig], QDP_all);
	  u_shift_mat(P5sig, P3mu, mu, P5sigtmp[mu]);
	  /* The above shift is not needed if mu is backwards */
	}
	//side_link_forces(mu, sig, ThreeSt, NULL, P3[sig], Pmu, P3mu,
	//mid[sig], midmu[sig]);
	side_link_forces(mu, sig, ThreeSt, mid[sig], P3[sig], midmu[sig], P3mu,
			 NULL, NULL);
      }

    /* Finally the OneLink and the Naik term */
    if(GOES_BACKWARDS(mu)) {
      /* Do only the forward terms in the Dslash */
      /* Because I have shifted with OPP_DIR(mu) Pmu is a forward shift. */
      /* The one link */
      //add_forces_to_mom(Pmu, NULL, mid[OPP_DIR(mu)], OPP_DIR(mu), OneLink);
      add_forces_to_mom(Pmu, mid[OPP_DIR(mu)], NULL, OPP_DIR(mu), OneLink);

      /* For the same reason Pmumu is the forward double link */
      /* Popmu is a backward shift */
      //u_shift_mat(unit, Popmu, mu, unittmp[mu]);
      //QDP_M_eq_M(Popmu, fblink[mu], QDP_all);
      QDP_M_eq_M(P5sig, mid_naik[OPP_DIR(mu)], QDP_all);
      u_shift_mat(P5sig, Popmu, mu, P5sigtmp[mu]);
      /* The Naik */
      /* link no 1: - */
      //add_forces_to_mom(Pmumu, Popmu, OPP_DIR(mu), mNaik, mid_naik[OPP_DIR(mu)]);
      add_forces_to_mom(Pmumu, Popmu, NULL, OPP_DIR(mu), mNaik);

      QDP_M_eq_M(P5sig, Popmu, QDP_all);
      u_shift_mat(P5sig, Popmu, mu, P5sigtmp[mu]);
      add_forces_to_mom(Pmu, Popmu, NULL, OPP_DIR(mu), Naik);

      /* Pmumumu can overwrite Popmu which is no longer needed */
      u_shift_mat(Pmumu, Pmumumu, OPP_DIR(mu), Pmumutmp[OPP_DIR(mu)]);
      /* link no 0: + */
      add_forces_to_mom(Pmumumu, NULL, mid_naik[OPP_DIR(mu)], OPP_DIR(mu), Naik);

    } else {
      /* The rest of the Naik terms */
      //u_shift_mat(unit, Popmu, mu, unittmp[mu]);
      //QDP_M_eq_M(Popmu, fblink[mu], QDP_all);
      //QDP_M_eq_Ma(P5sig, mid_naik[mu], QDP_all);
      //u_shift_mat(P5sig, Popmu, mu, P5sigtmp[mu]);
      /* link no 2: + */
      /* Pmumu is double backward shift */
      //add_forces_to_mom(Popmu, Pmumu, NULL, mu, Naik);
    }
    /* Here we have to do together the Naik term and the one link term */

  }/* mu */

  TRACE;
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M(tm, tempmom_qdp[mu], QDP_even);
    QDP_M_eqm_M(tm, tempmom_qdp[mu], QDP_odd);
    QDP_M_eq_antiherm_M(tempmom_qdp[mu], tm, QDP_all);
  }

  TRACE;
  /* Free temporary vectors */
  QDP_destroy_M(Pmu);
  QDP_destroy_M(Pnumu);
  QDP_destroy_M(Prhonumu);
  QDP_destroy_M(P7);
  QDP_destroy_M(P7rho);
  TRACE;
  for(int dir=0; dir<8; dir++) {
    QDP_destroy_M(Pmutmp[dir]);
    QDP_destroy_M(Pnumutmp[dir]);
    QDP_destroy_M(Prhonumutmp[dir]);
    QDP_destroy_M(P7tmp[dir]);
    QDP_destroy_M(midmu[dir]);
    QDP_destroy_M(midnumu[dir]);
    QDP_destroy_M(midrhonumu[dir]);
  }
  TRACE;
  for(int mu=0; mu<4; mu++) {
    QDP_destroy_M(P5s[mu]);
    QDP_destroy_M(mid[OPP_DIR(mu)]);
    for(int dir=0; dir<8; dir++) {
      QDP_destroy_M(P5tmps[mu][dir]);
    }
  }

  TRACE;
  for(int mu=0; mu<8; mu++) {
    QDP_destroy_M(P3[mu]);
  }

  TRACE;
  QDP_destroy_M(tm);

  TRACE;
  for(int i=4; i<8; i++) {
    QDP_destroy_M(fblink[i]);
  }
}

#undef Pmu          
#undef Pnumu        
#undef Prhonumu     
#undef P7           
#undef P7rho        
#undef P7rhonu      
//#undef P5
#undef P3           
#undef P5nu         
#undef P3mu         
#undef Popmu        
#undef Pmumumu      

static void
u_shift_mat(QDP_ColorMatrix *src, QDP_ColorMatrix *dest,
	    int dir, QDP_ColorMatrix *tmpmat)
{
  QDP_M_eq_sM(tmpmat, src, fbshift[dir], fbshiftdir[dir], QDP_all);
  QDP_M_eq_M_times_M(dest, fblink[dir], tmpmat, QDP_all);
  QDP_discard_M(tmpmat);
}

static void
shift_mat(QDP_ColorMatrix *src, QDP_ColorMatrix *dest, int dir)
{
  //QDP_M_eq_sM(dest[i], src[i], fbshift[dir], fbshiftdir[dir], QDP_all);
  QDP_M_eq_sM(tm, src, fbshift[dir], fbshiftdir[dir], QDP_all);
  QDP_M_eq_M(dest, tm, QDP_all);
}

static void
shift_mid(QDP_ColorMatrix *src[], QDP_ColorMatrix *dest[], int dir)
{
  for(int i=0; i<8; i++) {
    //QDP_M_eq_sM(dest[i], src[i], fbshift[dir], fbshiftdir[dir], QDP_all);
    QDP_M_eq_sM(tm, src[i], fbshift[dir], fbshiftdir[dir], QDP_all);
    QDP_M_eq_M(dest[i], tm, QDP_all);
  }
}

/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
// mid should have followed the same path as forw
static void
add_forces_to_mom(QDP_ColorMatrix *back, QDP_ColorMatrix *forw, 
		  QDP_ColorMatrix *mid, int dir, REAL coeff)
{
  REAL tmp_coeff;
  if(GOES_BACKWARDS(dir)) {
    dir = OPP_DIR(dir); 
    tmp_coeff = -coeff;
  } else {
    tmp_coeff = coeff;
  }
  if(back==NULL) {
    QDP_M_eq_Ma_times_Ma(tm, mid, forw, QDP_all);
    QDP_M_peq_r_times_M(tempmom_qdp[dir], &tmp_coeff, tm, QDP_all);
  } else if(forw==NULL) {
    QDP_M_eq_r_times_M(tm, &tmp_coeff, mid, QDP_all);
    QDP_M_peq_M_times_Ma(tempmom_qdp[dir], back, tm, QDP_all);
  } else if(mid==NULL) {
    QDP_M_eq_r_times_M(tm, &tmp_coeff, forw, QDP_all);
    QDP_M_peq_M_times_Ma(tempmom_qdp[dir], back, tm, QDP_all);
  } else {
    QDP_M_eq_M_times_M(tm, forw, mid, QDP_all);
    QDP_M_eq_r_times_M(tm, &tmp_coeff, tm, QDP_all);
    QDP_M_peq_M_times_Ma(tempmom_qdp[dir], back, tm, QDP_all);
  }
}

/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */
static void
side_link_forces(int mu, int nu, REAL coeff,
		 QDP_ColorMatrix *Path, QDP_ColorMatrix *Path_nu,
		 QDP_ColorMatrix *Path_mu, QDP_ColorMatrix *Path_numu,
		 QDP_ColorMatrix *mid, QDP_ColorMatrix *mid_mu)
{
  REAL m_coeff = -coeff;

  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu))
	add_forces_to_mom(Path_numu, Path, mid, mu, coeff);
      else
	add_forces_to_mom(Path_numu, Path, mid, mu, m_coeff);
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      if(GOES_FORWARDS(nu))
	add_forces_to_mom(Path_nu, Path_mu, mid_mu, mu, m_coeff);
      else
	add_forces_to_mom(Path_nu, Path_mu, mid_mu, mu, coeff);
    }
}


/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
