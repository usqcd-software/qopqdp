/* gauge action for Symanzik improved 1x1 + 1x2 + 1x1x1 */

#include <qop_internal.h>

#define EQMTM (QLA_Nc*QLA_Nc*(8*QLA_Nc-2))
#define PEQMTM (QLA_Nc*QLA_Nc*(8*QLA_Nc))

void 
QOP_symanzik_1loop_gauge_action(QOP_info_t *info, QOP_GaugeField *gauge,
				REAL *acts, REAL *actt, QOP_gauge_coeffs_t *coeffs)
{
#define NC QDP_get_nc(gauge->links[0])
  double dtime = QOP_time();
  double nflops = 0;

  QLA_Real fac = 1./QLA_Nc;
  QLA_Real plaq = fac*coeffs->plaquette;
  QLA_Real rect = fac*coeffs->rectangle;
  QLA_Real pgm  = fac*coeffs->parallelogram;
  QLA_Real adpl = fac*fac*coeffs->adjoint_plaquette;
  QLA_Real plaqs=0, plaqt=0;
  QLA_Real rects=0, rectt=0;
  QLA_Real pgms=0, pgmt=0;
  QLA_Real adpls=0, adplt=0;

  QDP_ColorMatrix *U[4], *Uf[4][4];
  for(int mu=0; mu<4; mu++) {
    //QDP_create_M(U[mu]);
    //QDP_M_eq_M(U[mu], gauge->links[mu], QDP_all);
    U[mu] = gauge->links[mu];
    for(int nu=0; nu<4; nu++) {
      if(nu==mu) continue;
      Uf[mu][nu] = QDP_create_M();
      QDP_M_eq_sM(Uf[mu][nu], U[mu], QDP_neighbor[nu], QDP_forward, QDP_all);
    }
  }

  QDP_ColorMatrix *UUf[4][4], *fstpl[4][4], *bstpl0[4][4], *bstpl[4][4];
  for(int mu=1; mu<4; mu++) {
    for(int nu=0; nu<mu; nu++) {
      if(pgm) {
	UUf[mu][nu] = QDP_create_M();
	UUf[nu][mu] = QDP_create_M();
	fstpl[mu][nu] = QDP_create_M();
	fstpl[nu][mu] = QDP_create_M();
	bstpl0[mu][nu] = QDP_create_M();
	bstpl0[nu][mu] = QDP_create_M();
	bstpl[mu][nu] = QDP_create_M();
	bstpl[nu][mu] = QDP_create_M();
      } else {
	if(mu==1) {
	  UUf[mu][nu] = QDP_create_M();
	  UUf[nu][mu] = QDP_create_M();
	  if(rect) {
	    fstpl[mu][nu] = QDP_create_M();
	    fstpl[nu][mu] = QDP_create_M();
	    bstpl0[mu][nu] = QDP_create_M();
	    bstpl0[nu][mu] = QDP_create_M();
	    bstpl[mu][nu] = QDP_create_M();
	    bstpl[nu][mu] = QDP_create_M();
	  }
	} else {
	  UUf[mu][nu] = UUf[1][0];
	  UUf[nu][mu] = UUf[0][1];
	  if(rect) {
	    fstpl[mu][nu] = fstpl[1][0];
	    fstpl[nu][mu] = fstpl[0][1];
	    bstpl0[mu][nu] = bstpl0[1][0];
	    bstpl0[nu][mu] = bstpl0[0][1];
	    bstpl[mu][nu] = bstpl[1][0];
	    bstpl[nu][mu] = bstpl[0][1];
	  }
	}
      }
      QDP_M_eq_M_times_M(UUf[mu][nu], U[mu], Uf[nu][mu], QDP_all);
      QDP_M_eq_M_times_M(UUf[nu][mu], U[nu], Uf[mu][nu], QDP_all);
      nflops += 2*EQMTM;
      if(adpl) {
	QDP_Complex *tc = QDP_create_C();
	QDP_C_eq_M_dot_M(tc, UUf[mu][nu], UUf[nu][mu], QDP_all);
	QLA_Complex z;
	if(plaq) QDP_c_eq_sum_C(&z, tc, QDP_all);
	else QLA_c_eq_r(z, 0);
	QLA_Real r;
	QDP_r_eq_norm2_C(&r, tc, QDP_all);
	nflops += 8*QLA_Nc*QLA_Nc-2 +2 +4;
	if(mu==3) {
	  plaqt += QLA_real(z);
	  adplt += r;
	} else {
	  plaqs += QLA_real(z);
	  adpls += r;
	}
	QDP_destroy_C(tc);
	//QOP_printf0("adpl[%i][%i] = %g\n", mu, nu, r);
	//QOP_printf0("adpls: %g  adplt: %g\n", adpls, adplt);
      } else
      if(plaq) {
	QLA_Real t;
	QDP_r_eq_re_M_dot_M(&t, UUf[mu][nu], UUf[nu][mu], QDP_all);
	nflops += 4*QLA_Nc*QLA_Nc;
	if(mu==3) {
	  plaqt += t;
	} else {
	  plaqs += t;
	}
      }
      if(rect||pgm) {
	QDP_M_eq_Ma_times_M(bstpl0[mu][nu], U[nu], UUf[mu][nu], QDP_all);
	QDP_M_eq_sM(bstpl[mu][nu], bstpl0[mu][nu], QDP_neighbor[nu],QDP_backward,QDP_all);
	QDP_M_eq_Ma_times_M(bstpl0[nu][mu], U[mu], UUf[nu][mu], QDP_all);
	QDP_M_eq_sM(bstpl[nu][mu], bstpl0[nu][mu], QDP_neighbor[mu],QDP_backward,QDP_all);
	QDP_M_eq_M_times_Ma(fstpl[mu][nu], UUf[nu][mu], Uf[nu][mu], QDP_all);
	QDP_M_eq_M_times_Ma(fstpl[nu][mu], UUf[mu][nu], Uf[mu][nu], QDP_all);
	nflops += 4*EQMTM;
	if(rect) {
	  QLA_Real t, tr=0;
	  QDP_r_eq_re_M_dot_M(&t, bstpl[mu][nu], fstpl[mu][nu], QDP_all);
	  tr += t;
	  QDP_r_eq_re_M_dot_M(&t, bstpl[nu][mu], fstpl[nu][mu], QDP_all);
	  tr += t;
	  nflops += 2*4*QLA_Nc*QLA_Nc;
	  if(mu==3) {
	    rectt += tr;
	  } else {
	    rects += tr;
	  }
	}
      }
    }
  }

#define combineb(x,a,b,c) {			\
  QLA_Real t; \
  QDP_M_eq_Ma_times_M(UUf[0][1], bstpl[a][c], bstpl[b][c], QDP_all); \
  QDP_M_eq_M_times_M(UUf[1][0], UUf[0][1], Uf[a][b], QDP_all); \
  QDP_r_eq_re_M_dot_M(&t, Uf[b][a], UUf[1][0], QDP_all); \
  x += t; \
  nflops += 2*EQMTM+4*QLA_Nc*QLA_Nc; \
}
#define combineb2(x,a,b,c,d) {			\
  QLA_Real t; \
  QDP_M_eq_Ma_times_M(UUf[0][1], bstpl[a][c], bstpl[b][c], QDP_all); \
  QDP_M_peq_Ma_times_M(UUf[0][1], bstpl[a][d], bstpl[b][d], QDP_all); \
  QDP_M_eq_M_times_M(UUf[1][0], UUf[0][1], Uf[a][b], QDP_all); \
  QDP_r_eq_re_M_dot_M(&t, Uf[b][a], UUf[1][0], QDP_all); \
  x += t; \
  nflops += 2*EQMTM+PEQMTM+4*QLA_Nc*QLA_Nc; \
}
#define combinefb(x,a,b,c) {			\
  QLA_Real t; \
  QDP_M_eq_Ma_times_M(UUf[0][1], fstpl[a][c], fstpl[b][c], QDP_all); \
  QDP_M_peq_Ma_times_M(UUf[0][1], bstpl[a][c], bstpl[b][c], QDP_all); \
  QDP_M_eq_M_times_M(UUf[1][0], UUf[0][1], Uf[a][b], QDP_all); \
  QDP_r_eq_re_M_dot_M(&t, Uf[b][a], UUf[1][0], QDP_all); \
  x += t; \
  nflops += 2*EQMTM+PEQMTM+4*QLA_Nc*QLA_Nc; \
}
#define combinefb2(x,a,b,c,d) {			\
  QLA_Real t; \
  QDP_M_eq_Ma_times_M(UUf[0][1], fstpl[a][c], fstpl[b][c], QDP_all); \
  QDP_M_peq_Ma_times_M(UUf[0][1], bstpl[a][c], bstpl[b][c], QDP_all); \
  QDP_M_peq_Ma_times_M(UUf[0][1], fstpl[a][d], fstpl[b][d], QDP_all); \
  QDP_M_peq_Ma_times_M(UUf[0][1], bstpl[a][d], bstpl[b][d], QDP_all); \
  QDP_M_eq_M_times_M(UUf[1][0], UUf[0][1], Uf[a][b], QDP_all); \
  QDP_r_eq_re_M_dot_M(&t, Uf[b][a], UUf[1][0], QDP_all); \
  x += t; \
  nflops += 2*EQMTM+3*PEQMTM+4*QLA_Nc*QLA_Nc; \
}
  if(pgm) {
    // 0,1,2
    combinefb(pgms,0,1,2);
    combineb(pgms,0,2,1);
    combineb(pgms,1,2,0);

    // rest
    //combinefb(pgmt,0,3,1);
    //combinefb(pgmt,0,3,2);
    combinefb2(pgmt,0,3,1,2);
    combineb(pgmt,0,1,3);
    combineb(pgmt,0,2,3);
    combinefb(pgmt,1,2,3);
    //combineb(pgmt,1,3,0);
    //combineb(pgmt,1,3,2);
    combineb2(pgmt,1,3,0,2);
    //combineb(pgmt,2,3,0);
    //combineb(pgmt,2,3,1);
    combineb2(pgmt,2,3,0,1);
  }

  *acts = plaq*plaqs + rect*rects + pgm*pgms + adpl*adpls;
  *actt = plaq*plaqt + rect*rectt + pgm*pgmt + adpl*adplt;

  for(int mu=0; mu<4; mu++) {
    for(int nu=0; nu<4; nu++) {
      if(nu==mu) continue;
      QDP_destroy_M(Uf[mu][nu]);
    }
  }
  if(pgm) {
    for(int mu=0; mu<4; mu++) {
      for(int nu=0; nu<4; nu++) {
	if(nu==mu) continue;
	QDP_destroy_M(UUf[mu][nu]);
	QDP_destroy_M(fstpl[mu][nu]);
	QDP_destroy_M(bstpl0[mu][nu]);
	QDP_destroy_M(bstpl[mu][nu]);
      }
    }
  } else {
    QDP_destroy_M(UUf[1][0]);
    QDP_destroy_M(UUf[0][1]);
    if(rect) {
      QDP_destroy_M(fstpl[1][0]);
      QDP_destroy_M(fstpl[0][1]);
      QDP_destroy_M(bstpl0[1][0]);
      QDP_destroy_M(bstpl0[0][1]);
      QDP_destroy_M(bstpl[1][0]);
      QDP_destroy_M(bstpl[0][1]);
    }
  }

  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
#undef NC
}
