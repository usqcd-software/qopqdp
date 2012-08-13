//#define DO_TRACE
#include <qop_internal.h>

static void
project2(QDP_DiracFermion *xx2, QDP_DiracFermion *yy2, QDP_DiracFermion *xx, QDP_DiracFermion *yy, int mu)
{
  int i;
  QDP_loop_sites(i, QDP_all, {
      QLA_DiracFermion *x2 = QDP_site_ptr_readwrite_D(xx2, i);
      QLA_DiracFermion *y2 = QDP_site_ptr_readwrite_D(yy2, i);
      QLA_DiracFermion *x = QDP_site_ptr_readonly_D(xx, i);
      QLA_DiracFermion *y = QDP_site_ptr_readonly_D(yy, i);
      QLA_D_eq_spproj_D(y2, y, mu, -1);
      QLA_D_peq_spproj_D(y2, x, mu, 1);
      QLA_D_eq_spproj_D(x2, x, mu, -1);
      QLA_D_peq_spproj_D(x2, y, mu, 1);
    });
}

static void
add_force(QDP_ColorMatrix *ff, QDP_DiracFermion *xx, QDP_DiracFermion *yy, QLA_Real eps)
{
  int i;
  QDP_loop_sites(i, QDP_all, {
      QLA_ColorMatrix *f = QDP_site_ptr_readwrite_M(ff, i);
      QLA_DiracFermion *x = QDP_site_ptr_readonly_D(xx, i);
      QLA_DiracFermion *y = QDP_site_ptr_readonly_D(yy, i);
      for(int ic=0; ic<QLA_Nc; ic++) {
	for(int jc=0; jc<QLA_Nc; jc++) {
	  QLA_Complex z;
	  QLA_c_eq_r(z, 0);
	  for(int is=0; is<QLA_Ns; is++) {
	    QLA_c_peq_c_times_ca(z, QLA_elem_D(*x,ic,is), QLA_elem_D(*y,jc,is));
	  }
	  QLA_c_peq_r_times_c(QLA_elem_M(*f,ic,jc), eps, z);
	}
      }
    });
}

// A = 1 - kappa^2 H_eo H_oe
// A = 1 - 4 kappa^2 D_eo D_oe

#if 0
void
QOPPC(wilson_force_prec_qdp)(QOP_info_t *info,
			     QOP_FermionLinksWilson *flw,
			     QOP_Force *force, 
			     REAL kappa,
			     REAL eps, 
			     QDP_DiracVector *xx,
			     QDP_DiracVector *yy,
			     int sign)
{
  QDP_DiracFermion *x, *y;
  QDP_HalfVector *yh, *ys;
  x = QDP_create_D();
  y = QDP_create_D();
  yh = QDP_create_H();
  ys = QDP_create_H();
  if(sign>0) {
    // x^+ dA/dU y
    QDP_D_eq_D(x, xx, QDP_even);
    QOP_wilson_dslash(info, flw, kappa, -1, x, xx, QOP_ODD, QOP_EVEN);
    QDP_D_eq_D(y, yy, QDP_even);
    QOP_wilson_dslash(info, flw, kappa, 1, y, yy, QOP_ODD, QOP_EVEN);
    for(mu=0; mu<4; mu++) {
      QDP_H_eq_spproj_D(yh, y, mu, -1, QDP_all);
      QDP_H_eq_sH(ys, yh, QDP_neighbor[mu], QDP_forward, QDP_all);
      add_force(force->force[mu], ys, x, eps*kappa);
      QDP_discard_H(ys);
    }
  } else {
    // x^+ dA^+/dU y
    QDP_D_eq_D(x, xx, QDP_even);
    QOP_wilson_dslash(info, flw, kappa, 1, x, xx, QOP_ODD, QOP_EVEN);
    QDP_D_eq_D(y, yy, QDP_even);
    QOP_wilson_dslash(info, flw, kappa, -1, y, yy, QOP_ODD, QOP_EVEN);
    for(mu=0; mu<4; mu++) {
      QDP_H_eq_spproj_D(yh, y, mu, 1, QDP_all);
      QDP_H_eq_sH(ys, yh, QDP_neighbor[mu], QDP_forward, QDP_all);
      add_force(force->force[mu], ys, x, eps*kappa);
      QDP_discard_H(ys);
    }
  }
  QDP_destroy_D(x);
  QDP_destroy_D(y);
  QDP_destroy_H(yh);
  QDP_destroy_H(ys);
}

void
QOPPC(wilson_force_multi_qdp)(QOP_info_t *info,
			      QOP_FermionLinksWilson *flw,
			      QOP_Force *force, 
			      REAL kappa[],
			      REAL eps[], 
			      QDP_DiracFermion *x[],
			      int n)
{
  QDP_DiracFermion *xs, *xp, *y, *ys, *yp;
  xs = QDP_create_D();
  xp = QDP_create_D();
  y = QDP_create_D();
  ys = QDP_create_D();
  yp = QDP_create_D();
  for(int i=0; i<n; i++) {
    for(int mu=0; mu<4; mu++) {
      QDP_H_eq_spproj_D(xs, x[i], mu, -1, QDP_all);
      QDP_H_eq_sH(xp, xs, QDP_neighbor[mu], QDP_forward, QDP_all);
      QOP_wilson_dslash(&tinfo, flw, kappa[i], 1, y, x[i], QOP_EVENODD, QOP_EVENODD);
      QDP_H_eq_spproj_D(ys, y, mu, 1, QDP_all);
      QDP_H_eq_sH(yp, ys, QDP_neighbor[mu], QDP_forward, QDP_all);
      add_force(force->force[mu], xp, y, eps[i]);
      add_force(force->force[mu], yp, x, eps[i]);
      QDP_discard_H(xp);
      QDP_discard_H(yp);
    }
  }
  QDP_destroy_D(xs);
  QDP_destroy_D(xp);
  QDP_destroy_D(y);
  QDP_destroy_D(ys);
  QDP_destroy_D(yp);
}
#endif


// x^+ [(d/dX) D] y + y^+ [(d/dX) D^+] x
void
QOPPC(wilson_deriv_multi_qdp)(QOP_info_t *info,
			   QOP_FermionLinksWilson *flw,
			   QDP_ColorMatrix *deriv[],
			   QLA_Real eps[],
			   QDP_DiracFermion *x[],
			   QDP_DiracFermion *y[],
			   int n)
{
  QDP_DiracFermion *x2[4], *y2[4], *ys[4];
  for(int mu=0; mu<4; mu++) {
    x2[mu] = QDP_create_D();
    y2[mu] = QDP_create_D();
    ys[mu] = QDP_create_D();
  }
  for(int i=0; i<n; i++) {
    for(int mu=0; mu<4; mu++) {
      project2(x2[mu], y2[mu], x[i], y[i], mu);
      QDP_D_eq_sD(ys[mu], y2[mu], QDP_neighbor[mu], QDP_forward, QDP_all);
    }
    for(int mu=0; mu<4; mu++) {
      // 0.5 compensates for 2 in project2
      add_force(deriv[mu], ys[mu], x2[mu], 0.5*eps[i]);
      QDP_discard_D(ys[mu]);
    }
  }
  info->final_flop = (4.*n*(60+126))*QDP_sites_on_node; 
  for(int mu=0; mu<4; mu++) {
    QDP_destroy_D(x2[mu]);
    QDP_destroy_D(y2[mu]);
    QDP_destroy_D(ys[mu]);
  }
}

// x^+ [U(d/dU) D] y + y^+ [U(d/dU) D^+] x
void
QOPPC(wilson_force_multi_qdp)(QOP_info_t *info,
			   QOP_FermionLinksWilson *flw,
			   QOP_Force *force,
			   QLA_Real eps[],
			   QDP_DiracFermion *x[],
			   QDP_DiracFermion *y[],
			   int n)
{
  QDP_ColorMatrix *deriv[4], *t, *t2;
  for(int mu=0; mu<4; mu++) {
    deriv[mu] = QDP_create_M();
    QDP_M_eq_zero(deriv[mu], QDP_all);
  }
  // factor of -2 for GL -> U
  QLA_Real teps[n];
  for(int i=0; i<n; i++) teps[i] = -2*eps[i];
  QOPPC(wilson_deriv_multi_qdp)(info, flw, deriv, teps, x, y, n);
  t = QDP_create_M();
  t2 = QDP_create_M();
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M_times_M(t, flw->links[mu], deriv[mu], QDP_all);
    QDP_M_eq_antiherm_M(t2, t, QDP_all);
    QDP_M_peq_M(force->force[mu], t2, QDP_all);
    QDP_destroy_M(deriv[mu]);
  }
  info->final_flop += (4.*(198+24+18))*QDP_sites_on_node; 
  QDP_destroy_M(t);
  QDP_destroy_M(t2);
}

// U(d/dU) [ x^+ A y + y^+ A^+ x ]
// -4 kappa^2 [ x^+, x^+ D_eo ] [U(d/dU) D] [ y; D_oe y ] + h.c.
void
QOP_wilson_force_prec_multi_qdp(QOP_info_t *info,
				QOP_FermionLinksWilson *flw,
				QOP_Force *force,
				QLA_Real kappa[],
				QLA_Real eps[],
				QDP_DiracFermion *x[],
				QDP_DiracFermion *y[],
				int n)
{
  double dtime = QOP_time();

  QLA_Real teps[n];
  for(int i=0; i<n; i++) {
    QOP_wilson_dslash_qdp(info, flw, kappa[i], -1, x[i], x[i], QOP_ODD, QOP_EVEN);
    QOP_wilson_dslash_qdp(info, flw, kappa[i],  1, y[i], y[i], QOP_ODD, QOP_EVEN);
    teps[i] = -4*kappa[i]*kappa[i]*eps[i];
  }
  QOPPC(wilson_force_multi_qdp)(info, flw, force, teps, x, y, n);

  info->final_flop += ((144+168*7)+48)*n*QDP_sites_on_node; 
  info->final_sec = QOP_time() - dtime;
}
