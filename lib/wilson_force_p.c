//#define DO_TRACE
#include <qop_internal.h>

static void
project2(QDP_DiracFermion *xx2, QDP_DiracFermion *yy2, QDP_DiracFermion *xx, QDP_DiracFermion *yy, int mu)
{
#define NC QDP_get_nc(xx)
  int i;
  QDP_loop_sites(i, QDP_all, {
      QLA_DiracFermion(*x2) = QDP_site_ptr_readwrite_D(xx2, i);
      QLA_DiracFermion(*y2) = QDP_site_ptr_readwrite_D(yy2, i);
      QLA_DiracFermion(*x) = QDP_site_ptr_readonly_D(xx, i);
      QLA_DiracFermion(*y) = QDP_site_ptr_readonly_D(yy, i);
      QLA_D_eq_spproj_D(y2, y, mu, -1);
      QLA_D_peq_spproj_D(y2, x, mu, 1);
      QLA_D_eq_spproj_D(x2, x, mu, -1);
      QLA_D_peq_spproj_D(x2, y, mu, 1);
    });
#undef NC
}

static void
add_force(QDP_ColorMatrix *ff, QDP_DiracFermion *xx, QDP_DiracFermion *yy, QLA_Real eps)
{
#define NC QDP_get_nc(xx)
  int i;
  QDP_loop_sites(i, QDP_all, {
      QLA_ColorMatrix(*f) = QDP_site_ptr_readwrite_M(ff, i);
      QLA_DiracFermion(*x) = QDP_site_ptr_readonly_D(xx, i);
      QLA_DiracFermion(*y) = QDP_site_ptr_readonly_D(yy, i);
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
#undef NC
}

static void
wilson_deriv_multi_qdp(QOP_info_t *info,
		       QOP_FermionLinksWilson *flw,
		       QDP_ColorMatrix *deriv[],
		       QLA_Real eps[],
		       QDP_DiracFermion *x[],
		       QDP_DiracFermion *y[],
		       int n)
{
#define NC QDP_get_nc(flw->links[0])
  QDP_DiracFermion *x2[4], *y2[4], *ys[4];
  QDP_Lattice *lat = QDP_get_lattice_D(x[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  for(int mu=0; mu<4; mu++) {
    x2[mu] = QDP_create_D_L(lat);
    y2[mu] = QDP_create_D_L(lat);
    ys[mu] = QDP_create_D_L(lat);
  }
  for(int i=0; i<n; i++) {
    for(int mu=0; mu<4; mu++) {
      project2(x2[mu], y2[mu], x[i], y[i], mu);
      QDP_D_eq_sD(ys[mu], y2[mu], QDP_neighbor[mu], QDP_forward, QDP_all);
    }
    for(int mu=0; mu<4; mu++) {
      // 0.5 compensates for 2 in project2
      // -0.5 for dslash normalization
      add_force(deriv[mu], x2[mu], ys[mu], -0.25*eps[i]);
      QDP_discard_D(ys[mu]);
    }
  }
  info->final_flop = (4.*n*(32*QLA_Nc+34*QLA_Nc*QLA_Nc))*sites_on_node; 
  for(int mu=0; mu<4; mu++) {
    QDP_destroy_D(x2[mu]);
    QDP_destroy_D(y2[mu]);
    QDP_destroy_D(ys[mu]);
  }
#undef NC
}


// x^+ [(d/dX) D] y + y^+ [(d/dX) D^+] x
void
QOP_wilson_deriv_multi_qdp(QOP_info_t *info,
			   QOP_FermionLinksWilson *flw,
			   QDP_ColorMatrix *deriv[],
			   QLA_Real eps[],
			   QDP_DiracFermion *x[],
			   QDP_DiracFermion *y[],
			   int n)
{
#if 0
  QDP_ColorMatrix *d[4];
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  if(0) {
    //if(flw->gauge && flw->gauge->chained &&
    //(flw->gauge->nparents || doLastScale)) { // apply chain rule
    for(int mu=0; mu<4; mu++) {
      d[mu] = QDP_create_M_L(lat);
      QDP_M_eq_zero(d[mu], QDP_all);
    }
  } else {
    for(int mu=0; mu<4; mu++) {
      d[mu] = deriv[mu];
    }
  }
#endif

  wilson_deriv_multi_qdp(info, flw, deriv, eps, x, y, n);

#if 0
  if(0) {
    //if(flw->gauge && flw->gauge->chained &&
    //(flw->gauge->nparents || doLastScale)) { // apply chain rule
    //warning: passing argument 4 of 'QOP_F3_gauge_deriv_multi_qdp' from incompatible pointer type [enabled by default]
    //note: expected 'struct QDP_F3_ColorMatrix ***' but argument is of type 'struct QDP_F3_ColorMatrix * (*)[4]'
    // QOP_gauge_deriv_multi_qdp(info, deriv, &flw->gauge, &d, 1);
    for(int mu=0; mu<4; mu++) {
      QDP_destroy_M(d[mu]);
    }
  }
#endif
}

// antiherm x^+ [U(d/dU) D] y + y^+ [U(d/dU) D^+] x
void
QOP_wilson_force_multi_qdp(QOP_info_t *info,
			   QOP_FermionLinksWilson *flw,
			   QDP_ColorMatrix *force[],
			   QLA_Real eps[],
			   QDP_DiracFermion *x[],
			   QDP_DiracFermion *y[],
			   int n)
{
#define NC QDP_get_nc(flw->links[0])
  QDP_ColorMatrix *deriv[4];
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  for(int mu=0; mu<4; mu++) {
    deriv[mu] = QDP_create_M_L(lat);
    QDP_M_eq_zero(deriv[mu], QDP_all);
  }
  // factor of -2 for GL -> U
  QLA_Real s = -2;
  QDP_ColorMatrix **links;
  //if(flw->gauge && flw->gauge->chained && flw->gauge->nparents) {
  if(0) {
    QOP_GaugeField *top = flw->gauge;
    while(top->chained && top->nparents) top = top->parents[0];
    links = top->links;
  } else {
    links = flw->links;
    s *= -2; // compensate for -0.5 dslash normalization already in links
  }
  QLA_Real teps[n];
  for(int i=0; i<n; i++) teps[i] = s*eps[i];
  QOP_wilson_deriv_multi_qdp(info, flw, deriv, teps, x, y, n);
  QDP_ColorMatrix *t = QDP_create_M_L(lat);
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M_times_Ma(t, links[mu], deriv[mu], QDP_all);
    QDP_M_eq_antiherm_M(deriv[mu], t, QDP_all);
    QDP_M_peq_M(force[mu], deriv[mu], QDP_all);
    QDP_destroy_M(deriv[mu]);
  }
  info->final_flop += (4.*(QLA_Nc*QLA_Nc*(8*QLA_Nc+2)))*sites_on_node; 
  QDP_destroy_M(t);
#undef NC
}


// A = 1 - kappa^2 H_eo H_oe
// A = 1 - 4 kappa^2 D_eo D_oe

// (d/dU) [ x^+ A y + y^+ A^+ x ]
// -4 kappa^2 [ x^+, x^+ D_eo ] [(d/dU) D] [ y; D_oe y ] + h.c.
void
QOP_wilson_deriv_prec_multi_qdp(QOP_info_t *info,
				QOP_FermionLinksWilson *flw,
				QDP_ColorMatrix *deriv[],
				QLA_Real kappa[],
				QLA_Real eps[],
				QDP_DiracFermion *x[],
				QDP_DiracFermion *y[],
				int n)
{
#define NC QDP_get_nc(flw->links[0])
  double dtime = QOP_time();
  QDP_Lattice *lat = QDP_get_lattice_D(x[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);

  QLA_Real teps[n];
  for(int i=0; i<n; i++) {
    QOP_wilson_dslash_qdp(info, flw, kappa[i], -1, x[i], x[i], QOP_ODD, QOP_EVEN);
    QOP_wilson_dslash_qdp(info, flw, kappa[i],  1, y[i], y[i], QOP_ODD, QOP_EVEN);
    teps[i] = -4*kappa[i]*kappa[i]*eps[i];
  }
  QOP_wilson_deriv_multi_qdp(info, flw, deriv, teps, x, y, n);

  double nflop = 8*(16*QLA_Nc+9)*QLA_Nc;
  //if(flw->clov!=NULL) nflop += 32*(2*QLA_Nc-1)*QLA_Nc;
  info->final_flop += nflop*n*sites_on_node; 
  info->final_sec = QOP_time() - dtime;
#undef NC
}

// antiherm U(d/dU) [ x^+ A y + y^+ A^+ x ]
// -4 kappa^2 [ x^+, x^+ D_eo ] [U(d/dU) D] [ y; D_oe y ] + h.c.
void
QOP_wilson_force_prec_multi_qdp(QOP_info_t *info,
				QOP_FermionLinksWilson *flw,
				QDP_ColorMatrix *force[],
				QLA_Real kappa[],
				QLA_Real eps[],
				QDP_DiracFermion *x[],
				QDP_DiracFermion *y[],
				int n)
{
#define NC QDP_get_nc(flw->links[0])
  double dtime = QOP_time();
  QDP_Lattice *lat = QDP_get_lattice_D(x[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);

  QLA_Real teps[n];
  for(int i=0; i<n; i++) {
    QOP_wilson_dslash_qdp(info, flw, kappa[i], -1, x[i], x[i], QOP_ODD, QOP_EVEN);
    QOP_wilson_dslash_qdp(info, flw, kappa[i],  1, y[i], y[i], QOP_ODD, QOP_EVEN);
    teps[i] = -4*kappa[i]*kappa[i]*eps[i];
  }
  QOP_wilson_force_multi_qdp(info, flw, force, teps, x, y, n);

  double nflop = 8*(16*QLA_Nc+9)*QLA_Nc;
  //if(flw->clov!=NULL) nflop += 32*(2*QLA_Nc-1)*QLA_Nc;
  info->final_flop += nflop*n*sites_on_node; 
  info->final_sec = QOP_time() - dtime;
#undef NC
}

#if 0

static void
fmunu_deriv(QOP_info_t *info,
	    QDP_ColorMatrix *links[],
	    QDP_ColorMatrix *deriv[],
	    QDP_ColorMatrix *mid,
	    QLA_Real scale,
	    int mu, int nu)
{
  // deriv[mu] += (scale/8) [ mid Unu Umufnu Unufmu+ + Unu mid+ Umufnu Unufmu+
  //  + Unu Umufnu mid Unufmu+ + Unu Umufnu Unufmu+ mid+
  //  + mid Unubnu+ Umubnu Unufmubnu + Unubnu+ mid+ Umubnu Unufmubnu
  //  + Unubnu+ Umubnu mid Unufmubnu + Unubnu+ Umubnu Unufmubnu mid+ ]
}

void
QOP_wilson_clover_deriv_multi_qdp(QOP_info_t *info,
				  QOP_FermionLinksWilson *flw,
				  QDP_ColorMatrix *deriv[],
				  QLA_Real kappa[],
				  QLA_Real eps[],
				  QDP_DiracFermion *x[],
				  QDP_DiracFermion *y[],
				  int n)
{
#define NC QDP_get_nc(flw->links[0])
  double dtime = QOP_time();
  QDP_Lattice *lat = QDP_get_lattice_D(x[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);

  QLA_Real teps[n];
  for(int i=0; i<n; i++) {
    QOP_wilson_dslash_qdp(info, flw, kappa[i], -1, x[i], x[i], QOP_ODD, QOP_EVEN);
    QOP_wilson_dslash_qdp(info, flw, kappa[i],  1, y[i], y[i], QOP_ODD, QOP_EVEN);
    teps[i] = -4*kappa[i]*kappa[i]*eps[i];
  }
  QOP_wilson_deriv_multi_qdp(info, flw, deriv, teps, x, y, n);

  double nflop = 8*(16*QLA_Nc+9)*QLA_Nc;
  //if(flw->clov!=NULL) nflop += 32*(2*QLA_Nc-1)*QLA_Nc;
  info->final_flop += nflop*n*sites_on_node; 
  info->final_sec = QOP_time() - dtime;
#undef NC
}

#endif
