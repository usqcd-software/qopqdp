#include <qop_internal.h>

static int
is_ancestor(QOP_GaugeField *a, QOP_GaugeField *g)
{
  if(g->chained) {
    for(int i=0; i<g->nparents; i++) {
      QOP_GaugeField *t = g->parents[i];
      if((a==t) || is_ancestor(a, t)) return 1;
    }
  }
  return 0;
}

// (d/dU) \sum_i < G_i C_i^+ + G_i^+ C_i >
// where all the G's should have a common top parent gauge U
// and the chain fields C are treated as independent of U
// the conventions for the derivative are (d/dU) < U C^+ + U^+ C > = C
void
QOP_gauge_deriv_multi_qdp(QOP_info_t *info, QDP_ColorMatrix *deriv[],
			  QOP_GaugeField *g[], QDP_ColorMatrix **chain[],
			  int n, int doLastScale)
{
#define NC QDP_get_nc(deriv[0])
  QDP_Lattice *lat = QDP_get_lattice_M(*chain[0]);
  QDP_Subset all = QDP_all_L(lat);
  int ndim = QDP_ndim_L(lat);
  int k=0;
  while(k<n && (g[k]->chained==0 || g[k]->nparents==0)) k++;
  if(k<n) { // has parents
    // find oldest ancestor of g[k] in g, if any
    for(int i=0; i<n; i++) {
      if(i==k) continue;
      if(is_ancestor(g[k], g[i])) k = i;
    }
    QOP_GaugeField *gk = g[k];
    // make copies of arrays
    int np = gk->nparents;
    int nmax = n - 1 + np;
    QOP_GaugeField *gg[nmax];
    QDP_ColorMatrix **cc[nmax];
    for(int i=0; i<n; i++) {
      gg[i] = g[i];
      cc[i] = chain[i];
    }
    // find duplicates of g[k] in g and combine
    for(int i=k+1; i<n; i++) {
      if(gg[i]==gk) {
	for(int mu=0; mu<ndim; mu++) {
	  QDP_M_peq_M(cc[k][mu], cc[i][mu], all);
	}
	for(int j=i+1; j<n; j++) {
	  gg[j-1] = gg[j];
	  cc[j-1] = cc[j];
	}
	i--;
	n--;
      }
    }
    // apply scaling, if any
    if(gk->scale) {
      gk->scale(cc[k], gk, 0);
    }
    // get derivatives wrt parents
    QDP_ColorMatrix **d[np];
    for(int i=0; i<np; i++) {
      if(gk->parents[i]->chained && 
	 (gk->parents[i]->nparents || doLastScale)) {
	QOP_malloc(d[i], QDP_ColorMatrix *, ndim);
	for(int mu=0; mu<ndim; mu++) {
	  d[i][mu] = QDP_create_M_L(lat);
	  QDP_M_eq_zero(d[i][mu], all);
	}
      } else {
	d[i] = deriv;
      }
    }
    gk->deriv(d, gk, cc[k]);
    // collect remaining terms for chain rule
    for(int i=k+1; i<n; i++) {
      gg[i-1] = gg[i];
      cc[i-1] = cc[i];
    }
    n--;
    for(int i=0; i<np; i++) {
      if(gk->parents[i]->chained && 
	 (gk->parents[i]->nparents || doLastScale)) {
	gg[n] = gk->parents[i];
	cc[n] = d[i];
	n++;
      }
    }
    if(n>0) QOP_gauge_deriv_multi_qdp(info, deriv, gg, cc, n, doLastScale);
    for(int i=0; i<np; i++) {
      if(gk->parents[i]->chained && 
	 (gk->parents[i]->nparents || doLastScale)) {
	for(int mu=0; mu<ndim; mu++) {
	  QDP_destroy_M(d[i][mu]);
	}
	QOP_free(d[i]);
      }
    }
  } else { // no parents
    if(doLastScale) {
      for(int i=0; i<n; i++) {
	if(g[i]->chained && g[i]->scale) {
	  g[i]->scale(chain[i], g[i], 0);
	}
      }
    }
    for(int mu=0; mu<ndim; mu++) {
      for(int i=0; i<n; i++) {
	QDP_M_peq_M(deriv[mu], chain[i][mu], all);
      }
    }
  }
#undef NC
}

// convert gauge derivative to force (in SU(N) lie algebra)
void
QOP_gauge_force_multi_qdp(QOP_info_t *info, QDP_ColorMatrix *f[],
			  QOP_GaugeField *g[], QDP_ColorMatrix **chain[],
			  int n)
{
#define NC QDP_get_nc(f[0])
  QDP_Lattice *lat = QDP_get_lattice_M(*chain[0]);
  QDP_Subset all = QDP_all_L(lat);
  int ndim = QDP_ndim_L(lat);
  QDP_ColorMatrix *d[ndim];
  int cr = 0;
  for(int i=0; i<n; i++) {
    if(g[i]->chained && g[i]->nparents) {
      cr = 1;
      break;
    }
  }
  if(cr) { // apply chain rule
    for(int mu=0; mu<ndim; mu++) {
      d[mu] = QDP_create_M_L(lat);
      QDP_M_eq_zero(d[mu], all);
    }
    QOP_gauge_deriv_multi_qdp(info, d, g, chain, n, 0);
  } else{
    for(int mu=0; mu<ndim; mu++) {
      d[mu] = chain[0][mu];
    }
    for(int mu=0; mu<ndim; mu++) {
      for(int i=1; i<n; i++) {
	QDP_M_peq_M(d[mu], chain[i][mu], all);
      }
    }
  }

  QOP_GaugeField *top = g[0];
  while(top->chained && top->nparents) top = top->parents[0];

  QDP_ColorMatrix *m = QDP_create_M_L(lat);
  QLA_Real s = -2;
  for(int i=0; i<ndim; i++) {
    QDP_M_eq_M_times_Ma(m, top->links[i], d[i], all);
    QDP_M_eq_antiherm_M(d[i], m, all);
    QDP_M_peq_r_times_M(f[i], &s, d[i], all);
    QDP_destroy_M(d[i]);
  }
  QDP_destroy_M(m);
  if(cr) { // apply chain rule
    for(int mu=0; mu<ndim; mu++) {
      QDP_destroy_M(d[mu]);
    }
  }
#undef NC
}
