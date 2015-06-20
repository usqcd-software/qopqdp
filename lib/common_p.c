#include <stdlib.h>
#include <string.h>
#include <qop_internal.h>

// helpers

static QOP_GaugeField *
QOP_new_G(void)
{
  QOP_GaugeField *qopgf;
  QOP_malloc(qopgf, QOP_GaugeField, 1);
  QOP_malloc(qopgf->links, QDP_ColorMatrix *, QOP_common.ndim);
  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = NULL;
  }
  qopgf->raw = NULL;
  qopgf->parents = NULL;
  qopgf->nparents = 0;
  qopgf->deriv = NULL;
  qopgf->chained = 0;
  qopgf->scale = NULL;
  qopgf->r0 = NULL;
  qopgf->bc.phase = NULL;
  qopgf->sign.signmask = NULL;
  return qopgf;
}

/*******************/
/* public routines */
/*******************/

/* raw routines */

#define NC int nc
QOP_ColorVector *
QOP_create_V_from_raw(QDP_Lattice *lat, QOP_Real *src, QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_ColorVector *qopcv;

  QOP_malloc(qopcv, QOP_ColorVector, 1);
  qopcv->cv = QDP_create_V_L(lat);
  QOP_qdp_eq_raw(V, qopcv->cv, src, evenodd);
  qopcv->raw = NULL;
  return qopcv;
}
#undef NC

#define NC int nc
QOP_DiracFermion *
QOP_create_D_from_raw(QDP_Lattice *lat, QOP_Real *src, QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_DiracFermion *qopdf;

  QOP_malloc(qopdf, QOP_DiracFermion, 1);
  qopdf->df = QDP_create_D_L(lat);
  QOP_qdp_eq_raw(D, qopdf->df, src, evenodd);
  qopdf->raw = NULL;
  return qopdf;
}
#undef NC

#define NC int nc
QOP_GaugeField *
QOP_create_G_from_raw(QDP_Lattice *lat, QOP_Real *links[], QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_GaugeField *qopgf = QOP_new_G();

  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDP_create_M_L(lat);
    QOP_qdp_eq_raw(M, qopgf->links[i], links[i], evenodd);
  }
  return qopgf;
}
#undef NC

#define NC int nc
QOP_Force *
QOP_create_F_from_raw(QDP_Lattice *lat, QOP_Real *force[], QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_Force *qopf;
  int i;
  QOP_malloc(qopf, QOP_Force, 1);
  QOP_malloc(qopf->force, QDP_ColorMatrix *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDP_create_M_L(lat);
    QOP_qdp_eq_raw(M, qopf->force[i], force[i], evenodd);
  }
  qopf->raw = NULL;
  return qopf;
}
#undef NC

void
QOP_extract_V_to_raw(QOP_Real *dest, QOP_ColorVector *src,
			QOP_evenodd_t evenodd)
{
  QOP_raw_eq_qdp(V, dest, src->cv, evenodd);
}

void
QOP_extract_D_to_raw(QOP_Real *dest, QOP_DiracFermion *src,
			QOP_evenodd_t evenodd)
{
  QOP_raw_eq_qdp(D, dest, src->df, evenodd);
}

void
QOP_extract_G_to_raw(QOP_Real *dest[], QOP_GaugeField *src,
		     QOP_evenodd_t evenodd)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QOP_raw_eq_qdp(M, dest[i], src->links[i], evenodd);
  }
}

void
QOP_extract_F_to_raw(QOP_Real *dest[], QOP_Force *src, QOP_evenodd_t evenodd)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QOP_raw_eq_qdp(M, dest[i], src->force[i], evenodd);
  }
}

void
QOP_destroy_V(QOP_ColorVector *field)
{
  QDP_destroy_V(field->cv);
  free(field);
}

void
QOP_destroy_D(QOP_DiracFermion *field)
{
  QDP_destroy_D(field->df);
  free(field);
}

void
QOP_destroy_G(QOP_GaugeField *field)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QDP_destroy_M(field->links[i]);
  }
  free(field->links);
  free(field);
}

void
QOP_destroy_F(QOP_Force *field)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QDP_destroy_M(field->force[i]);
  }
  free(field->force);
  free(field);
}

#define NC int nc
QOP_ColorVector *
QOP_convert_V_from_raw(QDP_Lattice *lat, QOP_Real *src, QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_ColorVector *ret =
    QOP_create_V_from_raw(lat, src, evenodd);
  ret->raw = src;
  return ret;
}
#undef NC

#define NC int nc
QOP_DiracFermion *
QOP_convert_D_from_raw(QDP_Lattice *lat, QOP_Real *src, QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_DiracFermion *ret =
    QOP_create_D_from_raw(lat, src, evenodd);
  ret->raw = src;
  return ret;
}
#undef NC

#define NC int nc
QOP_GaugeField *
QOP_convert_G_from_raw(QDP_Lattice *lat, QOP_Real *links[], QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_GaugeField *ret =
    QOP_create_G_from_raw(lat, links, evenodd);
  ret->raw = links;
  return ret;
}
#undef NC

#define NC int nc
QOP_Force *
QOP_convert_F_from_raw(QDP_Lattice *lat, QOP_Real *force[], QOP_evenodd_t evenodd)
#undef NC
#define NC nc
{
  QOP_Force *ret =
    QOP_create_F_from_raw(lat, force, evenodd);
  ret->raw = force;
  return ret;
}
#undef NC

QOP_Real *
QOP_convert_V_to_raw(QDP_Lattice *lat, QOP_ColorVector *src, QOP_evenodd_t evenodd)
{
  QOP_Real *ret = src->raw;
  if(!ret) QOP_malloc(ret, QOP_Real, QOP_sites_on_node_raw_V(evenodd));
  QOP_extract_V_to_raw(ret, src, evenodd);
  QOP_destroy_V(src);
  return ret;
}

QOP_Real *
QOP_convert_D_to_raw(QDP_Lattice *lat, QOP_DiracFermion *src, QOP_evenodd_t evenodd)
{
  QOP_Real *ret = src->raw;
  if(!ret) QOP_malloc(ret, QOP_Real, QOP_sites_on_node_raw_D(evenodd));
  QOP_extract_D_to_raw(ret, src, evenodd);
  QOP_destroy_D(src);
  return ret;
}

QOP_Real **
QOP_convert_G_to_raw(QDP_Lattice *lat, QOP_GaugeField *src, QOP_evenodd_t evenodd)
{
  QOP_Real **ret = src->raw;
  if(!ret) {
    int i;
    QOP_malloc(ret, QOP_Real *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(ret[i], QOP_Real, QOP_sites_on_node_raw_G(evenodd));
    }
  }
  QOP_extract_G_to_raw(ret, src, evenodd);
  QOP_destroy_G(src);
  return ret;
}

QOP_Real **
QOP_convert_F_to_raw(QDP_Lattice *lat, QOP_Force *src, QOP_evenodd_t evenodd)
{
  QOP_Real **ret = src->raw;
  if(!ret) {
    int i;
    QOP_malloc(ret, QOP_Real *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(ret[i], QOP_Real, QOP_sites_on_node_raw_F(evenodd));
    }
  }
  QOP_extract_F_to_raw(ret, src, evenodd);
  QOP_destroy_F(src);
  return ret;
}


/* qdp routines */

#define NC QDP_get_nc(src)
QOP_ColorVector *
QOP_create_V_from_qdp(QDP_ColorVector *src)
{
  QOP_ColorVector *qopcv;
  QDP_Lattice *lat = QDP_get_lattice_V(src);

  QOP_malloc(qopcv, QOP_ColorVector, 1);
  qopcv->cv = QDP_create_V_L(lat);
  QDP_V_eq_V(qopcv->cv, src, QDP_all);
  qopcv->raw = NULL;
  return qopcv;
}

QOP_DiracFermion *
QOP_create_D_from_qdp(QDP_DiracFermion *src)
{
  QOP_DiracFermion *qopdf;
  QDP_Lattice *lat = QDP_get_lattice_D(src);

  QOP_malloc(qopdf, QOP_DiracFermion, 1);
  qopdf->df = QDP_create_D_L(lat);
  QDP_D_eq_D(qopdf->df, src, QDP_all);
  qopdf->raw = NULL;
  return qopdf;
}
#undef NC

#define NC QDP_get_nc(links[i])
QOP_GaugeField *
QOP_create_G_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_GaugeField *qopgf = QOP_new_G();
  QDP_Lattice *lat = QDP_get_lattice_M(links[0]);
  QDP_Subset all = QDP_all_L(lat);
  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDP_create_M_L(lat);
    QDP_M_eq_M(qopgf->links[i], links[i], all);
  }
  return qopgf;
}
#undef NC

#define NC QDP_get_nc(links[i])
QOP_GaugeField *
QOPPO(create_G_from_qdp)(QDPO(ColorMatrix) *links[])
{
  QOP_GaugeField *qopgf = QOP_new_G();
  QDP_Lattice *lat = QDPO(get_lattice_M)(links[0]);

  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDP_create_M_L(lat);
    QDPPO(M_eq_M)(qopgf->links[i], links[i], QDP_all);
  }
  return qopgf;
}
#undef NC

#define NC QDP_get_nc(force[i])
QOP_Force *
QOP_create_F_from_qdp(QDP_ColorMatrix *force[])
{
  QOP_Force *qopf;
  QDP_Lattice *lat = QDP_get_lattice_M(force[0]);
  int i;
  QOP_malloc(qopf, QOP_Force, 1);
  QOP_malloc(qopf->force, QDP_ColorMatrix *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDP_create_M_L(lat);
    QDP_M_eq_M(qopf->force[i], force[i], QDP_all);
  }
  qopf->raw = NULL;
  return qopf;
}
#undef NC

void
QOP_extract_V_to_qdp(QDP_ColorVector *dest, QOP_ColorVector *src)
{
  QDP_V_eq_V(dest, src->cv, QDP_all);
}

void
QOP_extract_D_to_qdp(QDP_DiracFermion *dest, QOP_DiracFermion *src)
{
  QDP_D_eq_D(dest, src->df, QDP_all);
}

void
QOP_extract_G_to_qdp(QDP_ColorMatrix *d[], QOP_GaugeField *src)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QDP_M_eq_M(d[i], src->links[i], QDP_all);
  }
}

void
QOP_extract_F_to_qdp(QDP_ColorMatrix *d[], QOP_Force *src)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QDP_M_eq_M(d[i], src->force[i], QDP_all);
  }
}

QOP_ColorVector *
QOP_convert_V_from_qdp(QDP_ColorVector *src)
{
  QOP_ColorVector *qopcv;
  QOP_malloc(qopcv, QOP_ColorVector, 1);
  qopcv->cv = src;
  qopcv->raw = NULL;
  return qopcv;
}

QOP_DiracFermion *
QOP_convert_D_from_qdp(QDP_DiracFermion *src)
{
  QOP_DiracFermion *qopdf;
  QOP_malloc(qopdf, QOP_DiracFermion, 1);
  qopdf->df = src;
  qopdf->raw = NULL;
  return qopdf;
}

QOP_GaugeField *
QOP_convert_G_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_GaugeField *qopgf = QOP_new_G();
  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = links[i];
  }
  return qopgf;
}

QOP_Force *
QOP_convert_F_from_qdp(QDP_ColorMatrix *force[])
{
  QOP_Force *qopf;
  int i;
  QOP_malloc(qopf, QOP_Force, 1);
  QOP_malloc(qopf->force, QDP_ColorMatrix *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = force[i];
  }
  qopf->raw = NULL;
  return qopf;
}

QDP_ColorVector *
QOP_convert_V_to_qdp(QOP_ColorVector *src)
{
  QDP_ColorVector *ret;
  ret = src->cv;
  free(src);
  return ret;
}

QDP_DiracFermion *
QOP_convert_D_to_qdp(QOP_DiracFermion *src)
{
  QDP_DiracFermion *ret;
  ret = src->df;
  free(src);
  return ret;
}

QDP_ColorMatrix **
QOP_convert_G_to_qdp(QOP_GaugeField *src)
{
  QDP_ColorMatrix **ret;
  ret = src->links;
  free(src);
  return ret;
}

QDP_ColorMatrix **
QOP_convert_F_to_qdp(QOP_Force *src)
{
  QDP_ColorMatrix **ret;
  ret = src->force;
  free(src);
  return ret;
}

/* rephase */

static int nd;
static int bc_dir;
static int bc_coord;
static int bc_origin;
static QLA_Complex bc_phase;
static int staggered_sign_bits;

#define NC nc
static void
rephase_func(NCPROT1 QLA_ColorMatrix(*m), int coords[])
{
  if(bc_dir>=0) {
    int rshift = (coords[bc_dir] + bc_coord - bc_origin) % bc_coord;
    if(rshift == bc_coord-1) {
      QLA_ColorMatrix(t);
      QLA_M_eq_c_times_M(&t, &bc_phase, m);
      QLA_M_eq_M(m, &t);
    }
  }
  if(staggered_sign_bits) {
    int s=0;
    for(int i=0; i<nd; i++) s += ((staggered_sign_bits>>i)&1) * coords[i];
    if(s&1) {
      QLA_M_eqm_M(m, m);
    }
  }
}
#undef NC

static void
scale(QDP_ColorMatrix *l[], QOP_GaugeField *g, int inv)
{
  nd = QDP_ndim();
  for(int i=0; i<nd; i++) {
    bc_dir = -1;
    staggered_sign_bits = 0;
    if(g->bc.phase && (g->bc.phase[i].re!=1. || g->bc.phase[i].im!=0.)) {
      bc_dir = i;
      bc_coord = QDP_coord_size(i);
      bc_origin = g->r0[i];
      if(inv) {
	QLA_Complex z;
	QLA_c_eq_r_plus_ir(z, g->bc.phase[i].re, g->bc.phase[i].im);
	QLA_c_eq_r_div_c(bc_phase, 1, z);
      } else {
	QLA_c_eq_r_plus_ir(bc_phase, g->bc.phase[i].re, g->bc.phase[i].im);
      }
    }
    if(g->sign.signmask) {
      staggered_sign_bits = g->sign.signmask[i];
    }
    QDP_M_eq_funct(l[i], rephase_func, QDP_all);
  }
}

void
QOP_rephase_G(QOP_GaugeField *links,
	      int *r0,
	      QOP_bc_t *bc,
	      QOP_staggered_sign_t *sign)
{
  nd = QDP_ndim();
  //links->scale = scale;
  //links->chained = 1;
#define copyarray(d,s,t,n) do{ QOP_malloc(d,t,n); memcpy(d,s,(n)*sizeof(t)); }while(0)
  copyarray(links->r0, r0, int, nd);
  copyarray(links->bc.phase, bc->phase, QOP_Complex, nd);
  copyarray(links->sign.signmask, sign->signmask, int, nd);
  scale(links->links, links, 0);
}

void
QOP_rephase_G_qdp(QDP_ColorMatrix *links[],
		  int *r0,
		  QOP_bc_t *bc,
		  QOP_staggered_sign_t *sign)
{
  nd = QDP_ndim();
  QOP_GaugeField g;
  //links->scale = scale;
  //links->chained = 1;
#define copyarray(d,s,t,n) do{ QOP_malloc(d,t,n); memcpy(d,s,(n)*sizeof(t)); }while(0)
  copyarray(g.r0, r0, int, nd);
  copyarray(g.bc.phase, bc->phase, QOP_Complex, nd);
  copyarray(g.sign.signmask, sign->signmask, int, nd);
  scale(links, &g, 0);
}

// qll routines
#ifdef HAVE_QLL

#include <qll.h>

#if QOP_Precision == 'F'
#define P(x) x ## F
#if QOP_Colors == 'N'
#define PC(x) APPLY(CAT3,x,F,N)
#else
#define PC(x) APPLY(CAT3,x,F,QOP_Colors)
#endif
typedef float real;
#else
#define P(x) x ## D
#if QOP_Colors == 'N'
#define PC(x) APPLY(CAT3,x,D,N)
#else
#define PC(x) APPLY(CAT3,x,D,QOP_Colors)
#endif
typedef double real;
#endif

#if QOP_Colors == 3

static int setup = 0;
static Layout layout;

static void
get_rankGeom(QDP_Lattice *lat, int myrank, int nd, int *ls, int *rg)
{
  int n = QDP_numsites_L(lat, myrank);
  int xmin[nd], xmax[nd];
  for(int i=0; i<nd; i++) {
    xmin[i] = ls[i];
    xmax[i] = 0;
  }
  for(int i=0; i<n; i++) {
    int x[nd];
    QDP_get_coords(x, myrank, i);
    for(int j=0; j<nd; j++) {
      if(x[j]<xmin[j]) xmin[j] = x[j];
      if(x[j]>xmax[j]) xmax[j] = x[j];
    }
  }
  for(int i=0; i<nd; i++) {
    rg[i] = ls[i]/(1+xmax[i]-xmin[i]);
  }
}

void
setup_qll(QDP_Lattice *lat)
{
  if(!setup) {
    setup = 1;
    layout.nranks = QMP_get_number_of_nodes();
    layout.myrank = QMP_get_node_number();
    int nd = QDP_ndim_L(lat);
    int ls[nd], rg[nd];
    QDP_latsize_L(lat, ls);
    get_rankGeom(lat, layout.myrank, nd, ls, rg);
    P(stagDslashSetup)(&layout, nd, ls, rg);
    //printf("done %s\n", __func__);
  }
}

void *
get_qll_layout(void)
{
  return (void *)&layout;
}

void
toQDP(real *xx, QDP_Lattice *lat, real *yy, void *ll, int nelem)
{
  Layout *l = (Layout *)ll;
  int nd = l->nDim;
  int nl = l->nSitesInner;
  int ysize = nelem * nl;
  //double n2=0;
  int nsites = l->nSites;
  if(nsites!=QDP_sites_on_node_L(lat)) {
    printf("%s: nsites(%i) != QDP_sites_on_node_L(lat)(%i)\n",
           __func__, nsites, QDP_sites_on_node_L(lat));
    QDP_abort(-1);
  }
#pragma omp parallel for
  for(int i=0; i<nsites; i++) { // QDP sites
    int x[nd];
    LayoutIndex li;
    QDP_get_coords_L(lat, x, QDP_this_node, i);
    layoutIndex(l, &li, x);
    if(li.rank==l->myrank) {
      int oi = li.index / l->nSitesInner;
      int ii = li.index % l->nSitesInner;
      real *yi = yy + (ysize*oi + ii);
      for(int e=0; e<nelem; e++) {
        xx[i*nelem+e] = yi[nl*e];
        //n2 += er*er + ei*ei;
      }
    } else {
      printf("unpack: site on wrong node!\n");
      QDP_abort(-1);
    }
  }
  //printf("unpack2: %g\n", n2);
}

void
fromQDP(real *yy, void *ll, real *xx, QDP_Lattice *lat, int nelem)
{
  Layout *l = (Layout *)ll;
  int nd = l->nDim;
  int nl = l->nSitesInner;
  int ysize = nelem * nl;
  //double n2=0;
  int nsites = l->nSites;
  if(nsites!=QDP_sites_on_node_L(lat)) {
    printf("%s: nsites(%i) != QDP_sites_on_node_L(lat)(%i)\n",
           __func__, nsites, QDP_sites_on_node_L(lat));
    QDP_abort(-1);
  }
#pragma omp parallel for
  for(int j=0; j<nsites; j++) { // qll sites
    int x[nd];
    LayoutIndex li;
    li.rank = l->myrank;
    li.index = j;
    layoutCoord(l, x, &li);
    int r = QDP_node_number_L(lat, x);
    if(r==QDP_this_node) {
      int i = QDP_index_L(lat, x);
      int oi = j / l->nSitesInner;
      int ii = j % l->nSitesInner;
      real *yi = yy + (ysize*oi + ii);
      for(int e=0; e<nelem; e++) {
        yi[nl*e] = xx[i*nelem+e];
        //n2 += er*er + ei*ei;
      }
    } else {
      printf("unpack: site on wrong node!\n");
      QDP_abort(-1);
    }
  }
  //printf("unpack2: %g\n", n2);
}

#endif // QOP_Colors == 3

void
PC(unpackV)(Layout *l, QDP_ColorVector *xx, real *y)
{
#define NC QDP_get_nc(xx)
  QDP_Lattice *lat = QDP_get_lattice_V(xx);
  real *x = QDP_expose_V(xx);
  P(toQDP)(x, lat, y, l, 2*QDP_Nc);
  QDP_reset_V(xx);
#undef NC
}

void
PC(unpackM)(Layout *l, QDP_ColorMatrix *m, real *y)
{
#define NC QDP_get_nc(m)
  QDP_Lattice *lat = QDP_get_lattice_M(m);
  QDP_ColorMatrix *xx = QDP_create_M_L(lat);
  real *x = QDP_expose_M(xx);
  int nelem = 2*QDP_Nc*QDP_Nc;
  P(toQDP)(x, lat, y, l, nelem);
  QDP_reset_M(xx);
  //QDP_M_eq_M(m, xx, QDP_all_L(lat));
  QDP_M_eq_transpose_M(m, xx, QDP_all_L(lat));
  QDP_destroy_M(xx);
#undef NC
}

void
PC(packV)(Layout *l, real *x, QDP_ColorVector *yy)
{
#define NC QDP_get_nc(yy)
  QDP_Lattice *lat = QDP_get_lattice_V(yy);
  real *y = QDP_expose_V(yy);
  P(fromQDP)(x, l, y, lat, 2*QDP_Nc);
  QDP_reset_V(yy);
#undef NC
}

void
PC(packM)(Layout *l, real *y, QDP_ColorMatrix *m)
{
#define NC QDP_get_nc(m)
  QDP_Lattice *lat = QDP_get_lattice_M(m);
  QDP_ColorMatrix *xx = QDP_create_M_L(lat);
  //QDP_M_eq_M(xx, m, QDP_all_L(lat));
  QDP_M_eq_transpose_M(xx, m, QDP_all_L(lat));
  real *x = QDP_expose_M(xx);
  int nelem = 2*QDP_Nc*QDP_Nc;
  P(fromQDP)(y, l, x, lat, nelem);
  QDP_reset_M(xx);
  QDP_destroy_M(xx);
#undef NC
}

void *
create_qll_gauge(int nc)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  int nd = layout->nDim;
  int nelemM = 2*nc*nc;
  P(Field) *f;
  QOP_malloc(f, P(Field), nd);
  for(int i=0; i<nd; i++) {
    P(fieldNew)(&f[i], layout, nelemM);
  }
  return (void *)f;
}

void *
create_qll_from_gauge(QDP_ColorMatrix *g[])
{
#define NC QDP_get_nc(g[0])
  Layout *layout = (Layout *)P(get_qll_layout)();
  int nd = layout->nDim;
  int nelemM = 2*QDP_Nc*QDP_Nc;
  P(Field) *f;
  QOP_malloc(f, P(Field), nd);
  for(int i=0; i<nd; i++) {
    P(fieldNew)(&f[i], layout, nelemM);
    PC(packM)(layout, f[i].f, g[i]);
  }
  return (void *)f;
#undef NC
}

void
copy_gauge_from_qll(QDP_ColorMatrix *g[], void *ff)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  int nd = layout->nDim;
  P(Field) *f = (P(Field)*)ff;
  for(int i=0; i<nd; i++) {
    PC(unpackM)(layout, g[i], f[i].f);
  }
}

void
free_qll_gauge(void *ff)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  int nd = layout->nDim;
  P(Field) *f = (P(Field)*)ff;
  for(int i=0; i<nd; i++) {
    P(fieldFree)(&f[i]);
  }
}

void
fat7_qll(void *qllfl, void *qllll, QOP_asqtad_coeffs_t *coef,
	 void *qllu, void *qllul)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  int nd = layout->nDim;
  P(Field) *fl[nd], *ll[nd], *u[nd], *ul[nd], **llp=NULL, **ulp=NULL;
  for(int i=0; i<nd; i++) {
    fl[i] = i + ((P(Field)*)qllfl);
    u[i] = i + ((P(Field)*)qllu);
    ll[i] = i + ((P(Field)*)qllll);
    ul[i] = i + ((P(Field)*)qllul);
  }
  if(qllll) {
    llp = ll;
    ulp = ul;
  }
  P(smearFat7)(fl, llp, (Fat7Coeffs*)coef, u, ulp);
}

static P(StaggeredSolver) sse;
static QLA_Real **u=NULL, **u3=NULL;
static int nl = 0;

void
setup_qll_solver(QOP_FermionLinksAsqtad *fla)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  if(!fla->dblstored) {
    if(layout->myrank==0) printf("error: links not double stored\n");
  }
  nl = fla->nlinks;
  //printf("nl: %i\n", nl);
  if(nl==8) { // no naik
    u = myalloc(nl*sizeof(QLA_Real*));
    u3 = NULL;
    for(int i=0; i<nl; i++) {
      u[i] = myalloc(layout->nSites*18*sizeof(QLA_Real));
      PC(packM)(layout, u[i], fla->dbllinks[i]);
    }
  } else { // naik
    nl /= 2;
    u = myalloc(nl*sizeof(QLA_Real*));
    u3 = myalloc(nl*sizeof(QLA_Real*));
    for(int i=0; i<nl; i++) {
      u[i] = myalloc(layout->nSites*18*sizeof(QLA_Real));
      u3[i] = myalloc(layout->nSites*18*sizeof(QLA_Real));
      PC(packM)(layout, u[i], fla->dbllinks[2*i]);
      PC(packM)(layout, u3[i], fla->dbllinks[2*i+1]);
    }
  }
  P(createStaggeredSolver)(&sse, layout, u, u3, "even");
  sse.noscale = 1;
}

void
free_qll_solver(void)
{
  P(freeStaggeredSolver)(&sse);
  //printf("nl: %i\n", nl);
  //printf("u: %p\n", u);
  if(u) {
    for(int i=0; i<nl; i++) {
      free(u[i]);
    }
    free(u);
    u = NULL;
  }
  if(u3) {
    for(int i=0; i<nl; i++) {
      free(u3[i]);
    }
    free(u3);
    u3 = NULL;
  }
  //printf("done %s\n", __func__);
}

void
solve_qll(QDP_ColorVector *dest, QDP_ColorVector *src, double mass,
	  QOP_invert_arg_t *inv_arg, QOP_resid_arg_t *res_arg)
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  QLA_Real *s = myalloc(layout->nSites*6*sizeof(QLA_Real));
  QLA_Real *d = myalloc(layout->nSites*6*sizeof(QLA_Real));
  PC(packV)(layout, s, src);
  PC(packV)(layout, d, dest);
  double rsq = res_arg->rsqmin;
  int maxits = inv_arg->max_iter;
  int its = P(solve2)(&sse, d, mass, s, rsq, maxits, "even");
  res_arg->final_iter = its;
  PC(unpackV)(layout, dest, d);
  free(s);
  free(d);
}

void
solveMulti_qll(QDP_ColorVector *dest[], QDP_ColorVector *src, double ms[],
	       int nm,  QOP_invert_arg_t *invarg,
		  QOP_resid_arg_t *resargs[])
{
  Layout *layout = (Layout *)P(get_qll_layout)();
  QLA_Real *s = myalloc(layout->nSites*6*sizeof(QLA_Real));
  PC(packV)(layout, s, src);
  QLA_Real *d[nm];
  for(int i=0; i<nm; i++) {
    d[i] = myalloc(layout->nSites*6*sizeof(QLA_Real));
    PC(packV)(layout, d[i], dest[i]);
  }
  double rsq = resargs[0]->rsqmin;
  //int maxits = invarg->max_iter;
  int maxits = invarg->restart;
  int its = P(solveMulti2)(&sse, d, ms, nm, s, rsq, maxits, "even");
  resargs[0]->final_iter = its;
  free(s);
  for(int i=0; i<nm; i++) {
    PC(unpackV)(layout, dest[i], d[i]);
    free(d[i]);
  }
  //printf("done %s\n", __func__);
}

#endif // HAVE_QLL
