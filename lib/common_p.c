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

QOP_ColorVector *
QOP_create_V_from_raw(QOP_Real *src, QOP_evenodd_t evenodd)
{
  QOP_ColorVector *qopcv;
  QOP_malloc(qopcv, QOP_ColorVector, 1);
  qopcv->cv = QDP_create_V();
  QOP_qdp_eq_raw(V, qopcv->cv, src, evenodd);
  qopcv->raw = NULL;
  return qopcv;
}

QOP_DiracFermion *
QOP_create_D_from_raw(QOP_Real *src, QOP_evenodd_t evenodd)
{
  QOP_DiracFermion *qopdf;
  QOP_malloc(qopdf, QOP_DiracFermion, 1);
  qopdf->df = QDP_create_D();
  QOP_qdp_eq_raw(D, qopdf->df, src, evenodd);
  qopdf->raw = NULL;
  return qopdf;
}

QOP_GaugeField *
QOP_create_G_from_raw(QOP_Real *links[], QOP_evenodd_t evenodd)
{
  QOP_GaugeField *qopgf = QOP_new_G();
  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDP_create_M();
    QOP_qdp_eq_raw(M, qopgf->links[i], links[i], evenodd);
  }
  return qopgf;
}

QOP_Force *
QOP_create_F_from_raw(QOP_Real *force[], QOP_evenodd_t evenodd)
{
  QOP_Force *qopf;
  int i;
  QOP_malloc(qopf, QOP_Force, 1);
  QOP_malloc(qopf->force, QDP_ColorMatrix *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDP_create_M();
    QOP_qdp_eq_raw(M, qopf->force[i], force[i], evenodd);
  }
  qopf->raw = NULL;
  return qopf;
}

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

QOP_ColorVector *
QOP_convert_V_from_raw(QOP_Real *src, QOP_evenodd_t evenodd)
{
  QOP_ColorVector *ret =
    QOP_create_V_from_raw(src, evenodd);
  ret->raw = src;
  return ret;
}

QOP_DiracFermion *
QOP_convert_D_from_raw(QOP_Real *src, QOP_evenodd_t evenodd)
{
  QOP_DiracFermion *ret =
    QOP_create_D_from_raw(src, evenodd);
  ret->raw = src;
  return ret;
}

QOP_GaugeField *
QOP_convert_G_from_raw(QOP_Real *links[], QOP_evenodd_t evenodd)
{
  QOP_GaugeField *ret =
    QOP_create_G_from_raw(links, evenodd);
  ret->raw = links;
  return ret;
}

QOP_Force *
QOP_convert_F_from_raw(QOP_Real *force[], QOP_evenodd_t evenodd)
{
  QOP_Force *ret =
    QOP_create_F_from_raw(force, evenodd);
  ret->raw = force;
  return ret;
}

QOP_Real *
QOP_convert_V_to_raw(QOP_ColorVector *src, QOP_evenodd_t evenodd)
{
  QOP_Real *ret = src->raw;
  if(!ret) QOP_malloc(ret, QOP_Real, QOP_sites_on_node_raw_V(evenodd));
  QOP_extract_V_to_raw(ret, src, evenodd);
  QOP_destroy_V(src);
  return ret;
}

QOP_Real *
QOP_convert_D_to_raw(QOP_DiracFermion *src, QOP_evenodd_t evenodd)
{
  QOP_Real *ret = src->raw;
  if(!ret) QOP_malloc(ret, QOP_Real, QOP_sites_on_node_raw_D(evenodd));
  QOP_extract_D_to_raw(ret, src, evenodd);
  QOP_destroy_D(src);
  return ret;
}

QOP_Real **
QOP_convert_G_to_raw(QOP_GaugeField *src, QOP_evenodd_t evenodd)
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
QOP_convert_F_to_raw(QOP_Force *src, QOP_evenodd_t evenodd)
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

QOP_ColorVector *
QOP_create_V_from_qdp(QDP_ColorVector *src)
{
  QOP_ColorVector *qopcv;
  QOP_malloc(qopcv, QOP_ColorVector, 1);
  qopcv->cv = QDP_create_V();
  QDP_V_eq_V(qopcv->cv, src, QDP_all);
  qopcv->raw = NULL;
  return qopcv;
}

QOP_DiracFermion *
QOP_create_D_from_qdp(QDP_DiracFermion *src)
{
  QOP_DiracFermion *qopdf;
  QOP_malloc(qopdf, QOP_DiracFermion, 1);
  qopdf->df = QDP_create_D();
  QDP_D_eq_D(qopdf->df, src, QDP_all);
  qopdf->raw = NULL;
  return qopdf;
}

QOP_GaugeField *
QOP_create_G_from_qdp(QDP_ColorMatrix *links[])
{
  QOP_GaugeField *qopgf = QOP_new_G();
  for(int i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDP_create_M();
    QDP_M_eq_M(qopgf->links[i], links[i], QDP_all);
  }
  return qopgf;
}

QOP_Force *
QOP_create_F_from_qdp(QDP_ColorMatrix *force[])
{
  QOP_Force *qopf;
  int i;
  QOP_malloc(qopf, QOP_Force, 1);
  QOP_malloc(qopf->force, QDP_ColorMatrix *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDP_create_M();
    QDP_M_eq_M(qopf->force[i], force[i], QDP_all);
  }
  qopf->raw = NULL;
  return qopf;
}

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

static void
rephase_func(QLA_ColorMatrix *m, int coords[])
{
  if(bc_dir>=0) {
    int rshift = (coords[bc_dir] + bc_coord - bc_origin) % bc_coord;
    if(rshift == bc_coord-1) {
      QLA_ColorMatrix t;
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
    QDP_M_eq_func(l[i], rephase_func, QDP_all);
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
