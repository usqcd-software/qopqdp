#include <stdlib.h>
#include <sys/time.h>
#include <qop_internal.h>


/*******************/
/* public routines */
/*******************/

/* raw routines */

QOPPC(ColorVector) *
QOPPC(create_V_from_raw)(REAL *src, QOP_evenodd_t evenodd)
{
  QOPPC(ColorVector) *qopcv;
  QOP_malloc(qopcv, QOPPC(ColorVector), 1);
  qopcv->cv = QDPPC(create_V)();
  QOP_qdp_eq_raw(V, qopcv->cv, src, evenodd);
  qopcv->raw = NULL;
  return qopcv;
}

QOPPC(DiracFermion) *
QOPPC(create_D_from_raw)(REAL *src, QOP_evenodd_t evenodd)
{
  QOPPC(DiracFermion) *qopdf;
  QOP_malloc(qopdf, QOPPC(DiracFermion), 1);
  qopdf->df = QDPPC(create_D)();
  QOP_qdp_eq_raw(D, qopdf->df, src, evenodd);
  qopdf->raw = NULL;
  return qopdf;
}

QOPPC(GaugeField) *
QOPPC(create_G_from_raw)(REAL *links[], QOP_evenodd_t evenodd)
{
  QOPPC(GaugeField) *qopgf;
  int i;
  QOP_malloc(qopgf, QOPPC(GaugeField), 1);
  QOP_malloc(qopgf->links, QDPPC(ColorMatrix) *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDPPC(create_M)();
    QOP_qdp_eq_raw(M, qopgf->links[i], links[i], evenodd);
  }
  qopgf->raw = NULL;
  return qopgf;
}

QOPPC(Force) *
QOPPC(create_F_from_raw)(REAL *force[], QOP_evenodd_t evenodd)
{
  QOPPC(Force) *qopf;
  int i;
  QOP_malloc(qopf, QOPPC(Force), 1);
  QOP_malloc(qopf->force, QDPPC(ColorMatrix) *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDPPC(create_M)();
    QOP_qdp_eq_raw(M, qopf->force[i], force[i], evenodd);
  }
  qopf->raw = NULL;
  return qopf;
}

void
QOPPC(extract_V_to_raw)(REAL *dest, QOPPC(ColorVector) *src,
			QOP_evenodd_t evenodd)
{
  QOP_raw_eq_qdp(V, dest, src->cv, evenodd);
}

void
QOPPC(extract_D_to_raw)(REAL *dest, QOPPC(DiracFermion) *src,
			QOP_evenodd_t evenodd)
{
  QOP_raw_eq_qdp(D, dest, src->df, evenodd);
}

void
QOPPC(extract_G_to_raw)(REAL *dest[], QOPPC(GaugeField) *src,
			QOP_evenodd_t evenodd)
{
  int i;
  for(i=0; i<QOP_common.ndim; i++) {
    QOP_raw_eq_qdp(M, dest[i], src->links[i], evenodd);
  }
}

void
QOPPC(extract_F_to_raw)(REAL *dest[], QOPPC(Force) *src, QOP_evenodd_t evenodd)
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
QOP_convert_V_from_raw(REAL *src, QOP_evenodd_t evenodd)
{
  QOP_ColorVector *ret =
    QOP_create_V_from_raw(src, evenodd);
  ret->raw = src;
  return ret;
}

QOP_DiracFermion *
QOP_convert_D_from_raw(REAL *src, QOP_evenodd_t evenodd)
{
  QOP_DiracFermion *ret =
    QOP_create_D_from_raw(src, evenodd);
  ret->raw = src;
  return ret;
}

QOP_GaugeField *
QOP_convert_G_from_raw(REAL *links[], QOP_evenodd_t evenodd)
{
  QOP_GaugeField *ret =
    QOP_create_G_from_raw(links, evenodd);
  ret->raw = links;
  return ret;
}

QOP_Force *
QOP_convert_F_from_raw(REAL *force[], QOP_evenodd_t evenodd)
{
  QOP_Force *ret =
    QOP_create_F_from_raw(force, evenodd);
  ret->raw = force;
  return ret;
}

REAL *
QOP_convert_V_to_raw(QOP_ColorVector *src, QOP_evenodd_t evenodd)
{
  REAL *ret = src->raw;
  if(!ret) QOP_malloc(ret, REAL, QOP_sites_on_node_raw_V(evenodd));
  QOP_extract_V_to_raw(ret, src, evenodd);
  QOP_destroy_V(src);
  return ret;
}

REAL *
QOP_convert_D_to_raw(QOP_DiracFermion *src, QOP_evenodd_t evenodd)
{
  REAL *ret = src->raw;
  if(!ret) QOP_malloc(ret, REAL, QOP_sites_on_node_raw_D(evenodd));
  QOP_extract_D_to_raw(ret, src, evenodd);
  QOP_destroy_D(src);
  return ret;
}

REAL **
QOP_convert_G_to_raw(QOP_GaugeField *src, QOP_evenodd_t evenodd)
{
  REAL **ret = src->raw;
  if(!ret) {
    int i;
    QOP_malloc(ret, REAL *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(ret[i], REAL, QOP_sites_on_node_raw_G(evenodd));
    }
  }
  QOP_extract_G_to_raw(ret, src, evenodd);
  QOP_destroy_G(src);
  return ret;
}

REAL **
QOP_convert_F_to_raw(QOP_Force *src, QOP_evenodd_t evenodd)
{
  REAL **ret = src->raw;
  if(!ret) {
    int i;
    QOP_malloc(ret, REAL *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(ret[i], REAL, QOP_sites_on_node_raw_F(evenodd));
    }
  }
  QOP_extract_F_to_raw(ret, src, evenodd);
  QOP_destroy_F(src);
  return ret;
}


/* qdp routines */

QOPPC(ColorVector) *
QOPPC(create_V_from_qdp)(QDP_ColorVector *src)
{
  QOPPC(ColorVector) *qopcv;
  QOP_malloc(qopcv, QOPPC(ColorVector), 1);
  qopcv->cv = QDP_create_V();
  QDP_V_eq_V(qopcv->cv, src, QDP_all);
  return qopcv;
}

QOPPC(DiracFermion) *
QOPPC(create_D_from_qdp)(QDP_DiracFermion *src)
{
  QOPPC(DiracFermion) *qopdf;
  QOP_malloc(qopdf, QOPPC(DiracFermion), 1);
  qopdf->df = QDP_create_D();
  QDP_D_eq_D(qopdf->df, src, QDP_all);
  return qopdf;
}

QOPPC(GaugeField) *
QOPPC(create_G_from_qdp)(QDP_ColorMatrix *links[])
{
  QOPPC(GaugeField) *qopgf;
  int i;
  QOP_malloc(qopgf, QOPPC(GaugeField), 1);
  QOP_malloc(qopgf->links, QDPPC(ColorMatrix) *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = QDPPC(create_M)();
    QDP_M_eq_M(qopgf->links[i], links[i], QDP_all);
  }
  qopgf->raw = NULL;
  return qopgf;
}

QOPPC(Force) *
QOPPC(create_F_from_qdp)(QDP_ColorMatrix *force[])
{
  QOPPC(Force) *qopf;
  int i;
  QOP_malloc(qopf, QOPPC(Force), 1);
  QOP_malloc(qopf->force, QDPPC(ColorMatrix) *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopf->force[i] = QDPPC(create_M)();
    QDP_M_eq_M(qopf->force[i], force[i], QDP_all);
  }
  qopf->raw = NULL;
  return qopf;
}

QOPPC(GaugeField) *
QOPPC(convert_G_from_qdp)(QDP_ColorMatrix *links[])
{
  QOPPC(GaugeField) *qopgf;
  int i;
  QOP_malloc(qopgf, QOPPC(GaugeField), 1);
  QOP_malloc(qopgf->links, QDPPC(ColorMatrix) *, QOP_common.ndim);
  for(i=0; i<QOP_common.ndim; i++) {
    qopgf->links[i] = links[i];
  }
  qopgf->raw = NULL;
  return qopgf;
}

void
QOPPC(extract_V_to_qdp)(QDP_ColorVector *dest, QOPPC(ColorVector) *src)
{
  QDP_V_eq_V(dest, src->cv, QDP_all);
}

void
QOPPC(extract_D_to_qdp)(QDP_DiracFermion *dest, QOPPC(DiracFermion) *src)
{
  QDP_D_eq_D(dest, src->df, QDP_all);
}
