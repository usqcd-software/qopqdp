#ifndef _QOP_INTERNALP_H
#define _QOP_INTERNALP_H

#if QOP_Precision == 1
#define QOPP(x) QOP_F_##x
#define QOPPC(x) QOP_F3_##x
#define QDPP(x) QDP_F_##x
#define QDPPC(x) QDP_F3_##x
#define REAL float
#else
#define QOPP(x) QOP_D_##x
#define QOPPC(x) QOP_D3_##x
#define QDPP(x) QDP_D_##x
#define QDPPC(x) QDP_D3_##x
#define REAL double
#endif

#define QLATYPE_V QLA_ColorVector
#define QLATYPE_D QLA_DiracFermion
#define QLATYPE_M QLA_ColorMatrix

#define QOP_qdp_eq_raw(abbr, qdp, raw, evenodd) \
  QDPPC(insert_##abbr)(qdp, (QLATYPE_##abbr *)raw, QDP_all)

#define QOP_raw_eq_qdp(abbr, raw, qdp, evenodd) \
  QDP_extract_##abbr((QLATYPE_##abbr *)raw, qdp, QDP_all)

typedef struct {
  int tmp;
} QOPPC(common_t);
extern QOPPC(common_t) QOPPC(common);

struct QOPPC(ColorVector_struct) {
  QDPPC(ColorVector) *cv;
  REAL *raw;
};

struct QOPPC(DiracFermion_struct) {
  QDPPC(DiracFermion) *df;
  REAL *raw;
};

struct QOPPC(GaugeField_struct) {
  QDPPC(ColorMatrix) **links;
  REAL **raw;
};

struct QOPPC(Force_struct) {
  QDPPC(ColorMatrix) **force;
  REAL **raw;
};

  /* Asqtad datatypes */

struct QOPPC(FermionLinksAsqtad_struct) {
  int dblstored;
  QDPPC(ColorMatrix) **fatlinks;
  QDPPC(ColorMatrix) **longlinks;
  QDPPC(ColorMatrix) **fwdlinks;
  QDPPC(ColorMatrix) **bcklinks;
  QDPPC(ColorMatrix) **dbllinks;
};

  /* Wilson datatypes */

struct QOPPC(FermionLinksWilson_struct) {
  QDPPC(ColorMatrix) **links;
};

/* internal routines */

void QOPPC(qdpM_eq_raw)(QDP_ColorMatrix *cm, REAL *lnk);

#endif /* _QOP_INTERNALP_H */
