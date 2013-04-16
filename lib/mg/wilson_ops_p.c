//#define DO_TRACE
#include <qop_internal.h>
#include <qdp_fn.h>
#include <qdp_dn.h>

#if QOP_Precision == 'F'
#define QCDPC(x) QOP_F3_##x
#define QDPN(x) QDP_FN_##x
#else
#define QCDPC(x) QOP_D3_##x
#define QDPN(x) QDP_DN_##x
#endif

//#define printf0 QOP_printf0
#define printf0(...)

static QDP_DiracFermion *glv=NULL, *glv2=NULL/*, *et1=NULL, *et2=NULL*/;

#define check_glv if(glv==NULL) { glv = QDP_create_D(); glv2 = QDP_create_D(); }
#define check_et if(et1==NULL) { et1 = QDP_create_D(); et2 = QDP_create_D(); }

#define allVN(v) QDP_all_L(QDPN(get_lattice_V)(v))
#define evenVN(v) QDP_even_L(QDPN(get_lattice_V)(v))
#define oddVN(v) QDP_odd_L(QDPN(get_lattice_V)(v))
#define allD3(v) QDP_all_L(QDPN(get_lattice_D)(v))
#define evenD3(v) QDP_even_L(QDPN(get_lattice_D)(v))
#define oddD3(v) QDP_odd_L(QDPN(get_lattice_D)(v))

void
QOP_wilsonDslash(QDP_DiracFermion *out, QDP_DiracFermion *in, QOP_FermionLinksWilson *wil,
		 QLA_Real kappa, int sign, QOP_evenodd_t pout, QOP_evenodd_t pin)
{
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, out, in, pout, pin);
}

void
QOP_wilsonDiaginv(QDP_DiracFermion *out, QDP_DiracFermion *in, QOP_FermionLinksWilson *wil,
		  QLA_Real kappa, QOP_evenodd_t pout)
{
  QOP_wilson_diaginv_qdp(NULL, wil, kappa, out, in, pout);
}

void
QOP_wilsonDslashEO(QDP_DiracFermion *out, QDP_DiracFermion *in,
		   QOP_FermionLinksWilson *wil, QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  check_glv;
  QOP_evenodd_t paro = QOP_ODD;
  QDP_Subset sub = QDP_even;
  QDP_DiracFermion *tv = glv;
  QDP_DiracFermion *tvo = glv2;
  if(par == QOP_ODD) {
    paro = QOP_EVEN;
    sub = QDP_odd;
    tv = glv2;
    tvo = glv;
  }
  if(sign>0) {
    QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, tv, in, paro, par);
    QOP_wilson_diaginv_qdp(NULL, wil, kappa, tvo, tv, paro);
    QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, out, tvo, par, paro);
    QOP_wilson_diaginv_qdp(NULL, wil, kappa, tvo, out, par);
    QDP_D_eq_D_minus_D(out, in, tvo, sub);
  } else {
    QOP_wilson_diaginv_qdp(NULL, wil, kappa, tvo, in, par);
    QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, tv, tvo, paro, par);
    QOP_wilson_diaginv_qdp(NULL, wil, kappa, tvo, tv, paro);
    QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, out, tvo, par, paro);
    QDP_D_eq_D_minus_D(out, in, out, sub);
  }
}

void
QOP_wilsonDslashEOS(QDP_DiracFermion *out, QDP_DiracFermion *in,
		    QOP_FermionLinksWilson *wil, QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  check_glv;
  QOP_evenodd_t paro = QOP_ODD;
  QDP_Subset sub = QDP_even;
  QDP_DiracFermion *tv = glv;
  QDP_DiracFermion *tvo = glv2;
  if(par == QOP_ODD) {
    paro = QOP_EVEN;
    sub = QDP_odd;
    tv = glv2;
    tvo = glv;
  }
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, tv, in, paro, par);
  QOP_wilson_diaginv_qdp(NULL, wil, kappa, tvo, tv, paro);
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, tv, tvo, par, paro);
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, out, in, par, par);
  QDP_D_meq_D(out, tv, sub);
}

void
QOP_wilsonDslashEOH(QDP_DiracFermion *out, QDP_DiracFermion *in,
		    QOP_FermionLinksWilson *wil, QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  QDP_Subset sub = QDP_even;
  if(par == QOP_ODD) sub = QDP_odd;
  QOP_wilsonDslashEOS(out, in, wil, kappa, sign, par);
  QDP_D_eq_gamma_times_D(out, out, 15, sub);
}

#if 0
void
QOP_wilsonDslashVec(vec *out, vec *in, QOP_FermionLinksWilson *wil,
		    QLA_Real kappa, int sign, QOP_evenodd_t pout, QOP_evenodd_t pin)
{
  check_et;
  QDP_Subset sout = (pout==QOP_EVEN) ? QDP_even : (pout==QOP_ALL) ? QDP_all : QDP_odd;
  QDP_Subset sin = (pin==QOP_EVEN) ? QDP_even : (pin==QOP_ALL) ? QDP_all : QDP_odd;
  QDP_insert_packed_D(et1, (QLA_DiracFermion *)(in->data), sin);
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, et2, et1, pout, pin);
  QDP_extract_packed_D((QLA_DiracFermion *)(out->data), et2, sout);
}

void
QOP_wilsonDslashEOVec(vec *out, vec *in, QOP_FermionLinksWilson *wil,
		      QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  check_et;
  QDP_Subset sub = (par==QOP_EVEN) ? QDP_even : QDP_odd;
  QDP_insert_packed_D(et1, (QLA_DiracFermion *)(in->data), sub);
  QOP_wilsonDslashEO(et2, et1, wil, kappa, sign, par);
  QDP_extract_packed_D((QLA_DiracFermion *)(out->data), et2, sub);
}

void
QOP_wilsonDslashEOSVec(vec *out, vec *in, QOP_FermionLinksWilson *wil,
		       QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  check_et;
  QDP_Subset sub = (par==QOP_EVEN) ? QDP_even : QDP_odd;
  QDP_insert_packed_D(et1, (QLA_DiracFermion *)(in->data), sub);
  QOP_wilsonDslashEOS(et2, et1, wil, kappa, sign, par);
  QDP_extract_packed_D((QLA_DiracFermion *)(out->data), et2, sub);
}

void
QOP_wilsonDslashEOHVec(vec *out, vec *in, QOP_FermionLinksWilson *wil,
		       QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  check_et;
  QDP_Subset sub = (par==QOP_EVEN) ? QDP_even : QDP_odd;
  QDP_insert_packed_D(et1, (QLA_DiracFermion *)(in->data), sub);
  QOP_wilsonDslashEOH(et2, et1, wil, kappa, sign, par);
  QDP_extract_packed_D((QLA_DiracFermion *)(out->data), et2, sub);
}
#endif

static int h2func(int x[], void *args)
{
  int k=0;
  for(int i=3; i>=0; i--) k = 4*k + x[i]%4;
  return k;
}

void
QOP_getDiagEOS(QDP_DiracPropagator *p, QOP_FermionLinksWilson *wil,
	       QLA_Real kappa, int sign, QOP_evenodd_t par)
{
  QDP_DiracFermion *in = QDP_create_D();
  QDP_DiracFermion *out = QDP_create_D();
  QDP_Subset *sub;
  QDP_Subset eosub = QDP_even;
  if(par == QOP_ODD) eosub = QDP_odd;

  int ns = 256;
  sub = QDP_create_subset(h2func, NULL, 0, ns);

  QDP_P_eq_zero(p, eosub);
  for(int s=0; s<ns; s++) {
    for(int i=0; i<12; i++) {
      int ic = i%3;
      int is = i/3;
      QLA_DiracFermion d;
      QLA_D_eq_zero(&d);
      QLA_c_eq_r(QLA_elem_D(d, ic, is), 1);
      QDP_D_eq_zero(in, eosub);
      QDP_D_eq_d(in, &d, sub[s]);
      QOP_wilsonDslashEOS(out, in, wil, kappa, sign, par);
      QDP_P_eq_diracvec_D(p, out, ic, is, sub[s]);
    }
  }

  QDP_destroy_subset(sub);
  QDP_destroy_D(in);
  QDP_destroy_D(out);
}

// MG operators

void
QCDPC(V1eqD)(QDPN(ColorVector) *v[1], QDP_DiracFermion *d, QDP_Subset sub)
{
  QDP_D_eq_D((QDP_DiracFermion *)v[0], d, sub);
}

void
QCDPC(DeqV1)(QDP_DiracFermion *d, QDPN(ColorVector) *v[1], QDP_Subset sub)
{
  QDP_D_eq_D(d, (QDP_DiracFermion *)v[0], sub);
}

void
QCDPC(V1eqDg5)(QDPN(ColorVector) *v[1], QDP_DiracFermion *d, QDP_Subset sub)
{
  QDP_D_eq_gamma_times_D((QDP_DiracFermion *)v[0], d, 15 ,sub);
}

void
QCDPC(DeqV1g5)(QDP_DiracFermion *d, QDPN(ColorVector) *v[1], QDP_Subset sub)
{
  QDP_D_eq_gamma_times_D(d, (QDP_DiracFermion *)v[0], 15, sub);
}

void
QCDPC(V2eqD)(QDPN(ColorVector) *v[2], QDP_DiracFermion *d, QDP_Subset sub)
{
  int g5[2] = {4,4};
  int pm[2] = {1,-1};
  QDP_DiracFermion *d2[2] = {d,d};
  QDP_H_veq_spproj_D((QDP_HalfFermion **)v, d2, g5, pm, sub, 2);
  //QDP_H_eq_spproj_D((QDP_HalfFermion *)v[0], d, 4, 1, QDP_all);
  //QDP_H_eq_spproj_D((QDP_HalfFermion *)v[1], d, 4, -1, QDP_all);
}

void
QCDPC(DeqV2)(QDP_DiracFermion *d, QDPN(ColorVector) *v[2], QDP_Subset sub)
{
  //int g5[2] = {4,4};
  //int pm[2] = {1,-1};
  //QDP_DiracFermion *d2[2] = {d,d};
  //QDP_D_veq_sprecon_H(d2, (QDP_HalfFermion **)v, g5, pm, QDP_all, 2);
  QDP_D_eq_sprecon_H(d, (QDP_HalfFermion *)v[0], 4, 1, sub);
  QDP_D_peq_sprecon_H(d, (QDP_HalfFermion *)v[1], 4, -1, sub);
}

void
QCDPC(V2eqDg5)(QDPN(ColorVector) *v[2], QDP_DiracFermion *d, QDP_Subset sub)
{
  //int g5[2] = {4,4};
  //int pm[2] = {1,-1};
  //QDP_DiracFermion *d2[2] = {d,d};
  //QDP_H_veq_spproj_D((QDP_HalfFermion **)v, d2, g5, pm, QDP_all, 2);
  QDP_H_eq_spproj_D((QDP_HalfFermion *)v[0], d, 4, 1, sub);
  QDP_H_eqm_spproj_D((QDP_HalfFermion *)v[1], d, 4, -1, sub);
}

void
QCDPC(DeqV2g5)(QDP_DiracFermion *d, QDPN(ColorVector) *v[2], QDP_Subset sub)
{
  //int g5[2] = {4,4};
  //int pm[2] = {1,-1};
  //QDP_DiracFermion *d2[2] = {d,d};
  //QDP_D_veq_sprecon_H(d2, (QDP_HalfFermion **)v, g5, pm, QDP_all, 2);
  QDP_D_eq_sprecon_H(d, (QDP_HalfFermion *)v[0], 4, 1, sub);
  QDP_D_meq_sprecon_H(d, (QDP_HalfFermion *)v[1], 4, -1, sub);
}

#define SETUP if(din==NULL) setup_temps()

static QDP_DiracFermion *din=NULL, *dout=NULL;

static void
setup_temps()
{
  din = QDP_create_D();
  dout = QDP_create_D();
}

void
QCDPC(wilDV1)(QDPN(ColorVector) *out[1], QDPN(ColorVector) *in[1],
	      int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1)(din, in, allVN(in[0]));
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, sign, QOP_EVENODD, QOP_EVENODD);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilDV2)(QDPN(ColorVector) *out[2], QDPN(ColorVector) *in[2],
	      int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2)(din, in, allVN(in[0]));
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, sign, QOP_EVENODD, QOP_EVENODD);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilPV1)(QDPN(ColorVector) *out[1], QDPN(ColorVector) *in[1],
	      int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1)(dout, in, allVN(in[0]));
  if(sign>0) {
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_ODD);
    QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVENODD, QOP_EVENODD);
  } else {
    QOP_wilsonDslash(din, dout, w->wil, w->kappa, -1, QOP_EVENODD, QOP_EVENODD);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  }
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilPV2)(QDPN(ColorVector) *out[2], QDPN(ColorVector) *in[2],
	      int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2)(dout, in, allVN(in[0]));
  if(sign>0) {
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_ODD);
    QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVENODD, QOP_EVENODD);
  } else {
    QOP_wilsonDslash(din, dout, w->wil, w->kappa, -1, QOP_EVENODD, QOP_EVENODD);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  }
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilPNEV2)(QDPN(ColorVector) *out[2], QDPN(ColorVector) *in[2], int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2)(dout, in, allVN(in[0]));
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVENODD, QOP_EVENODD);
  QDP_D_eq_D(din, dout, QDP_all);
  QOP_wilsonDslash(din, dout, w->wil, w->kappa, -1, QOP_EVENODD, QOP_EVENODD);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilEoV1)(QDPN(ColorVector) *out[1], QDPN(ColorVector) *in[1], int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1g5)(din, in, evenVN(in[0]));
  if(sign>0) {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, -1, QOP_EVEN);
  } else {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, 1, QOP_EVEN);
  }
  QCDPC(V1eqDg5)(out, dout, evenVN(out[0]));
}

void
QCDPC(wilEoV2)(QDPN(ColorVector) *out[2], QDPN(ColorVector) *in[2], int sign, void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2g5)(din, in, evenVN(in[0]));
  if(sign>0) {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, -1, QOP_EVEN);
  } else {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, 1, QOP_EVEN);
  }
  QCDPC(V2eqDg5)(out, dout, evenVN(out[0]));
}

void
QCDPC(wilEoProjectD)(QDPPC(DiracFermion) *ineo, QDPPC(DiracFermion) *in, QCDPC(WilArgs) *w)
{
  SETUP;
#if 0
  QOP_wilsonDiaginv(din, in, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(dout, in, dout, QDP_even);
  QOP_wilsonDiaginv(ineo, dout, w->wil, w->kappa, QOP_EVEN);
#else
  QOP_wilsonDiaginv(din, in, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(ineo, din, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(ineo, in, ineo, QDP_even);
#endif
}

void
QCDPC(wilEoReconstructD)(QDPPC(DiracFermion) *out, QDPPC(DiracFermion) *outeo,
			 QDPPC(DiracFermion) *in, QCDPC(WilArgs) *w)
{
  SETUP;
#if 0
  QOP_wilsonDslash(dout, out, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QDP_D_eq_D_minus_D(dout, in, dout, QDP_odd);
  QOP_wilsonDiaginv(out, dout, w->wil, w->kappa, QOP_ODD);
#else
  QOP_wilsonDiaginv(out, outeo, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(dout, out, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QDP_D_eq_D_minus_D(dout, in, dout, QDP_odd);
  QOP_wilsonDiaginv(out, dout, w->wil, w->kappa, QOP_ODD);
#endif
}

void
QCDPC(wilEoProjectV1)(QDPN(ColorVector) *ineo[1], QDPN(ColorVector) *in[1],
		      void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(din, in, allVN(in[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, dout, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_even);
  QCDPC(V1eqD)(ineo, dout, evenVN(ineo[0]));
}

void
QCDPC(wilEoReconstructV1)(QDPN(ColorVector) *out[1],
			  QDPN(ColorVector) *outeo[1],
			  QDPN(ColorVector) *in[1], void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(din, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(din, dout, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV1)(dout, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(din, dout, din, QDP_odd);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
}

// not usual conventions, used in setup for null vectors
void
QCDPC(wilEoReconstructPV1)(QDPN(ColorVector) *out[1],
			   QDPN(ColorVector) *outeo[1],
			   QDPN(ColorVector) *in[1], void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(dout, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV1)(din, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_odd);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
}

void
QCDPC(wilEoProjectV2)(QDPN(ColorVector) *ineo[2], QDPN(ColorVector) *in[2],
		      void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(din, in, allVN(in[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, dout, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_even);
  QCDPC(V2eqD)(ineo, dout, evenVN(ineo[0]));
}

void
QCDPC(wilEoReconstructV2)(QDPN(ColorVector) *out[2],
			  QDPN(ColorVector) *outeo[2],
			  QDPN(ColorVector) *in[2], void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(din, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(din, dout, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV2)(dout, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(din, dout, din, QDP_odd);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
}

// not usual conventions, used in setup for null vectors
void
QCDPC(wilEoReconstructPV2)(QDPN(ColorVector) *out[2],
			   QDPN(ColorVector) *outeo[2],
			   QDPN(ColorVector) *in[2], void *args)
{
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(dout, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV2)(din, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_odd);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
}
