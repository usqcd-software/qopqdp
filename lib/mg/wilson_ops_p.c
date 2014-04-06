//#define DO_TRACE
#include <qop_internal.h>

#define QCDPC(x) QOP_##x
#define QDPN(x) QDP_N_##x

//#define printf0 QOP_printf0
#define printf0(...)

// FIXME: avoid globals
static QDP_DiracFermion *glv=NULL, *glv2=NULL/*, *et1=NULL, *et2=NULL*/;

#define check_glv if(glv==NULL) { glv = QDP_create_D(); glv2 = QDP_create_D(); }
#define check_et if(et1==NULL) { et1 = QDP_create_D(); et2 = QDP_create_D(); }

#define allVN(v) QDP_all_L(QDPN(get_lattice_V)(v))
#define evenVN(v) QDP_even_L(QDPN(get_lattice_V)(v))
#define oddVN(v) QDP_odd_L(QDPN(get_lattice_V)(v))

void
QOP_wilsonDslash(QDP_DiracFermion *out, QDP_DiracFermion *in,
		 QOP_FermionLinksWilson *wil, QLA_Real kappa, int sign,
		 QOP_evenodd_t pout, QOP_evenodd_t pin)
{
  QOP_wilson_dslash_qdp(NULL, wil, kappa, sign, out, in, pout, pin);
}

void
QOP_wilsonDiaginv(QDP_DiracFermion *out, QDP_DiracFermion *in,
		  QOP_FermionLinksWilson *wil, QLA_Real kappa,
		  QOP_evenodd_t pout)
{
  QOP_wilson_diaginv_qdp(NULL, wil, kappa, out, in, pout);
}

void
QOP_wilsonDslashEO(QDP_DiracFermion *out, QDP_DiracFermion *in,
		   QOP_FermionLinksWilson *wil, QLA_Real kappa, int sign,
		   QOP_evenodd_t par)
{
#define NC QDP_get_nc(out)
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
#undef NC
}

// MG operators

#define SETUP if(din==NULL) setup_temps(NCARGVOID)

// FIXME: avoid globals
static QDP_DiracFermion *din=NULL, *dout=NULL;
static QDP_HalfFermion *hin=NULL, *hout=NULL;

static void
setup_temps(NCPROTVOID)
{
  QOP_printf0("nc: %i\n", QLA_Nc);
  din = QDP_create_D();
  dout = QDP_create_D();
  hin = QDP_create_H();
  hout = QDP_create_H();
}

void
QCDPC(wilEoProjectD)(QDP_DiracFermion *ineo, QDP_DiracFermion *in,
		     QCDPC(WilArgs) *w)
{
#define NC QDP_get_nc(ineo)
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
#undef NC
}

void
QCDPC(wilEoReconstructD)(QDP_DiracFermion *out, QDP_DiracFermion *outeo,
			 QDP_DiracFermion *in, QCDPC(WilArgs) *w)
{
#define NC QDP_get_nc(out)
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
#undef NC
}

#ifdef HAVE_NCN

#include <qdp_fn.h>
#include <qdp_dn.h>

void
QCDPC(V1eqD)(QDP_N_ColorVector *v[1], QDP_DiracFermion *d, QDP_Subset sub)
{
#if QDP_Colors == 2 || QDP_Colors == 3
  QDP_D_eq_D((QDP_DiracFermion *)v[0], d, sub);
#else
  QDP_N_V_eq_V(v[0], (QDP_N_ColorVector *)d, sub);
#endif
}

void
QCDPC(DeqV1)(QDP_DiracFermion *d, QDP_N_ColorVector *v[1], QDP_Subset sub)
{
  QDP_D_eq_D(d, (QDP_DiracFermion *)v[0], sub);
}

static void
QCDPC(V1eqDg5)(QDP_N_ColorVector *v[1], QDP_DiracFermion *d, QDP_Subset sub)
{
#if QDP_Colors == 2 || QDP_Colors == 3
  QDP_D_eq_gamma_times_D((QDP_DiracFermion *)v[0], d, 15, sub);
#else
#define NC QDP_get_nc(d)
  SETUP;
  QDP_D_eq_gamma_times_D(dout, d, 15, sub);
  QDP_N_V_eq_V(v[0], (QDP_N_ColorVector *)dout, sub);
#undef NC
#endif
}

static void
QCDPC(DeqV1g5)(QDP_DiracFermion *d, QDP_N_ColorVector *v[1], QDP_Subset sub)
{
  QDP_D_eq_gamma_times_D(d, (QDP_DiracFermion *)v[0], 15, sub);
}

void
QCDPC(V2eqD)(QDP_N_ColorVector *v[2], QDP_DiracFermion *d, QDP_Subset sub)
{
  QDP_DiracFermion *d2[2] = {d,d};
  int g5[2] = {4,4};
  int pm[2] = {1,-1};
#if QDP_Colors == 2 || QDP_Colors == 3
  QDP_H_veq_spproj_D((QDP_HalfFermion **)v, d2, g5, pm, sub, 2);
  //QDP_H_eq_spproj_D((QDP_HalfFermion *)v[0], d, 4, 1, QDP_all);
  //QDP_H_eq_spproj_D((QDP_HalfFermion *)v[1], d, 4, -1, QDP_all);
#else
#define NC QDP_get_nc(d)
  SETUP;
  QDP_HalfFermion *h2[2] = {hout,hin};
  QDP_H_veq_spproj_D(h2, d2, g5, pm, sub, 2);
  QDP_N_V_veq_V(v, (QDP_N_ColorVector **)h2, sub, 2);
#undef NC
#endif
}

void
QCDPC(DeqV2)(QDP_DiracFermion *d, QDP_N_ColorVector *v[2], QDP_Subset sub)
{
  QDP_D_eq_sprecon_H(d, (QDP_HalfFermion *)v[0], 4, 1, sub);
  QDP_D_peq_sprecon_H(d, (QDP_HalfFermion *)v[1], 4, -1, sub);
}

static void
QCDPC(V2eqDg5)(QDP_N_ColorVector *v[2], QDP_DiracFermion *d, QDP_Subset sub)
{
#if QDP_Colors == 2 || QDP_Colors == 3
  QDP_H_eq_spproj_D((QDP_HalfFermion *)v[0], d, 4, 1, sub);
  QDP_H_eqm_spproj_D((QDP_HalfFermion *)v[1], d, 4, -1, sub);
#else
#define NC QDP_get_nc(d)
  SETUP;
  QDP_H_eq_spproj_D(hout, d, 4, 1, sub);
  QDP_N_V_eq_V(v[0], (QDP_N_ColorVector *)hout, sub);
  QDP_H_eqm_spproj_D(hout, d, 4, -1, sub);
  QDP_N_V_eq_V(v[1], (QDP_N_ColorVector *)hout, sub);
#undef NC
#endif
}

static void
QCDPC(DeqV2g5)(QDP_DiracFermion *d, QDP_N_ColorVector *v[2], QDP_Subset sub)
{
  QDP_D_eq_sprecon_H(d, (QDP_HalfFermion *)v[0], 4, 1, sub);
  QDP_D_meq_sprecon_H(d, (QDP_HalfFermion *)v[1], 4, -1, sub);
}

void
QCDPC(wilDV1)(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1],
	      int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1)(din, in, allVN(in[0]));
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, sign, QOP_EVENODD,QOP_EVENODD);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
#undef NC
}

void
QCDPC(wilDV2)(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2],
	      int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/2)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2)(din, in, allVN(in[0]));
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, sign, QOP_EVENODD,QOP_EVENODD);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
#undef NC
}

void
QCDPC(wilPV1)(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1],
	      int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1)(dout, in, allVN(in[0]));
  if(sign>0) {
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_ODD);
    QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVENODD, QOP_EVENODD);
  } else {
    QOP_wilsonDslash(din, dout, w->wil, w->kappa, -1, QOP_EVENODD,QOP_EVENODD);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  }
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
#undef NC
}

void
QCDPC(wilPV2)(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2],
	      int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/2)
  TRACE;
  SETUP;
  TRACE;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2)(dout, in, allVN(in[0]));
  TRACE;
  if(sign>0) {
  TRACE;
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  TRACE;
    QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_ODD);
  TRACE;
    QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_EVENODD, QOP_EVENODD);
  TRACE;
  } else {
    QOP_wilsonDslash(din, dout, w->wil, w->kappa, -1, QOP_EVENODD,QOP_EVENODD);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
    QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  }
  TRACE;
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
  TRACE;
#undef NC
}

void
QCDPC(wilPNEV2)(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2],
		int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/2)
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
#undef NC
}

void
QCDPC(wilEoV1)(QDP_N_ColorVector *out[1], QDP_N_ColorVector *in[1],
	       int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV1g5)(din, in, evenVN(in[0]));
  if(sign>0) {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, -1, QOP_EVEN);
  } else {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, 1, QOP_EVEN);
  }
  QCDPC(V1eqDg5)(out, dout, evenVN(out[0]));
#undef NC
}

void
QCDPC(wilEoV2)(QDP_N_ColorVector *out[2], QDP_N_ColorVector *in[2],
	       int sign, void *args)
{
#define NC (QDP_get_nc(out[0])/2)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *)args;
  QCDPC(DeqV2g5)(din, in, evenVN(in[0]));
  if(sign>0) {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, -1, QOP_EVEN);
  } else {
    QOP_wilsonDslashEO(dout, din, w->wil, w->kappa, 1, QOP_EVEN);
  }
  QCDPC(V2eqDg5)(out, dout, evenVN(out[0]));
#undef NC
}

void
QCDPC(wilEoProjectV1)(QDP_N_ColorVector *ineo[1], QDP_N_ColorVector *in[1],
		      void *args)
{
#define NC (QDP_get_nc(ineo[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(din, in, allVN(in[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, dout, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_even);
  QCDPC(V1eqD)(ineo, dout, evenVN(ineo[0]));
#undef NC
}

void
QCDPC(wilEoReconstructV1)(QDP_N_ColorVector *out[1],
			  QDP_N_ColorVector *outeo[1],
			  QDP_N_ColorVector *in[1], void *args)
{
#define NC (QDP_get_nc(out[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(din, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(din, dout, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV1)(dout, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(din, dout, din, QDP_odd);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
#undef NC
}

// not usual conventions, used in setup for null vectors
void
QCDPC(wilEoReconstructPV1)(QDP_N_ColorVector *out[1],
			   QDP_N_ColorVector *outeo[1],
			   QDP_N_ColorVector *in[1], void *args)
{
#define NC (QDP_get_nc(out[0])/4)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV1)(dout, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV1)(din, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_odd);
  QCDPC(V1eqD)(out, dout, allVN(out[0]));
#undef NC
}

void
QCDPC(wilEoProjectV2)(QDP_N_ColorVector *ineo[2], QDP_N_ColorVector *in[2],
		      void *args)
{
#define NC (QDP_get_nc(ineo[0])/2)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(din, in, allVN(in[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QOP_wilsonDslash(dout, dout, w->wil, w->kappa, 1, QOP_EVEN, QOP_ODD);
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_even);
  QCDPC(V2eqD)(ineo, dout, evenVN(ineo[0]));
#undef NC
}

void
QCDPC(wilEoReconstructV2)(QDP_N_ColorVector *out[2],
			  QDP_N_ColorVector *outeo[2],
			  QDP_N_ColorVector *in[2], void *args)
{
#define NC (QDP_get_nc(out[0])/2)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(din, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(din, dout, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV2)(dout, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(din, dout, din, QDP_odd);
  QOP_wilsonDiaginv(dout, din, w->wil, w->kappa, QOP_ODD);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
#undef NC
}

// not usual conventions, used in setup for null vectors
void
QCDPC(wilEoReconstructPV2)(QDP_N_ColorVector *out[2],
			   QDP_N_ColorVector *outeo[2],
			   QDP_N_ColorVector *in[2], void *args)
{
#define NC (QDP_get_nc(out[0])/2)
  SETUP;
  QCDPC(WilArgs) *w = (QCDPC(WilArgs) *) args;
  QCDPC(DeqV2)(dout, outeo, evenVN(outeo[0]));
  QOP_wilsonDiaginv(din, dout, w->wil, w->kappa, QOP_EVEN);
  QOP_wilsonDslash(dout, din, w->wil, w->kappa, 1, QOP_ODD, QOP_EVEN);
  QCDPC(DeqV2)(din, in, oddVN(in[0]));
  QDP_D_eq_D_minus_D(dout, din, dout, QDP_odd);
  QCDPC(V2eqD)(out, dout, allVN(out[0]));
#undef NC
}

#endif // HAVE_NCN
