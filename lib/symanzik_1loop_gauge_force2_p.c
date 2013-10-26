/* gauge force for Symanzik improved 1x1 + 1x2 + 1x1x1 */
//#define DO_TRACE
#include <qop_internal.h>

#define NMTMP 6
#define NFTMP 4
#define NBTMP 2
typedef struct {
  QDP_ColorMatrix **mtmp, **ftmp0, **ftmp[NFTMP], **btmp0, **btmp[NBTMP];
  QDP_Lattice *lat;
  QDP_Subset sub;
  int nc, nd;
} tmpstruct;

#define PEQM (2*NC*NC)
#define EQMTM (NC*NC*(8*NC-2))
#define PEQMTM (NC*NC*(8*NC))

static void
set_temps(tmpstruct *t)
{
#define NC t->nc
  for(int i=0; i<NMTMP; i++) t->mtmp[i] = QDP_create_M_L(t->lat);
  for(int i=0; i<NFTMP; i++) {
    t->ftmp0[i] = QDP_create_M_L(t->lat);
    for(int j=0; j<t->nd; j++) t->ftmp[i][j] = QDP_create_M_L(t->lat);
  }
  for(int i=0; i<NBTMP; i++) {
    t->btmp0[i] = QDP_create_M_L(t->lat);
    for(int j=0; j<t->nd; j++) t->btmp[i][j] = QDP_create_M_L(t->lat);
  }
#undef NC
}

static void
free_temps(tmpstruct *t)
{
  for(int i=0; i<NMTMP; i++) QDP_destroy_M(t->mtmp[i]);
  for(int i=0; i<NFTMP; i++) {
    QDP_destroy_M(t->ftmp0[i]);
    for(int j=0; j<t->nd; j++) QDP_destroy_M(t->ftmp[i][j]);
  }
  for(int i=0; i<NBTMP; i++) {
    QDP_destroy_M(t->btmp0[i]);
    for(int j=0; j<t->nd; j++) QDP_destroy_M(t->btmp[i][j]);
  }
}

#if 0
static void
staplefb(QDP_ColorMatrix *outf, QDP_ColorMatrix *outb, QDP_ColorMatrix *in0,
	 QDP_ColorMatrix *link0, int mu, int nu)
{
#define link     ftmp0[0]
#define linkmu   ftmp[0][mu]
#define in       ftmp0[1]
#define innu     ftmp[1][nu]
#define linkin   mtmp[0]
#define back     btmp0[0]
#define backnu   btmp[0][nu]
#define linkinnu mtmp[1]

  // fmu += Unu * Umufnu * Unufmu+
  // fmu += Unubnu+ * Umubnu * Unufmubnu
  QDP_M_eq_M(link, link0, t->sub);
  QDP_M_eq_sM(linkmu, link, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(in, in0, t->sub);
  QDP_M_eq_sM(innu, in, QDP_neighbor[nu], QDP_forward, t->sub);
  QDP_M_eq_Ma_times_M(linkin, link, in, t->sub);
  QDP_M_eq_M_times_M(back, linkin, linkmu, t->sub);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, t->sub);
  QDP_M_eq_M_times_M(linkinnu, link, innu, t->sub);
  QDP_discard_M(innu);
  QDP_M_eq_M_times_Ma(outf, linkinnu, linkmu, t->sub);
  QDP_discard_M(linkmu);
  QDP_M_eq_M(outb, backnu, t->sub);
  QDP_discard_M(backnu);
#define STAPLEFB_FLOPS (4*EQMTM)

#undef link
#undef linkmu
#undef in
#undef innu
#undef linkin
#undef back
#undef backnu
#undef linkinnu
}
#endif

static void
staples(QDP_ColorMatrix *out, QDP_ColorMatrix *top0, QDP_ColorMatrix *bot,
	QDP_ColorMatrix *left, QDP_ColorMatrix *right0,
	int mu, int nu, tmpstruct *t)
{
#define right     t->ftmp0[0]
#define rightmu   t->ftmp[0][mu]
#define top       t->ftmp0[1]
#define topnu     t->ftmp[1][nu]
#define leftbot   t->mtmp[0]
#define back      t->btmp0[0]
#define backnu    t->btmp[0][nu]
#define lefttopnu t->mtmp[1]

  QDP_M_eq_M(right, right0, t->sub);
  QDP_M_eq_sM(rightmu, right, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(top, top0, t->sub);
  QDP_M_eq_sM(topnu, top, QDP_neighbor[nu], QDP_forward, t->sub);
  QDP_M_eq_Ma_times_M(leftbot, left, bot, t->sub);
  QDP_M_eq_M_times_M(back, leftbot, rightmu, t->sub);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, t->sub);
  QDP_M_eq_M_times_M(lefttopnu, left, topnu, t->sub);
  QDP_discard_M(topnu);
  QDP_M_peq_M_times_Ma(out, lefttopnu, rightmu, t->sub);
  QDP_discard_M(rightmu);
  QDP_M_peq_M(out, backnu, t->sub);
  QDP_discard_M(backnu);
#define STAPLES_FLOPS (3*EQMTM+PEQMTM+PEQM)

#undef right
#undef rightmu
#undef top
#undef topnu
#undef leftbot
#undef back
#undef backnu
#undef lefttopnu
}

static void
staple2fb(QDP_ColorMatrix *outmuf, QDP_ColorMatrix *outmub,
	  QDP_ColorMatrix *outnuf, QDP_ColorMatrix *outnub,
	  QDP_ColorMatrix *Umu0, QDP_ColorMatrix *Unu0,
	  int mu, int nu, tmpstruct *t)
{
#define Unu        t->ftmp0[0]
#define Unufmu     t->ftmp[0][mu]
#define Umu        t->ftmp0[1]
#define Umufnu     t->ftmp[1][nu]
#define UmuUnufmu  t->mtmp[0]
#define backmu     t->btmp0[0]
#define backmubnu  t->btmp[0][nu]
#define UnuUmufnu  t->mtmp[1]
#define backnu     t->btmp0[1]
#define backnubmu  t->btmp[1][mu]

  // fmu += Unu * Umufnu * Unufmu+
  // fmu += Unubnu+ * Umubnu * Unufmubnu
  // fnu += Umu * Unufmu * Umufnu+
  // fnu += Umubmu+ * Unubmu * Umufnubmu
  QDP_M_eq_M(Unu, Unu0, t->sub);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(Umu, Umu0, t->sub);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t->sub);
  QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, t->sub);
  QDP_M_eq_Ma_times_M(backmu, Unu, UmuUnufmu, t->sub);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t->sub);
  QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, t->sub);
  QDP_M_eq_Ma_times_M(backnu, Umu, UnuUmufnu, t->sub);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t->sub);
  QDP_M_eq_M_times_Ma(outmuf, UnuUmufnu, Unufmu, t->sub);
  QDP_M_eq_M_times_Ma(outnuf, UmuUnufmu, Umufnu, t->sub);
  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_M_eq_M(outmub, backmubnu, t->sub);
  QDP_discard_M(backmubnu);
  QDP_M_eq_M(outnub, backnubmu, t->sub);
  QDP_discard_M(backnubmu);
#define STAPLE2FB_FLOPS (6*EQMTM)

#undef Unu
#undef Unufmu
#undef Umu
#undef Umufnu
#undef UmuUnufmu
#undef backmu
#undef backmubnu
#undef UnuUmufnu
#undef backnu
#undef backnubmu
}

static void
stapler(QDP_ColorMatrix *outmu, QDP_ColorMatrix *outnu,
	QDP_ColorMatrix *Umu0, QDP_ColorMatrix *Unu0,
	QDP_ColorMatrix *Fmu0, QDP_ColorMatrix *Fnu0,
	QDP_ColorMatrix *Bmu, QDP_ColorMatrix *Bnu,
	int mu, int nu, tmpstruct *t)
{
#define Unu        t->ftmp0[0]
#define Unufmu     t->ftmp[0][mu]
#define Umu        t->ftmp0[1]
#define Umufnu     t->ftmp[1][nu]
#define Fnu        t->ftmp0[2]
#define Fnufmu     t->ftmp[2][mu]
#define Fmu        t->ftmp0[3]
#define Fmufnu     t->ftmp[3][nu]
#define UmuUnufmu  t->mtmp[0]
#define BmuUnufmu  t->mtmp[1]
  //#define UmuFnufmu  mtmp[4]
#define backmu     t->btmp0[0]
#define backmubnu  t->btmp[0][nu]
#define UnuUmufnu  t->mtmp[2]
#define BnuUmufnu  t->mtmp[3]
  //#define UnuFmufnu  mtmp[5]
#define backnu     t->btmp0[1]
#define backnubmu  t->btmp[1][mu]

  // fmu += Unu * Umufnu * Fnufmu+
  // fmu += Unu * Fmufnu * Unufmu+
  // fmu += Bnu * Umufnu * Unufmu+
  // fmu += Unubnu+ * Umubnu * Fnufmubnu
  // fmu += Unubnu+ * Bmubnu * Unufmubnu
  // fmu += Bnubnu+ * Umubnu * Unufmubnu
  // fnu += Umu * Unufmu * Fmufnu+
  // fnu += Umu * Fnufmu * Umufnu+
  // fnu += Bmu * Unufmu * Umufnu+
  // fnu += Umubmu+ * Unubmu * Fmufnubmu
  // fnu += Umubmu+ * Bnubmu * Umufnubmu
  // fnu += Bmubmu+ * Unubmu * Umufnubmu
  QDP_M_eq_M(Unu, Unu0, t->sub);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(Umu, Umu0, t->sub);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t->sub);
  QDP_M_eq_M(Fnu, Fnu0, t->sub);
  QDP_M_eq_sM(Fnufmu, Fnu, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(Fmu, Fmu0, t->sub);
  QDP_M_eq_sM(Fmufnu, Fmu, QDP_neighbor[nu], QDP_forward, t->sub);

  QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, t->sub);
  QDP_M_eq_Ma_times_M(backmu, Bnu, UmuUnufmu, t->sub);
  QDP_M_eq_M_times_M(BmuUnufmu, Bmu, Unufmu, t->sub);
  QDP_M_peq_M_times_M(BmuUnufmu, Umu, Fnufmu, t->sub);
  QDP_M_peq_Ma_times_M(backmu, Unu, BmuUnufmu, t->sub);
  //QDP_M_eq_M_times_M(UmuFnufmu, Umu, Fnufmu, t->sub);
  //QDP_M_peq_Ma_times_M(backmu, Unu, UmuFnufmu, t->sub);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t->sub);

  QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, t->sub);
  QDP_M_eq_Ma_times_M(backnu, Bmu, UnuUmufnu, t->sub);
  QDP_M_eq_M_times_M(BnuUmufnu, Bnu, Umufnu, t->sub);
  QDP_M_peq_M_times_M(BnuUmufnu, Unu, Fmufnu, t->sub);
  QDP_M_peq_Ma_times_M(backnu, Umu, BnuUmufnu, t->sub);
  //QDP_M_eq_M_times_M(UnuFmufnu, Unu, Fmufnu, t->sub);
  //QDP_M_peq_Ma_times_M(backnu, Umu, UnuFmufnu, t->sub);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t->sub);

  QDP_M_peq_M_times_Ma(outmu, UnuUmufnu, Fnufmu, t->sub);
  //QDP_M_peq_M_times_Ma(outmu, UnuFmufnu, Unufmu, t->sub);
  //QDP_M_peq_M(BnuUmufnu, UnuFmufnu, t->sub);
  QDP_M_peq_M_times_Ma(outmu, BnuUmufnu, Unufmu, t->sub);
  QDP_M_peq_M(outmu, backmubnu, t->sub);

  QDP_M_peq_M_times_Ma(outnu, UmuUnufmu, Fmufnu, t->sub);
  //QDP_M_peq_M_times_Ma(outnu, UmuFnufmu, Umufnu, t->sub);
  //QDP_M_peq_M(BmuUnufmu, UmuFnufmu, t->sub);
  QDP_M_peq_M_times_Ma(outnu, BmuUnufmu, Umufnu, t->sub);
  QDP_M_peq_M(outnu, backnubmu, t->sub);

  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_discard_M(Fmufnu);
  QDP_discard_M(Fnufmu);
  QDP_discard_M(backmubnu);
  QDP_discard_M(backnubmu);
  //#define STAPLER_FLOPS (8*198+10*216+2*18)
  //#define STAPLER_FLOPS (8*198+8*216+4*18)
#define STAPLER_FLOPS (6*EQMTM+8*PEQMTM+2*PEQM)

#undef Unu
#undef Unufmu
#undef Umu
#undef Umufnu
#undef Fnu
#undef Fnufmu
#undef Fmu
#undef Fmufnu
#undef UmuUnufmu
#undef BmuUnufmu
    //#undef UmuFnufmu
#undef backmu
#undef backmubnu
#undef UnuUmufnu
#undef BnuUmufnu
    //#undef UnuFmufnu
#undef backnu
#undef backnubmu
}

static void
staplep(QDP_ColorMatrix *outmu, QDP_ColorMatrix *outnu,
	QDP_ColorMatrix *Umu0, QDP_ColorMatrix *Unu0,
	QDP_ColorMatrix *Xmu0[], QDP_ColorMatrix *Xnu0[],
	int nx, int mu, int nu, tmpstruct *t)
{
#define XXnumu     t->mtmp[0]
#define XXmunu     t->mtmp[1]
#define XXdnumu    t->mtmp[2]
#define XXdmunu    t->mtmp[3]
#define XdXnumu    t->mtmp[4]
#define XdXmunu    t->mtmp[5]
#define Xnu        t->ftmp0[0]
#define Xnufmu     t->ftmp[0][mu]
#define Xmu        t->ftmp0[1]
#define Xmufnu     t->ftmp[1][nu]
#define Unu        t->ftmp0[2]
#define Unufmu     t->ftmp[2][mu]
#define Umu        t->ftmp0[3]
#define Umufnu     t->ftmp[3][nu]
#define backmu     t->btmp0[0]
#define backmubnu  t->btmp[0][nu]
#define backnu     t->btmp0[1]
#define backnubmu  t->btmp[1][mu]

  // fmu += Unu * Xmufnu * Xnufmu+
  // fmu += Unubnu+ * Xmubnu * Xnufmubnu
  // fmu += Xnu * Xmufnu * Unufmu+
  // fmu += Xnubnu+ * Xmubnu * Unufmubnu
  // fnu += Umu * Xnufmu * Xmufnu+
  // fnu += Umubmu+ * Xnubmu * Xmufnubmu
  // fnu += Xmu * Xnufmu * Umufnu+
  // fnu += Xmubmu+ * Xnubmu * Umufnubmu
  QDP_M_eq_zero(XXnumu, t->sub);
  QDP_M_eq_zero(XXmunu, t->sub);
  QDP_M_eq_zero(XXdnumu, t->sub);
  QDP_M_eq_zero(XXdmunu, t->sub);
  QDP_M_eq_zero(XdXnumu, t->sub);
  QDP_M_eq_zero(XdXmunu, t->sub);
  for(int i=0; i<nx; i++) {
    QDP_M_eq_M(Xnu, Xnu0[i], t->sub);
    QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, t->sub);
    QDP_M_eq_M(Xmu, Xmu0[i], t->sub);
    QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, t->sub);
    QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, t->sub);
    QDP_M_peq_Ma_times_M(XdXmunu, Xmu, Xnu, t->sub);
    QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, t->sub);
    QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, t->sub);
    QDP_M_peq_M_times_Ma(XXdnumu, Xnufmu, Xmufnu, t->sub);
    QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t->sub);
    QDP_discard_M(Xmufnu);
    QDP_discard_M(Xnufmu);
  }
  QDP_M_eq_M(Unu, Unu0, t->sub);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t->sub);
  QDP_M_eq_M(Umu, Umu0, t->sub);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t->sub);

  QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, t->sub);
  QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, t->sub);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t->sub);
  QDP_M_eq_Ma_times_M(backnu, Umu, XXnumu, t->sub);
  QDP_M_peq_M_times_M(backnu, XdXmunu, Umufnu, t->sub);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t->sub);

  QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, t->sub);
  QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, t->sub);
  QDP_M_peq_M(outmu, backmubnu, t->sub);
  QDP_M_peq_M_times_M(outnu, Umu, XXdnumu, t->sub);
  QDP_M_peq_M_times_Ma(outnu, XXmunu, Umufnu, t->sub);
  QDP_M_peq_M(outnu, backnubmu, t->sub);

  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_discard_M(backmubnu);
  QDP_discard_M(backnubmu);
#define STAPLEP_FLOPS(n) (2*EQMTM+(6*n+6)*PEQMTM+2*PEQM)

#undef XXnumu
#undef XXmunu
#undef XXdnumu
#undef XXdmunu
#undef XdXnumu
#undef XdXmunu
#undef Xnu
#undef Xnufmu
#undef Xmu
#undef Xmufnu
#undef Unu
#undef Unufmu
#undef Umu
#undef Umufnu
#undef backmu
#undef backmubnu
#undef backnu
#undef backnubmu
}

static void
staplep2(QDP_ColorMatrix *outmu, QDP_ColorMatrix *Unu0,
	 QDP_ColorMatrix *Xmu0[], QDP_ColorMatrix *Xnu0[],
	 int nx, int mu, int nu, tmpstruct *t)
{
#define XXnumu     t->mtmp[0]
#define XXmunu     t->mtmp[1]
#define XXdmunu    t->mtmp[2]
#define XdXnumu    t->mtmp[3]
#define Xnu        t->ftmp0[0]
#define Xnufmu     t->ftmp[0][mu]
#define Xmu        t->ftmp0[1]
#define Xmufnu     t->ftmp[1][nu]
#define Unu        t->ftmp0[2]
#define Unufmu     t->ftmp[2][mu]
#define backmu     t->btmp0[0]
#define backmubnu  t->btmp[0][nu]

  // fmu += Unu * Xmufnu * Xnufmu+
  // fmu += Unubnu+ * Xmubnu * Xnufmubnu
  // fmu += Xnu * Xmufnu * Unufmu+
  // fmu += Xnubnu+ * Xmubnu * Unufmubnu
  //QDP_M_eq_zero(XXdmunu, t->sub);
  //QDP_M_eq_zero(XXmunu, t->sub);
  //QDP_M_eq_zero(XXnumu, t->sub);
  //QDP_M_eq_zero(XdXnumu, t->sub);
  for(int i=0; i<nx; i++) {
    QDP_M_eq_M(Xnu, Xnu0[i], t->sub);
    QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, t->sub);
    QDP_M_eq_M(Xmu, Xmu0[i], t->sub);
    QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, t->sub);
    if(i==0) {
      QDP_M_eq_Ma_times_M(XdXnumu, Xnu, Xmu, t->sub);
      QDP_M_eq_M_times_M(XXnumu, Xnu, Xmufnu, t->sub);
      QDP_M_eq_M_times_M(XXmunu, Xmu, Xnufmu, t->sub);
      QDP_M_eq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t->sub);
    } else {
      QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, t->sub);
      QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, t->sub);
      QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, t->sub);
      QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t->sub);
    }
    QDP_discard_M(Xmufnu);
    QDP_discard_M(Xnufmu);
  }
  QDP_M_eq_M(Unu, Unu0, t->sub);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t->sub);

  QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, t->sub);
  QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, t->sub);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t->sub);

  QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, t->sub);
  QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, t->sub);
  QDP_M_peq_M(outmu, backmubnu, t->sub);

  QDP_discard_M(Unufmu);
  QDP_discard_M(backmubnu);
  //#define STAPLEP2_FLOPS(n) (198+(4*n+3)*216+18)
#define STAPLEP2_FLOPS(n) (5*EQMTM+(4*(n-1)+3)*PEQMTM+PEQM)

#undef XXnumu
#undef XXmunu
#undef XXdmunu
#undef XdXnumu
#undef Xnu
#undef Xnufmu
#undef Xmu
#undef Xmufnu
#undef Unu
#undef Unufmu
#undef backmu
#undef backmubnu
}

#define G(x) gauge->links[x]

// TODO: implement chain rule
void
QOP_symanzik_1loop_gauge_deriv_qdp(QOP_info_t *info, QOP_GaugeField *gauge, 
				   QDP_ColorMatrix *deriv[],
				   QOP_gauge_coeffs_t *coeffs,
				   REAL eps, int doLastScale)
{
#define NC QDP_get_nc(gauge->links[0])
  double dtime = QOP_time();
  double nflops = 0;
  QDP_Lattice *lat = QDP_get_lattice_M(gauge->links[0]);
  int nd = QDP_ndim_L(lat);
  QDP_Subset sub = QDP_all_L(lat);

  QLA_Real fac = -eps/QLA_Nc;
  QLA_Real plaq = fac*coeffs->plaquette;
  QLA_Real rect = fac*coeffs->rectangle;
  QLA_Real pgm  = fac*coeffs->parallelogram;
  QLA_Real adpl = 2*(fac/QLA_Nc)*coeffs->adjoint_plaquette;

  QDP_ColorMatrix *mtmp[NMTMP], *ftmp0[NFTMP], *ftmp[NFTMP][nd],
    *btmp0[NBTMP], *btmp[NBTMP][nd];
  tmpstruct t;
  t.nc = QOP_Nc;
  t.nd = nd;
  t.sub = sub;
  t.lat = lat;
  t.mtmp = mtmp;
  t.ftmp0 = ftmp0;
  for(int i=0; i<NFTMP; i++) t.ftmp[i] = ftmp[i];
  t.btmp0 = btmp0;
  for(int i=0; i<NBTMP; i++) t.btmp[i] = btmp[i];
  set_temps(&t);
  TRACE;

  QDP_ColorMatrix *tforce[nd], *stplf[nd][nd], *stplb[nd][nd], *tmat[nd];
  for(int mu=0; mu<nd; mu++) {
    tforce[mu] = QDP_create_M_L(lat);
    QDP_M_eq_zero(tforce[mu], sub);
  }
  if(rect||pgm) {
    for(int mu=0; mu<nd; mu++) {
      tmat[mu] = QDP_create_M_L(lat);
    }
  }
  TRACE;
  for(int mu=1; mu<nd; mu++) {
    for(int nu=0; nu<mu; nu++) {
      stplf[mu][nu] = QDP_create_M_L(lat);
      stplb[mu][nu] = QDP_create_M_L(lat);
      stplf[nu][mu] = QDP_create_M_L(lat);
      stplb[nu][mu] = QDP_create_M_L(lat);
      staple2fb(stplf[mu][nu], stplb[mu][nu], stplf[nu][mu], stplb[nu][mu],
		G(mu), G(nu), mu, nu, &t);
      nflops += STAPLE2FB_FLOPS;
      if(adpl) {
	QDP_Complex *tc = QDP_create_C_L(lat);
	QLA_Complex z;
	QLA_c_eq_r(z, plaq/adpl);

	QDP_C_eq_c(tc, &z, sub);
	QDP_C_peq_M_dot_M(tc, stplf[mu][nu], G(mu), sub);
	QDP_C_eq_r_times_C(tc, &adpl, tc, sub);
	QDP_M_peq_C_times_M(tforce[mu], tc, stplf[mu][nu], sub);

	QDP_C_eq_c(tc, &z, sub);
	QDP_C_peq_M_dot_M(tc, stplb[mu][nu], G(mu), sub);
	QDP_C_eq_r_times_C(tc, &adpl, tc, sub);
	QDP_M_peq_C_times_M(tforce[mu], tc, stplb[mu][nu], sub);

	QDP_C_eq_c(tc, &z, sub);
	QDP_C_peq_M_dot_M(tc, stplf[nu][mu], G(nu), sub);
	QDP_C_eq_r_times_C(tc, &adpl, tc, sub);
	QDP_M_peq_C_times_M(tforce[nu], tc, stplf[nu][mu], sub);

	QDP_C_eq_c(tc, &z, sub);
	QDP_C_peq_M_dot_M(tc, stplb[nu][mu], G(nu), sub);
	QDP_C_eq_r_times_C(tc, &adpl, tc, sub);
	QDP_M_peq_C_times_M(tforce[nu], tc, stplb[nu][mu], sub);

	QDP_destroy_C(tc);
	nflops += 4*(16*NC*NC+2);
      } else
      if(plaq) {
	QDP_M_peq_r_times_M(tforce[mu], &plaq, stplf[mu][nu], sub);
	QDP_M_peq_r_times_M(tforce[mu], &plaq, stplb[mu][nu], sub);
	QDP_M_peq_r_times_M(tforce[nu], &plaq, stplf[nu][mu], sub);
	QDP_M_peq_r_times_M(tforce[nu], &plaq, stplb[nu][mu], sub);
	nflops += 4*4*NC*NC;
      }
    }
  }
  TRACE;

  if(rect) {
    for(int mu=0; mu<nd; mu++) {
      QDP_M_eq_zero(tmat[mu], sub);
    }
    for(int mu=1; mu<nd; mu++) {
      for(int nu=0; nu<mu; nu++) {
	//if(nu==mu) continue;
	//staples(tmat[mu], G(mu), G(mu), stplb[nu][mu], G(nu), mu, nu);
	//staples(tmat[mu], G(mu), G(mu), G(nu), stplf[nu][mu], mu, nu);
	//staples(tmat[mu], stplf[mu][nu], stplb[mu][nu], G(nu), G(nu), mu,nu);
	stapler(tmat[mu], tmat[nu], G(mu), G(nu), stplf[mu][nu], stplf[nu][mu],
		stplb[mu][nu], stplb[nu][mu], mu, nu, &t);
	nflops += STAPLER_FLOPS;
      }
    }
    for(int mu=0; mu<nd; mu++) {
      QDP_M_peq_r_times_M(tforce[mu], &rect, tmat[mu], sub);
    }
    nflops += nd*4*NC*NC;
  }
  TRACE;

  if(pgm) {
    if(nd!=4) {

      for(int mu=0; mu<nd; mu++) {
	QDP_M_eq_zero(tmat[0], sub);
	for(int nu=0; nu<nd; nu++) {
	  if(nu==mu) continue;
	  for(int rho=nu+1; rho<nd; rho++) {
	    if(rho==mu||rho==nu) continue;
	    staples(tmat[0], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu, &t);
	    staples(tmat[0], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu, &t);
	    staples(tmat[0], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu, &t);
	    staples(tmat[0], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu, &t);
	    nflops += 4*STAPLES_FLOPS;
	  }
	}
	QDP_M_peq_r_times_M(tforce[mu], &pgm, tmat[0], sub);
	nflops += 4*NC*NC;
      }

    } else { // nd == 4

      for(int mu=0; mu<4; mu++) {
	QDP_M_eq_zero(tmat[mu], sub);
      }
      {
	int mu, nu, rho;
	QDP_ColorMatrix *Xmu[4], *Xnu[4];
#if 1
	mu=0; nu=1;
	Xmu[0] = stplf[mu][2];
	Xmu[1] = stplf[mu][3];
	Xmu[2] = stplb[mu][2];
	Xmu[3] = stplb[mu][3];
	Xnu[0] = stplf[nu][2];
	Xnu[1] = stplf[nu][3];
	Xnu[2] = stplb[nu][2];
	Xnu[3] = stplb[nu][3];
	staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu, &t);
	nflops += STAPLEP_FLOPS(4);
#else
	mu=0; nu=1; rho=2;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu, &t);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu, &t);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu, &t);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu, &t);
	mu=0; nu=1; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=1; nu=0; rho=2;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=1; nu=0; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	nflops += 16*STAPLES_FLOPS;
#endif
#if 1
	mu=2; nu=3;
	Xmu[0] = stplf[mu][0];
	Xmu[1] = stplf[mu][1];
	Xmu[2] = stplb[mu][0];
	Xmu[3] = stplb[mu][1];
	Xnu[0] = stplf[nu][0];
	Xnu[1] = stplf[nu][1];
	Xnu[2] = stplb[nu][0];
	Xnu[3] = stplb[nu][1];
	staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu, &t);
	nflops += STAPLEP_FLOPS(4);
#else
	mu=2; nu=0; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=2; nu=1; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=3; nu=0; rho=2;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=3; nu=1; rho=2;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	nflops += 16*STAPLES_FLOPS;
#endif
#if 1
	mu=0; nu=2; rho=3;
	Xmu[0] = stplf[mu][rho];
	Xmu[1] = stplb[mu][rho];
	Xnu[0] = stplf[nu][rho];
	Xnu[1] = stplb[nu][rho];
	staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, &t);
	mu=1; nu=2; rho=3;
	Xmu[0] = stplf[mu][rho];
	Xmu[1] = stplb[mu][rho];
	Xnu[0] = stplf[nu][rho];
	Xnu[1] = stplb[nu][rho];
	staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, &t);
	mu=2; nu=0; rho=1;
	Xmu[0] = stplf[mu][rho];
	Xmu[1] = stplb[mu][rho];
	Xnu[0] = stplf[nu][rho];
	Xnu[1] = stplb[nu][rho];
	staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, &t);
	mu=3; nu=0; rho=1;
	Xmu[0] = stplf[mu][rho];
	Xmu[1] = stplb[mu][rho];
	Xnu[0] = stplf[nu][rho];
	Xnu[1] = stplb[nu][rho];
	staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, &t);
	nflops += 4*STAPLEP2_FLOPS(2);
#else
	mu=0; nu=2; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=1; nu=2; rho=3;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=2; nu=0; rho=1;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	mu=3; nu=0; rho=1;
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	nflops += 16*STAPLES_FLOPS;
#endif
      }
      for(int mu=0; mu<4; mu++) {
	QDP_M_peq_r_times_M(tforce[mu], &pgm, tmat[mu], sub);
      }
      nflops += 4*4*NC*NC;

    } // nd!=4
  } // pgm
  TRACE;

  for(int mu=0; mu<nd; mu++) {
    QDP_M_peq_M(deriv[mu], tforce[mu], sub);
  }
  nflops += nd*PEQM;

  if(rect||pgm) {
    for(int mu=0; mu<nd; mu++) {
      QDP_destroy_M(tmat[mu]);
    }
  }
  for(int mu=0; mu<nd; mu++) {
    QDP_destroy_M(tforce[mu]);
    for(int nu=0; nu<nd; nu++) {
      if(nu==mu) continue;
      QDP_destroy_M(stplf[mu][nu]);
      QDP_destroy_M(stplb[mu][nu]);
    }
  }
  free_temps(&t);

  //double nflop = 96720 - 4*(24+18);
  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
#undef NC
} 

// TODO: use top gauge field in force
void
QOP_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, QOP_GaugeField *gauge,
				   QDP_ColorMatrix *force[],
				   QOP_gauge_coeffs_t *coeffs,
				   QLA_Real eps)
{
#define NC QDP_get_nc(gauge->links[0])
  double dtime = QOP_time();
  QDP_Lattice *lat = QDP_get_lattice_M(gauge->links[0]);
  int nd = QDP_ndim_L(lat);
  QDP_Subset sub = QDP_all_L(lat);

  QDP_ColorMatrix *deriv[nd];
  for(int mu=0; mu<nd; mu++) {
    deriv[mu] = QDP_create_M_L(lat);
    QDP_M_eq_zero(deriv[mu], sub);
  }
  QOP_symanzik_1loop_gauge_deriv_qdp(info, gauge, deriv, coeffs, eps, 0);

  QDP_ColorMatrix *mtmp = QDP_create_M_L(lat);
#ifdef CHKSUM
  QLA_ColorMatrix(qcm);
  QLA_Complex det, chk;
  QLA_c_eq_r(chk, 0);
#endif
  QLA_Real trace=0;
  for(int mu=0; mu<nd; mu++) {
    QDP_M_eq_M_times_Ma(mtmp, G(mu), deriv[mu], sub);
    if(QOP_common.verbosity==QOP_VERB_DEBUG) {
      QLA_ColorMatrix(tm);
      QLA_Real tr;
      QDP_m_eq_sum_M(&tm, mtmp, sub);
      QLA_R_eq_re_trace_M(&tr, &tm);
      trace += tr;
    }
    QDP_M_eq_antiherm_M(deriv[mu], mtmp, sub);
    QDP_M_peq_M(force[mu], deriv[mu], sub);
#ifdef CHKSUM
    QDP_m_eq_sum_M(&qcm, force->force[mu], sub);
    QLA_C_eq_det_M(&det, &qcm);
    QLA_c_peq_c(chk, det);
#endif
  }
#ifdef CHKSUM
  QOP_printf0("chksum: %g %g\n", QLA_real(chk), QLA_imag(chk));
#endif
  if(QOP_common.verbosity==QOP_VERB_DEBUG) {
    QOP_printf0("re trace: %g\n", trace);
  }
  info->final_flop += nd*(EQMTM+2*NC*NC+PEQM)*QDP_subset_len(sub);

  QDP_destroy_M(mtmp);
  for(int mu=0; mu<nd; mu++) {
    QDP_destroy_M(deriv[mu]);
  }

  info->final_sec = QOP_time() - dtime;
#undef NC
}

void
QOP_symanzik_1loop_gauge_force(QOP_info_t *info, QOP_GaugeField *gauge, 
			       QOP_Force *force, QOP_gauge_coeffs_t *coeffs,
			       REAL eps)
{
  QOP_symanzik_1loop_gauge_force_qdp(info, gauge, force->force, coeffs, eps);
}
