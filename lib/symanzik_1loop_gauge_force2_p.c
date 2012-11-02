/* gauge force for Symanzik improved 1x1 + 1x2 + 1x1x1 */

#include <qop_internal.h>

#define NMTMP 6
#define NFTMP 4
#define NBTMP 2
static QDP_ColorMatrix *mtmp[NMTMP], *ftmp0[NFTMP], *ftmp[NFTMP][4],
  *btmp0[NBTMP], *btmp[NBTMP][4];
static int setcount=0;
#define set_temps() if(!setcount++) set_temps0()
#define free_temps() if(!--setcount) free_temps0()

static void
set_temps0(void)
{
  for(int i=0; i<NMTMP; i++) mtmp[i] = QDP_create_M();
  for(int i=0; i<NFTMP; i++) {
    ftmp0[i] = QDP_create_M();
    for(int j=0; j<4; j++) ftmp[i][j] = QDP_create_M();
  }
  for(int i=0; i<NBTMP; i++) {
    btmp0[i] = QDP_create_M();
    for(int j=0; j<4; j++) btmp[i][j] = QDP_create_M();
  }
}

static void
free_temps0(void)
{
  for(int i=0; i<NMTMP; i++) QDP_destroy_M(mtmp[i]);
  for(int i=0; i<NFTMP; i++) {
    QDP_destroy_M(ftmp0[i]);
    for(int j=0; j<4; j++) QDP_destroy_M(ftmp[i][j]);
  }
  for(int i=0; i<NBTMP; i++) {
    QDP_destroy_M(btmp0[i]);
    for(int j=0; j<4; j++) QDP_destroy_M(btmp[i][j]);
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
  QDP_M_eq_M(link, link0, QDP_all);
  QDP_M_eq_sM(linkmu, link, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(in, in0, QDP_all);
  QDP_M_eq_sM(innu, in, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(linkin, link, in, QDP_all);
  QDP_M_eq_M_times_M(back, linkin, linkmu, QDP_all);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(linkinnu, link, innu, QDP_all);
  QDP_discard_M(innu);
  QDP_M_eq_M_times_Ma(outf, linkinnu, linkmu, QDP_all);
  QDP_discard_M(linkmu);
  QDP_M_eq_M(outb, backnu, QDP_all);
  QDP_discard_M(backnu);
#define STAPLEFB_FLOPS (4*198)

#undef link
#undef linkmu
#undef in
#undef innu
#undef linkin
#undef back
#undef backnu
#undef linkinnu
}

static void
staples(QDP_ColorMatrix *out, QDP_ColorMatrix *top0, QDP_ColorMatrix *bot,
	QDP_ColorMatrix *left, QDP_ColorMatrix *right0, int mu, int nu)
{
#define right     ftmp0[0]
#define rightmu   ftmp[0][mu]
#define top       ftmp0[1]
#define topnu     ftmp[1][nu]
#define leftbot   mtmp[0]
#define back      btmp0[0]
#define backnu    btmp[0][nu]
#define lefttopnu mtmp[1]

  QDP_M_eq_M(right, right0, QDP_all);
  QDP_M_eq_sM(rightmu, right, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(top, top0, QDP_all);
  QDP_M_eq_sM(topnu, top, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(leftbot, left, bot, QDP_all);
  QDP_M_eq_M_times_M(back, leftbot, rightmu, QDP_all);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(lefttopnu, left, topnu, QDP_all);
  QDP_discard_M(topnu);
  QDP_M_peq_M_times_Ma(out, lefttopnu, rightmu, QDP_all);
  QDP_discard_M(rightmu);
  QDP_M_peq_M(out, backnu, QDP_all);
  QDP_discard_M(backnu);
#define STAPLES_FLOPS (3*198+216+18)

#undef right
#undef rightmu
#undef top
#undef topnu
#undef leftbot
#undef back
#undef backnu
#undef lefttopnu
}
#endif

static void
staple2fb(QDP_ColorMatrix *outmuf, QDP_ColorMatrix *outmub,
	  QDP_ColorMatrix *outnuf, QDP_ColorMatrix *outnub,
	  QDP_ColorMatrix *Umu0, QDP_ColorMatrix *Unu0, int mu, int nu)
{
#define Unu        ftmp0[0]
#define Unufmu     ftmp[0][mu]
#define Umu        ftmp0[1]
#define Umufnu     ftmp[1][nu]
#define UmuUnufmu  mtmp[0]
#define backmu     btmp0[0]
#define backmubnu  btmp[0][nu]
#define UnuUmufnu  mtmp[1]
#define backnu     btmp0[1]
#define backnubmu  btmp[1][mu]

  // fmu += Unu * Umufnu * Unufmu+
  // fmu += Unubnu+ * Umubnu * Unufmubnu
  // fnu += Umu * Unufmu * Umufnu+
  // fnu += Umubmu+ * Unubmu * Umufnubmu
  QDP_M_eq_M(Unu, Unu0, QDP_all);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(Umu, Umu0, QDP_all);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, QDP_all);
  QDP_M_eq_Ma_times_M(backmu, Unu, UmuUnufmu, QDP_all);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, QDP_all);
  QDP_M_eq_Ma_times_M(backnu, Umu, UnuUmufnu, QDP_all);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_Ma(outmuf, UnuUmufnu, Unufmu, QDP_all);
  QDP_M_eq_M_times_Ma(outnuf, UmuUnufmu, Umufnu, QDP_all);
  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_M_eq_M(outmub, backmubnu, QDP_all);
  QDP_discard_M(backmubnu);
  QDP_M_eq_M(outnub, backnubmu, QDP_all);
  QDP_discard_M(backnubmu);
#define STAPLE2FB_FLOPS (6*198)

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
	QDP_ColorMatrix *Bmu, QDP_ColorMatrix *Bnu, int mu, int nu)
{
#define Unu        ftmp0[0]
#define Unufmu     ftmp[0][mu]
#define Umu        ftmp0[1]
#define Umufnu     ftmp[1][nu]
#define Fnu        ftmp0[2]
#define Fnufmu     ftmp[2][mu]
#define Fmu        ftmp0[3]
#define Fmufnu     ftmp[3][nu]
#define UmuUnufmu  mtmp[0]
#define BmuUnufmu  mtmp[1]
  //#define UmuFnufmu  mtmp[4]
#define backmu     btmp0[0]
#define backmubnu  btmp[0][nu]
#define UnuUmufnu  mtmp[2]
#define BnuUmufnu  mtmp[3]
  //#define UnuFmufnu  mtmp[5]
#define backnu     btmp0[1]
#define backnubmu  btmp[1][mu]

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
  QDP_M_eq_M(Unu, Unu0, QDP_all);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(Umu, Umu0, QDP_all);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_M(Fnu, Fnu0, QDP_all);
  QDP_M_eq_sM(Fnufmu, Fnu, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(Fmu, Fmu0, QDP_all);
  QDP_M_eq_sM(Fmufnu, Fmu, QDP_neighbor[nu], QDP_forward, QDP_all);

  QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, QDP_all);
  QDP_M_eq_Ma_times_M(backmu, Bnu, UmuUnufmu, QDP_all);
  QDP_M_eq_M_times_M(BmuUnufmu, Bmu, Unufmu, QDP_all);
  QDP_M_peq_M_times_M(BmuUnufmu, Umu, Fnufmu, QDP_all);
  QDP_M_peq_Ma_times_M(backmu, Unu, BmuUnufmu, QDP_all);
  //QDP_M_eq_M_times_M(UmuFnufmu, Umu, Fnufmu, QDP_all);
  //QDP_M_peq_Ma_times_M(backmu, Unu, UmuFnufmu, QDP_all);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, QDP_all);

  QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, QDP_all);
  QDP_M_eq_Ma_times_M(backnu, Bmu, UnuUmufnu, QDP_all);
  QDP_M_eq_M_times_M(BnuUmufnu, Bnu, Umufnu, QDP_all);
  QDP_M_peq_M_times_M(BnuUmufnu, Unu, Fmufnu, QDP_all);
  QDP_M_peq_Ma_times_M(backnu, Umu, BnuUmufnu, QDP_all);
  //QDP_M_eq_M_times_M(UnuFmufnu, Unu, Fmufnu, QDP_all);
  //QDP_M_peq_Ma_times_M(backnu, Umu, UnuFmufnu, QDP_all);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, QDP_all);

  QDP_M_peq_M_times_Ma(outmu, UnuUmufnu, Fnufmu, QDP_all);
  //QDP_M_peq_M_times_Ma(outmu, UnuFmufnu, Unufmu, QDP_all);
  //QDP_M_peq_M(BnuUmufnu, UnuFmufnu, QDP_all);
  QDP_M_peq_M_times_Ma(outmu, BnuUmufnu, Unufmu, QDP_all);
  QDP_M_peq_M(outmu, backmubnu, QDP_all);

  QDP_M_peq_M_times_Ma(outnu, UmuUnufmu, Fmufnu, QDP_all);
  //QDP_M_peq_M_times_Ma(outnu, UmuFnufmu, Umufnu, QDP_all);
  //QDP_M_peq_M(BmuUnufmu, UmuFnufmu, QDP_all);
  QDP_M_peq_M_times_Ma(outnu, BmuUnufmu, Umufnu, QDP_all);
  QDP_M_peq_M(outnu, backnubmu, QDP_all);

  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_discard_M(Fmufnu);
  QDP_discard_M(Fnufmu);
  QDP_discard_M(backmubnu);
  QDP_discard_M(backnubmu);
  //#define STAPLER_FLOPS (8*198+10*216+2*18)
  //#define STAPLER_FLOPS (8*198+8*216+4*18)
#define STAPLER_FLOPS (6*198+8*216+2*18)

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
	int nx, int mu, int nu)
{
#define XXnumu     mtmp[0]
#define XXmunu     mtmp[1]
#define XXdnumu    mtmp[2]
#define XXdmunu    mtmp[3]
#define XdXnumu    mtmp[4]
#define XdXmunu    mtmp[5]
#define Xnu        ftmp0[0]
#define Xnufmu     ftmp[0][mu]
#define Xmu        ftmp0[1]
#define Xmufnu     ftmp[1][nu]
#define Unu        ftmp0[2]
#define Unufmu     ftmp[2][mu]
#define Umu        ftmp0[3]
#define Umufnu     ftmp[3][nu]
#define backmu     btmp0[0]
#define backmubnu  btmp[0][nu]
#define backnu     btmp0[1]
#define backnubmu  btmp[1][mu]

  // fmu += Unu * Xmufnu * Xnufmu+
  // fmu += Unubnu+ * Xmubnu * Xnufmubnu
  // fmu += Xnu * Xmufnu * Unufmu+
  // fmu += Xnubnu+ * Xmubnu * Unufmubnu
  // fnu += Umu * Xnufmu * Xmufnu+
  // fnu += Umubmu+ * Xnubmu * Xmufnubmu
  // fnu += Xmu * Xnufmu * Umufnu+
  // fnu += Xmubmu+ * Xnubmu * Umufnubmu
  QDP_M_eq_zero(XXnumu, QDP_all);
  QDP_M_eq_zero(XXmunu, QDP_all);
  QDP_M_eq_zero(XXdnumu, QDP_all);
  QDP_M_eq_zero(XXdmunu, QDP_all);
  QDP_M_eq_zero(XdXnumu, QDP_all);
  QDP_M_eq_zero(XdXmunu, QDP_all);
  for(int i=0; i<nx; i++) {
    QDP_M_eq_M(Xnu, Xnu0[i], QDP_all);
    QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, QDP_all);
    QDP_M_eq_M(Xmu, Xmu0[i], QDP_all);
    QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, QDP_all);
    QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, QDP_all);
    QDP_M_peq_Ma_times_M(XdXmunu, Xmu, Xnu, QDP_all);
    QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, QDP_all);
    QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, QDP_all);
    QDP_M_peq_M_times_Ma(XXdnumu, Xnufmu, Xmufnu, QDP_all);
    QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, QDP_all);
    QDP_discard_M(Xmufnu);
    QDP_discard_M(Xnufmu);
  }
  QDP_M_eq_M(Unu, Unu0, QDP_all);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(Umu, Umu0, QDP_all);
  QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, QDP_all);

  QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, QDP_all);
  QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, QDP_all);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_Ma_times_M(backnu, Umu, XXnumu, QDP_all);
  QDP_M_peq_M_times_M(backnu, XdXmunu, Umufnu, QDP_all);
  QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, QDP_all);

  QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, QDP_all);
  QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, QDP_all);
  QDP_M_peq_M(outmu, backmubnu, QDP_all);
  QDP_M_peq_M_times_M(outnu, Umu, XXdnumu, QDP_all);
  QDP_M_peq_M_times_Ma(outnu, XXmunu, Umufnu, QDP_all);
  QDP_M_peq_M(outnu, backnubmu, QDP_all);

  QDP_discard_M(Umufnu);
  QDP_discard_M(Unufmu);
  QDP_discard_M(backmubnu);
  QDP_discard_M(backnubmu);
#define STAPLEP_FLOPS(n) (2*198+(6*n+6)*216+2*18)

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
	 int nx, int mu, int nu)
{
#define XXnumu     mtmp[0]
#define XXmunu     mtmp[1]
#define XXdmunu    mtmp[2]
#define XdXnumu    mtmp[3]
#define Xnu        ftmp0[0]
#define Xnufmu     ftmp[0][mu]
#define Xmu        ftmp0[1]
#define Xmufnu     ftmp[1][nu]
#define Unu        ftmp0[2]
#define Unufmu     ftmp[2][mu]
#define backmu     btmp0[0]
#define backmubnu  btmp[0][nu]

  // fmu += Unu * Xmufnu * Xnufmu+
  // fmu += Unubnu+ * Xmubnu * Xnufmubnu
  // fmu += Xnu * Xmufnu * Unufmu+
  // fmu += Xnubnu+ * Xmubnu * Unufmubnu
  //QDP_M_eq_zero(XXdmunu, QDP_all);
  //QDP_M_eq_zero(XXmunu, QDP_all);
  //QDP_M_eq_zero(XXnumu, QDP_all);
  //QDP_M_eq_zero(XdXnumu, QDP_all);
  for(int i=0; i<nx; i++) {
    QDP_M_eq_M(Xnu, Xnu0[i], QDP_all);
    QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, QDP_all);
    QDP_M_eq_M(Xmu, Xmu0[i], QDP_all);
    QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, QDP_all);
    if(i==0) {
      QDP_M_eq_Ma_times_M(XdXnumu, Xnu, Xmu, QDP_all);
      QDP_M_eq_M_times_M(XXnumu, Xnu, Xmufnu, QDP_all);
      QDP_M_eq_M_times_M(XXmunu, Xmu, Xnufmu, QDP_all);
      QDP_M_eq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, QDP_all);
    } else {
      QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, QDP_all);
      QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, QDP_all);
      QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, QDP_all);
      QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, QDP_all);
    }
    QDP_discard_M(Xmufnu);
    QDP_discard_M(Xnufmu);
  }
  QDP_M_eq_M(Unu, Unu0, QDP_all);
  QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, QDP_all);

  QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, QDP_all);
  QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, QDP_all);
  QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, QDP_all);

  QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, QDP_all);
  QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, QDP_all);
  QDP_M_peq_M(outmu, backmubnu, QDP_all);

  QDP_discard_M(Unufmu);
  QDP_discard_M(backmubnu);
  //#define STAPLEP2_FLOPS(n) (198+(4*n+3)*216+18)
#define STAPLEP2_FLOPS(n) (5*198+(4*(n-1)+3)*216+18)

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
  double dtime = QOP_time();
  double nflops = 0;

  QLA_Real fac = -eps/QLA_Nc;
  QLA_Real plaq = fac*coeffs->plaquette;
  QLA_Real rect = fac*coeffs->rectangle;
  QLA_Real pgm  = fac*coeffs->parallelogram;

  QDP_ColorMatrix *tforce[4], *stplf[4][4], *stplb[4][4], *tmat[4];
  for(int mu=0; mu<4; mu++) {
    tforce[mu] = QDP_create_M();
    QDP_M_eq_zero(tforce[mu], QDP_all);
  }
  if(rect||pgm) {
    for(int mu=0; mu<4; mu++) {
      tmat[mu] = QDP_create_M();
    }
  }
  set_temps();
  for(int mu=1; mu<4; mu++) {
    for(int nu=0; nu<mu; nu++) {
      stplf[mu][nu] = QDP_create_M();
      stplb[mu][nu] = QDP_create_M();
      stplf[nu][mu] = QDP_create_M();
      stplb[nu][mu] = QDP_create_M();
      staple2fb(stplf[mu][nu], stplb[mu][nu], stplf[nu][mu], stplb[nu][mu],
		G(mu), G(nu), mu, nu);
      nflops += STAPLE2FB_FLOPS;
      if(plaq) {
	QDP_M_peq_M(tforce[mu], stplf[mu][nu], QDP_all);
	QDP_M_peq_M(tforce[mu], stplb[mu][nu], QDP_all);
	QDP_M_peq_M(tforce[nu], stplf[nu][mu], QDP_all);
	QDP_M_peq_M(tforce[nu], stplb[nu][mu], QDP_all);
	nflops += 4*18;
      }
    }
  }
  if(plaq) {
    for(int mu=0; mu<4; mu++) {
      QDP_M_eq_r_times_M(tforce[mu], &plaq, tforce[mu], QDP_all);
    }
    nflops += 4*18;
  }

  if(rect) {
    for(int mu=0; mu<4; mu++) {
      QDP_M_eq_zero(tmat[mu], QDP_all);
    }
    for(int mu=1; mu<4; mu++) {
      for(int nu=0; nu<mu; nu++) {
	//if(nu==mu) continue;
	//staples(tmat[mu], G(mu), G(mu), stplb[nu][mu], G(nu), mu, nu);
	//staples(tmat[mu], G(mu), G(mu), G(nu), stplf[nu][mu], mu, nu);
	//staples(tmat[mu], stplf[mu][nu], stplb[mu][nu], G(nu), G(nu), mu, nu);
	stapler(tmat[mu], tmat[nu], G(mu), G(nu), stplf[mu][nu], stplf[nu][mu],
		stplb[mu][nu], stplb[nu][mu], mu, nu);
	nflops += STAPLER_FLOPS;
      }
    }
    for(int mu=0; mu<4; mu++) {
      QDP_M_peq_r_times_M(tforce[mu], &rect, tmat[mu], QDP_all);
    }
    nflops += 4*36;
  }

  if(pgm) {
#if 0
    for(int mu=0; mu<4; mu++) {
      QDP_M_eq_zero(tmat[0], QDP_all);
      for(int nu=0; nu<4; nu++) {
	if(nu==mu) continue;
	for(int rho=nu+1; rho<4; rho++) {
	  if(rho==mu||rho==nu) continue;
	  staples(tmat[0], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
	  staples(tmat[0], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
	  staples(tmat[0], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
	  staples(tmat[0], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
	  nflops += 4*STAPLES_FLOPS;
	}
      }
      QDP_M_peq_r_times_M(tforce[mu], &pgm, tmat[0], QDP_all);
      nflops += 36;
    }
#else
    for(int mu=0; mu<4; mu++) {
      QDP_M_eq_zero(tmat[mu], QDP_all);
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
      staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu);
      nflops += STAPLEP_FLOPS(4);
#else
      mu=0; nu=1; rho=2;
      staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu), mu, nu);
      staples(tmat[mu], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho], mu, nu);
      staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu), mu, nu);
      staples(tmat[mu], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho], mu, nu);
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
      staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu);
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
      staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu);
      mu=1; nu=2; rho=3;
      Xmu[0] = stplf[mu][rho];
      Xmu[1] = stplb[mu][rho];
      Xnu[0] = stplf[nu][rho];
      Xnu[1] = stplb[nu][rho];
      staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu);
      mu=2; nu=0; rho=1;
      Xmu[0] = stplf[mu][rho];
      Xmu[1] = stplb[mu][rho];
      Xnu[0] = stplf[nu][rho];
      Xnu[1] = stplb[nu][rho];
      staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu);
      mu=3; nu=0; rho=1;
      Xmu[0] = stplf[mu][rho];
      Xmu[1] = stplb[mu][rho];
      Xnu[0] = stplf[nu][rho];
      Xnu[1] = stplb[nu][rho];
      staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu);
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
      QDP_M_peq_r_times_M(tforce[mu], &pgm, tmat[mu], QDP_all);
    }
    nflops += 4*36;
#endif
  }

  for(int mu=0; mu<4; mu++) {
    QDP_M_peq_M(deriv[mu], tforce[mu], QDP_all);
  }
  nflops += 4*18;

  if(rect||pgm) {
    for(int mu=0; mu<4; mu++) {
      QDP_destroy_M(tmat[mu]);
    }
  }
  for(int mu=0; mu<4; mu++) {
    QDP_destroy_M(tforce[mu]);
    for(int nu=0; nu<4; nu++) {
      if(nu==mu) continue;
      QDP_destroy_M(stplf[mu][nu]);
      QDP_destroy_M(stplb[mu][nu]);
    }
  }
  free_temps();

  //double nflop = 96720 - 4*(24+18);
  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
} 

// TODO: use top gauge field in force
void
QOP_symanzik_1loop_gauge_force_qdp(QOP_info_t *info, QOP_GaugeField *gauge,
				   QDP_ColorMatrix *force[], QOP_gauge_coeffs_t *coeffs,
				   QLA_Real eps)
{
  double dtime = QOP_time();

  QDP_ColorMatrix *deriv[4];
  for(int mu=0; mu<4; mu++) {
    deriv[mu] = QDP_create_M();
    QDP_M_eq_zero(deriv[mu], QDP_all);
  }
  QOP_symanzik_1loop_gauge_deriv_qdp(info, gauge, deriv, coeffs, eps, 0);

  QDP_ColorMatrix *mtmp = QDP_create_M();
#ifdef CHKSUM
  QLA_ColorMatrix qcm;
  QLA_Complex det, chk;
  QLA_c_eq_r(chk, 0);
#endif
  QLA_Real trace=0;
  for(int mu=0; mu<4; mu++) {
    QDP_M_eq_M_times_Ma(mtmp, G(mu), deriv[mu], QDP_all);
    if(QOP_common.verbosity==QOP_VERB_DEBUG) {
      QLA_ColorMatrix tm;
      QLA_Real tr;
      QDP_m_eq_sum_M(&tm, mtmp, QDP_all);
      QLA_R_eq_re_trace_M(&tr, &tm);
      trace += tr;
    }
    QDP_M_eq_antiherm_M(deriv[mu], mtmp, QDP_all);
    QDP_M_peq_M(force[mu], deriv[mu], QDP_all);
#ifdef CHKSUM
    QDP_m_eq_sum_M(&qcm, force->force[mu], QDP_all);
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
  info->final_flop += 4*(198+24+18)*QDP_sites_on_node; 

  QDP_destroy_M(mtmp);
  for(int mu=0; mu<4; mu++) {
    QDP_destroy_M(deriv[mu]);
  }

  info->final_sec = QOP_time() - dtime;
}

void
QOP_symanzik_1loop_gauge_force(QOP_info_t *info, QOP_GaugeField *gauge, 
			       QOP_Force *force, QOP_gauge_coeffs_t *coeffs, REAL eps)
{
  QOP_symanzik_1loop_gauge_force_qdp(info, gauge, force->force, coeffs, eps);
}
