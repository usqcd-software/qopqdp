//#define DO_TRACE
#include <qop.h>
#include <qop_internal.h>

#ifdef USE_MG

#include "solvers.h"

#define QCDPC(x) QOP_##x

#if QOP_Precision == 'F'
#define QDPN(x) QDP_FN_##x
#else
#define QDPN(x) QDP_DN_##x
#endif

//#define NV_BICGSTAB

#define EO_PREC
#define SPLIT_CHIRALITIES

//#define DO_TRACE
#include <math.h>
#include <float.h>
#include <string.h>
#include <qla_fn.h>
#include <qla_dn.h>
#include <qdp_fn.h>
#include <qdp_dn.h>
#include <qdp_dfn.h>

#if QDP_Precision == 1 || QDP_Precision == 'F'

static QDP_Lattice *theLat = NULL;

QOP_WilsonMg *
QOP_wilsonMgNew(QDP_Lattice *lat)
{
  QOP_WilsonMg *QOP_malloc(wmg, QOP_WilsonMg, 1);
  wmg->wilD = NULL;
  wmg->wilF = NULL;
  wmg->wilF_priv = NULL;
  wmg->kappa = 0;
  wmg->kappanv = 0;
  wmg->vcwaF.wil = NULL;
  wmg->vcwaF.kappa = wmg->kappa;
  wmg->nvwaF.wil = NULL;
  wmg->nvwaF.kappa = wmg->kappanv;
  wmg->nlevels = 0;
  wmg->mg = NULL;
  wmg->verbose = 0;
  wmg->profile = 0;
  wmg->itmax = 100;
  wmg->gcrF = NULL;
  wmg->gcrD = NULL;
  wmg->ngcr = 8;
  wmg->qdp_lattice = lat; /* Ugly */
#if QOP_Colors == 'N'
  wmg->nc = -1;
#else
  wmg->nc = QLA_Nc;
#endif
  return wmg;
}

static void
init_level(QOP_WilMgLevel l[], int n, QDP_Lattice *lat)
{
  QOP_printf0("init_level %i\n", n);
  l[n].ndim = 4;
  QOP_malloc(l[n].lattice_size, int, l[n].ndim);
  if(n==0) {
    l[n].lattice_size[0] = QDP_coord_size_L(lat, 0)/3;
    l[n].lattice_size[1] = QDP_coord_size_L(lat, 1)/3;
    l[n].lattice_size[2] = QDP_coord_size_L(lat, 2)/3;
    l[n].lattice_size[3] = QDP_coord_size_L(lat, 3)/8;
    l[n].nvecs = 24;
    l[n].npre = 0;
    l[n].npost = 5;
    l[n].scale = 1;
    l[n].fnc = -1;
  }
  if(n==1) {
    l[n].lattice_size[0] = l[n-1].lattice_size[0]/2;
    l[n].lattice_size[1] = l[n-1].lattice_size[1]/2;
    l[n].lattice_size[2] = l[n-1].lattice_size[2]/2;
    l[n].lattice_size[3] = l[n-1].lattice_size[3]/4;
    l[n].nvecs = 32;
    l[n].npre = 0;
    l[n].npost = 3;
    l[n].scale = 1;
    l[n].fnc = l[n-1].nvecs;
  }
  if(n>1) {
    l[n].lattice_size[0] = l[n-1].lattice_size[0]/4;
    l[n].lattice_size[1] = l[n-1].lattice_size[1]/4;
    l[n].lattice_size[2] = l[n-1].lattice_size[2]/4;
    l[n].lattice_size[3] = l[n-1].lattice_size[3]/4;
    l[n].nvecs = 10;
    l[n].npre = 0;
    l[n].npost = 2;
    l[n].scale = 1;
    l[n].fnc = l[n-1].nvecs;
  }
  for(int i=0; i<l[n].ndim; i++)
    if(l[n].lattice_size[i]<1) l[n].lattice_size[i] = 1;
#ifdef SPLIT_CHIRALITIES
  l[n].nv = 2;
#else
  l[n].nv = 1;
#endif
  l[n].cres = 0.2;
  l[n].itmax = 50;
  l[n].ngcr = 8;
  l[n].created = 0;
  l[n].setup_res = 0.4;
  l[n].setup_change_fac = 0.5;
  l[n].setup_maxit = 100;
  l[n].setup_nvecs = 0;
  l[n].verbose = 1;
  l[n].lattice = NULL;
  l[n].mgblock = NULL;
  l[n].mgargs = NULL;
  l[n].vca = NULL;
  l[n].cfoa = NULL;
  l[n].dargs = NULL;
  l[n].gcrc = NULL;
  l[n].sa = NULL;
  l[n].fcoa = NULL;
  l[n].gcrf = NULL;
}

static QLA_Real
norm2V(QDP_N_ColorVector **v, int nv)
{
  QDP_Lattice *lat = QDP_N_get_lattice_V(v[0]);
  QDP_Subset allf = QDP_all_L(lat);
  QLA_Real nrm2 = 0;
  for(int i=0; i<nv; i++) {
    QLA_Real tt;
    QDP_N_r_eq_norm2_V(&tt, v[i], allf);
    nrm2 += tt;
  }
  return nrm2;
}

static QLA_Complex
dotV(QDPN(ColorVector) **v1, QDPN(ColorVector) **v2, int nv)
{
  QDP_Lattice *lat = QDPN(get_lattice_V)(v1[0]);
  QDP_Subset allf = QDP_all_L(lat);
  QLA_Complex dot;
  QLA_c_eq_r(dot, 0);
  for(int i=0; i<nv; i++) {
    QLA_Complex tt;
    QDPN(c_eq_V_dot_V)(&tt, v1[i], v2[i], allf);
    QLA_c_peq_c(dot, tt);
  }
  return dot;
}

static void
print_norms(QDP_FN_ColorVector *vv[], int nv, int fnc, int cnc, 
            void op(QDP_FN_ColorVector *out[], QDP_FN_ColorVector *in[], int sign, void *args),
            void *opargs, int sign)
{
  QDP_FN_ColorVector *(*v)[cnc] = (QDP_FN_ColorVector *(*)[cnc]) vv;
  QDP_Lattice *fine = QDP_FN_get_lattice_V(v[0][0]);
  //QDP_Subset allf = QDP_all_L(fine);
  QDP_FN_ColorVector *cv[nv], *tv[nv];
  for(int j=0; j<nv; j++) {
    tv[j] = QDP_FN_create_V_L(fnc, fine);
  }
  for(int i=0; i<cnc; i++) {
    QLA_F_Real nrm, nrm2;
#if 0
    nrm = 0;
    for(int j=0; j<nv; j++) {
      QLA_F_Real tt;
      cv[j] = v[j][i];
      QDP_FN_r_eq_norm2_V(&tt, cv[j], allf);
      nrm += tt;
    }
#endif
    for(int j=0; j<nv; j++) cv[j] = v[j][i];
    nrm = norm2V(cv, nv);
    QOP_printf0("%-3i  %-10g", i, nrm);
    op(tv, cv, sign, opargs);
#if 0
    nrm2 = 0;
    for(int j=0; j<nv; j++) {
      QLA_F_Real tt;
      QDP_FN_r_eq_norm2_V(&tt, tv[j], allf);
      nrm2 += tt;
    }
#endif
    nrm2 = norm2V(tv, nv);
    QOP_printf0("  %-10g  %-10g\n", nrm2, sqrt(nrm2/nrm));
  }
  for(int i=0; i<nv; i++) QDP_FN_destroy_V(tv[i]);
}

static void
lex_int(QLA_Int *li, int coords[])
{
  int t = coords[0];
  for(int i=1; i<QDP_ndim_L(theLat); i++) {
    t = t*QDP_coord_size_L(theLat, i) + coords[i];
  }
  *li = t;
}

static void
QOP_randSeed(QDP_RandomState *rs, int seed)
{
  QDP_Int *li;
  QDP_Lattice *lat = QDP_get_lattice_S(rs);
  QDP_Subset all = QDP_all_L(lat);

  li = QDP_create_I_L(lat);

  TRACE;
  theLat = lat;
  QDP_I_eq_funct(li, lex_int, all);
  theLat = NULL;
  TRACE;
  QDP_S_eq_seed_i_I(rs, seed, li, all);
  TRACE;

  QDP_destroy_I(li);
}

static void
get_nullvecs(QOP_F_MgArgs *mgargs, QOPP(MgOp) *op, void *opargs, QLA_Real res,
	     QLA_Real change_fac, int maxit, int nvecs, int verbose,
	     int (smoother)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in,
			    int sign, void *args), void *sargs, int update)
{
  int nv = mgargs->nv;
  int fnc = mgargs->fnc;
  int cnc = mgargs->cnc;
  QDP_Lattice *fine = mgargs->mgb->fine;
  QDP_Subset allf = QDP_all_L(fine);
  TRACE;
  QDP_RandomState *rs = QDP_create_S_L(fine); // hack: lives on fine lattice
  TRACE;
  QOP_randSeed(rs, 987654321); // needs to be fixed to use any lattice
  TRACE;
  double dt = -QDP_time();
  int cnc0 = cnc;
  if(nvecs>cnc) cnc = nvecs;

  QDP_FN_ColorVector *pv[cnc][nv];
  QDP_FN_ColorVector *mpv[nv][cnc];
  for(int i=0; i<cnc; i++) {
    for(int j=0; j<nv; j++) {
      if(i<cnc0) {
	pv[i][j] = mgargs->pv[j][i];
      } else {
	pv[i][j] = QDP_FN_create_V_L(fnc, fine);
      }
      mpv[j][i] = pv[i][j];
    }
  }
  TRACE;

  QDP_FN_ColorVector *tv[nv], *tv2[nv];
  for(int j=0; j<nv; j++) tv[j] = QDP_FN_create_V_L(fnc, fine);
  for(int j=0; j<nv; j++) tv2[j] = QDP_FN_create_V_L(fnc, fine);
  TRACE;

#define prn(v) {QLA_Real n2 = norm2V(v, nv); QOP_printf0("norm2 " #v ": %g\n", n2); }

  int tnits=0;
  int maxhits=20;
  double svhist[maxhits];
  for(int i=0; i<cnc; i++) {
    if(i==0 && !update) {
      for(int j=0; j<nv; j++) QDP_FN_V_eq_gaussian_S(pv[i][j], rs, allf);
    }
    TRACE;
    //if(i==1) maxit /= 2;
    QLA_Real sv2=FLT_MAX, sv2o;
    int nhits=0;
    int nits=0;
    do {
      nhits++;
      //QOP_mgOrthoVn(mpv, nv, 0, i+1, allf);
      //prn(pv[i]);
      for(int k=0; k<i; k++) {
        QLA_Complex z = dotV(pv[k], pv[i], nv);
        for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(pv[i][j], &z, pv[k][j], allf); }
      }
      TRACE;
      //prn(pv[i]);
      nits += smoother(tv, pv[i], 1, sargs);
      TRACE;
      //prn(tv);
      QLA_Real nrm = norm2V(tv, nv);
      TRACE;
      op(pv[i], tv, 1, opargs);
      TRACE;
      QLA_Real nrm2 = norm2V(pv[i], nv);
      sv2o = sv2;
      sv2 = nrm2/nrm;
      svhist[nhits-1] = sv2;
      QLA_Real scale = 1/sqrt(nrm);
      if(i+1<cnc && !update) for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(pv[i+1][j], pv[i][j], allf); }
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_r_times_V(pv[i][j], &scale, tv[j], allf); }
      if(verbose>0) QOP_printf0("%i %g\n", i, sv2);
      if(i>0 && nhits==4) break;
      TRACE;
    } while(sv2<change_fac*sv2o);
    tnits += nits;
    QOP_printf0("%-3i %-3i %-5i  %-10g  %-10g", i, nhits, nits, sqrt(sv2), sqrt(sv2o));
    for(int i=nhits-3; i>0; i--) QOP_printf0("  %-10g", sqrt(svhist[i]));
    QOP_printf0("\n");
  }
  QOP_printf0("setup its = %i\n", tnits);

  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  QOP_F_mgOrthoVn(mgargs->pv, nv, 0, cnc0, allf);
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int j=0; j<nv; j++) {
    QOP_F_mgOrtho(mpv[j], cnc, mgargs->mgb);
  }
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int i=0; i<cnc0; i++) {
    for(int j=0; j<nv; j++) {
#ifndef SPLIT_CHIRALITIES
      if(fnc==4*QLA_Nc) { // hack -- top level without splitting chiralities
	QDP_F_D_eq_gamma_times_D((QDP_F_DiracFermion*)mgargs->rv[j][i], (QDP_F_DiracFermion*)mgargs->pv[j][i], 15, allf);
      } else
#endif
	QDP_FN_V_eq_V(mgargs->rv[j][i], mgargs->pv[j][i], allf);
    }
  }
  dt += QDP_time();
  QOP_printf0("%s: %g secs\n", __func__, dt);

  for(int i=0; i<nv; i++) {
    QDP_FN_V_eq_gaussian_S(tv[i], rs, allf);
  }
  QOP_F_mgTestCoarse(mgargs, tv);
  for(int i=0; i<nv; i++) {
    QDP_FN_destroy_V(tv[i]);
    QDP_FN_destroy_V(tv2[i]);
  }
  for(int i=cnc0; i<cnc; i++) {
    for(int j=0; j<nv; j++) {
      QDP_FN_destroy_V(pv[i][j]);
    }
  }
  QDP_destroy_S(rs);
}

static void
create_level(QOP_WilsonMg *wmg, int n)
{
#define NC (wmg->nc)
  QOP_printf0("creating MG level %i\n", n);

  QOP_WilMgLevel *l = wmg->mg;
  QDP_Lattice *lat0;
  int nv = l[n].nv;
  int update = 0;

  wmg->nvwaF.kappa = wmg->kappanv;
  wmg->vcwaF.kappa = wmg->kappa;

  void (*smoothop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  if(n==0) {
    if(l[n].fnc<=0) {
      if(wmg->nc<=0) {
	QOP_printf0("%s error: wmg->nc not set\n", __func__);
	QDP_abort(-1);
      } else {
#ifdef SPLIT_CHIRALITIES
	l[n].fnc = 2*wmg->nc;
#else
	l[n].fnc = 4*wmg->nc;
#endif
      }
    }

    lat0 = wmg->qdp_lattice;
#ifdef SPLIT_CHIRALITIES
#ifdef EO_PREC
    l[n].vcop = QOPFC(wilEoV2);
    l[n].nvop = QOPFC(wilPV2);
    l[n].cfop = QOPFC(wilPV2);
#else
    l[n].vcop = QOPFC(wilDV2);
    l[n].nvop = QOPFC(wilDV2);
    l[n].cfop = QOPFC(wilDV2);
#endif
    smoothop = QOPFC(wilEoV2);
#else
#ifdef EO_PREC
    l[n].vcop = QOPFC(wilEoV1);
    l[n].nvop = QOPFC(wilPV1);
    l[n].cfop = QOPFC(wilPV1);
#else
    l[n].vcop = QOPFC(wilDV1);
    l[n].nvop = QOPFC(wilDV1);
    l[n].cfop = QOPFC(wilDV1);
#endif
    smoothop = QOPFC(wilEoV1);
#endif
    l[n].vcopargs = &wmg->vcwaF;
    l[n].nvopargs = &wmg->nvwaF;
    l[n].cfopargs = &wmg->vcwaF;
  } else {
    lat0 = l[n-1].lattice;
#ifdef EO_PREC
    l[n].vcop = QOP_F_mgDslashEo;
    l[n].nvop = QOP_F_mgDslashPAll;
    l[n].cfop = QOP_F_mgDslashPAll;
#else
    l[n].vcop = QOP_F_mgDslashAll;
    l[n].nvop = QOP_F_mgDslashAll;
    l[n].cfop = QOP_F_mgDslashAll;
#endif
    smoothop = QOP_F_mgDslashEo;
    l[n].vcopargs = l[n-1].dargs;
    l[n].nvopargs = l[n-1].dargs;
    l[n].cfopargs = l[n-1].dargs;
  }

  // create coarse lattice, MgBlock and MgArgs
  if(l[n].lattice==NULL) {
    QOP_printf0("creating lattice:");
    for(int i=0; i<l[n].ndim; i++) QOP_printf0(" %i", l[n].lattice_size[i]);
    QOP_printf0("\n");
    l[n].lattice = QDP_create_lattice(QDP_get_default_layout(), NULL,
				      l[n].ndim, l[n].lattice_size);
  }
  if(l[n].mgblock==NULL) {
    QOP_printf0("\ncreating block\n");
    l[n].mgblock = QOP_mgCreateBlockFromLattice(lat0, l[n].lattice);
  }
  if(l[n].mgargs==NULL) {
    QOP_printf0("creating args\n");
    l[n].mgargs = QOP_mgCreateArgs(l[n].mgblock, l[n].nvecs, l[n].fnc, l[n].nv, NULL, NULL);
  } else {
    update = 1;
  }

  // create some vectors for the vcycle and also used in setup
  QDP_FN_ColorVector **fr, **fp, **fAp;
  if(l[n].vca==NULL) {
    QOP_malloc(fr, QDP_FN_ColorVector *, nv);
    QOP_malloc(fp, QDP_FN_ColorVector *, nv);
    QOP_malloc(fAp, QDP_FN_ColorVector *, nv);
    for(int i=0; i<nv; i++) {
      fr[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
      fp[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
      fAp[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
    }
  } else {
    fr = l[n].vca->r;
    fp = l[n].vca->p;
    fAp = l[n].vca->Ap;
  }

  // get nullvecs
  QOP_printf0("creating null vecs\n");
  {
    TRACE;
    QOP_Bicgstab *bcg = QOP_bicgstabInit(lat0, l[n].nv, l[n].fnc);
    TRACE;
    QOP_bicgstabSet(bcg, "verbose", l[n].verbose);
    QOP_bicgstabSet(bcg, "indent", 2);
    QOP_F_BicgstabSolveArgs *QOP_malloc(nvsa, QOP_F_BicgstabSolveArgs, 1);
    //nvsa->op = l[n].vcop;
    nvsa->op = smoothop;
    //nvsa->opargs = l[n].vcopargs;
    //nvsa->op = l[n].nvop;
    nvsa->opargs = l[n].nvopargs;
    nvsa->pop = NULL;
    nvsa->popargs = NULL;
    nvsa->res = l[n].setup_res;
    nvsa->itmax = l[n].setup_maxit;
    nvsa->bicgstab = bcg;
    //#ifdef EO_PREC
    nvsa->sub = QDP_even_L(lat0);
    //#else
    //nvsa->sub = QDP_all_L(lat0);
    //#endif
    nvsa->ineo = fr;
    nvsa->outeo = fp;
    if(n==0) {
#ifdef EO_PREC
#ifdef SPLIT_CHIRALITIES
      nvsa->project = QOP_wilEoProjectV2;
      nvsa->reconstruct = QOP_wilEoReconstructPV2;
#else
      nvsa->project = QOP_wilEoProjectV1;
      nvsa->reconstruct = QOP_wilEoReconstructPV1;
#endif
#else
#ifdef SPLIT_CHIRALITIES
      nvsa->project = QOP_wilEoProjectV2;
      nvsa->reconstruct = QOP_wilEoReconstructV2;
#else
      nvsa->project = QOP_wilEoProjectV1;
      nvsa->reconstruct = QOP_wilEoReconstructV1;
#endif
#endif
    } else {
#ifdef EO_PREC
      nvsa->project = QOP_F_mgDslashEoProject;
      nvsa->reconstruct = QOP_F_mgDslashEoReconstructP;
#else
      nvsa->project = QOP_F_mgDslashEoProject;
      nvsa->reconstruct = QOP_F_mgDslashEoReconstruct;
#endif
    }
    printf("%p %p\n", l[n].nvop, QOPFC(wilPV2));
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose);
    //#ifdef EO_PREC
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveEo, nvsa);
  TRACE;
  get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveEo, nvsa, update);
  TRACE;
    //#else
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveA, nvsa);
    //#endif
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveA, nvsa);
    //get_nvev(l[n].mgargs, l[n].nvop, l[n].nvopargs, QDP_all_L(lat0), l[n].vcop, l[n].vcopargs, QDP_even_L(lat0), l[n].setup_nvecs, l[n].setup_maxit, l[n].setup_change_fac);
    //l[n].nvecs -= 2;
    //l[n].mgargs->nv -= 2;
    QOP_free(nvsa);
    QOP_bicgstabFree(bcg);
  }
  TRACE;

  // create coarse Dslash
  if(l[n].cfoa==NULL) {
    QDP_FN_ColorVector **QOP_malloc(fin, QDP_FN_ColorVector *, nv);
    QDP_FN_ColorVector **QOP_malloc(fout, QDP_FN_ColorVector *, nv);
    for(int i=0; i<nv; i++) {
      fin[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
      fout[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
    }
    QOP_F_MgC2fOpArgs *QOP_malloc(cfoa, QOP_F_MgC2fOpArgs, 1);
    l[n].cfoa = cfoa;
    cfoa->mga = l[n].mgargs;
    cfoa->op = l[n].cfop;
    cfoa->opargs = l[n].cfopargs;
    cfoa->fin = fin;
    cfoa->fout = fout;
    cfoa->fpar = QOP_EVENODD;
  }

  TRACE;
  // create true coarse op
  if(l[n].dargs==NULL) {
    QOP_printf0("\ncreating true coarse operator\n");
#if 0
    QOP_F_MgDslashArgs *dargs =
      QOP_F_mgCreateDslash(l[n].nv, l[n].nvecs, l[n].lattice);
    QOP_F_mgCloneOp(QOP_F_mgC2fOp, cfoa, 1, dargs);
#else
    int npaths = 9;
    int paths[npaths][4];
    for(int j=0; j<4; j++) paths[0][j] = 0;
    for(int i=0; i<4; i++) {
      for(int j=0; j<4; j++) paths[i+1][j] = 0;
      paths[i+1][i] = 1;
      for(int j=0; j<4; j++) paths[i+5][j] = 0;
      paths[i+5][i] = -1;
    }
    QOP_F_MgDslashArgs *dargs =
      QOP_F_mgCreateDslash(l[n].nv, l[n].nv, l[n].nvecs,
			   npaths, paths[0], l[n].lattice);
#endif
    l[n].dargs = dargs;
  }
  TRACE;
  QOP_F_mgCloneOp(QOP_F_mgC2fOp, l[n].cfoa, l[n].dargs);
  TRACE;

  // coarse solver
  if(l[n].gcrc==NULL) {
    QOP_F_Gcr *gcrc = QOP_F_gcrInit(l[n].lattice, nv, l[n].nvecs, l[n].ngcr);
    QOP_F_gcrSet(gcrc, "verbose", l[n].verbose);
    QOP_F_gcrSet(gcrc, "indent", 2*(n+2));
    l[n].gcrc = gcrc;
  }
  if(l[n].sa==NULL) {
    QDP_FN_ColorVector **QOP_malloc(cv2ineo, QDP_FN_ColorVector *, nv);
    QDP_FN_ColorVector **QOP_malloc(cv2outeo, QDP_FN_ColorVector *, nv);
    for(int i=0; i<nv; i++) {
      cv2ineo[i] = QDP_FN_create_V_L(l[n].nvecs, l[n].lattice);
      cv2outeo[i] = QDP_FN_create_V_L(l[n].nvecs, l[n].lattice);
    }
    QOP_F_GcrSolveArgs *QOP_malloc(sa, QOP_F_GcrSolveArgs, 1);
    l[n].sa = sa;
#ifdef EO_PREC
    sa->op = QOP_F_mgDslashEo;
    sa->sub = QDP_even_L(l[n].lattice);
#else
    sa->op = QOP_F_mgDslashAll;
    sa->sub = QDP_all_L(l[n].lattice);
#endif
    sa->opargs = l[n].dargs;
    sa->pop = NULL;
    sa->popargs = NULL;
    sa->res = l[n].cres;
    sa->itmax = l[n].itmax;
    sa->ngcr = l[n].ngcr;
    sa->gcr = l[n].gcrc;
    sa->ineo = cv2ineo;
    sa->outeo = cv2outeo;
    sa->project = QOP_F_mgDslashEoProject;
    sa->reconstruct = QOP_F_mgDslashEoReconstruct;
  }

  // fine to coarse solver
  if(l[n].fcoa==NULL) {
    QDP_FN_ColorVector **QOP_malloc(cin, QDP_FN_ColorVector *, nv);
    QDP_FN_ColorVector **QOP_malloc(cout, QDP_FN_ColorVector *, nv);
    for(int i=0; i<nv; i++) {
      cin[i] = QDP_FN_create_V_L(l[n].nvecs, l[n].lattice);
      cout[i] = QDP_FN_create_V_L(l[n].nvecs, l[n].lattice);
    }
    QOP_MgF2cOpArgs *QOP_malloc(mgoa2, QOP_MgF2cOpArgs, 1);
    l[n].fcoa = mgoa2;
    mgoa2->mga = l[n].mgargs;
#ifdef EO_PREC
    mgoa2->op = QOP_F_gcrSolveEo;
    mgoa2->fpar = QOP_EVEN;
#else
    mgoa2->op = QOP_F_gcrSolveA;
    mgoa2->fpar = QOP_EVENODD;
#endif
    mgoa2->opargs = l[n].sa;
    mgoa2->cin = cin;
    mgoa2->cout = cout;
  }

  // smoother and vcycle
  if(l[n].gcrf==NULL) {
#define MAX(a,b) ((a)<(b)?(b):(a));
    int ngcrf = MAX(l[n].npre, l[n].npost);
    QOP_F_Gcr *gcrf = QOP_F_gcrInit(lat0, nv, l[n].fnc, ngcrf);
    QOP_F_gcrSet(gcrf, "verbose", l[n].verbose);
    QOP_F_gcrSet(gcrf, "indent", 2*(n+1));
    l[n].gcrf = gcrf;
  }
  if(l[n].vca==NULL) {
    QOP_F_MgVcycleArgs *QOP_malloc(vc, QOP_F_MgVcycleArgs, 1);
    l[n].vca = vc;
    vc->op = l[n].vcop;
    vc->opargs = l[n].vcopargs;
    vc->cop = QOP_F_mgF2cOp;
    vc->copargs = l[n].fcoa;
    vc->nv = l[n].nv;
    vc->npre = l[n].npre;
    vc->npost = l[n].npost;
    vc->verbose = l[n].verbose;
    vc->indent = 2*(n+1);
#ifdef EO_PREC
    vc->sub = QDP_even_L(lat0);
#else
    vc->sub = QDP_all_L(lat0);
#endif
    vc->gcr = l[n].gcrf;
    vc->r = fr;
    vc->p = fp;
    vc->Ap = fAp;
    vc->s = l[n].scale;
  }

  if(n>0) {
    // need to set vc as preconditioner for prev solver    
    l[n-1].sa->pop = QOP_F_mgVcycle;
    l[n-1].sa->popargs = l[n].vca;
  }
#undef NC
}

static void
free_level(QOP_WilMgLevel *l)
{
  //QOP_printf0("*** freeing level %p ***\n", l);
  QOP_free(l->lattice_size);
  // FIXME

  if (l->sa) {
    if (l->sa->ineo) {
      for(int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->sa->ineo[i]);
      QOP_free(l->sa->ineo);
    }
    if (l->sa->outeo) {
      for(int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->sa->outeo[i]);
      QOP_free(l->sa->outeo);
    }
    QOP_free(l->sa);
    l->sa = NULL;
  }

  if (l->cfoa) {
    if (l->cfoa->fin) {
      for(int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->cfoa->fin[i]);
      QOP_free(l->cfoa->fin);
    }
    if (l->cfoa->fout) {
      for(int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->cfoa->fout[i]);
      QOP_free(l->cfoa->fout);
    }
    QOP_free(l->cfoa);
    l->cfoa = NULL;
  }

  if (l->fcoa) {
    if (l->fcoa->cin) {
      for (int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->fcoa->cin[i]);
      QOP_free(l->fcoa->cin);
    }
    if (l->fcoa->cout) {
      for (int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->fcoa->cout[i]);
      QOP_free(l->fcoa->cout);
    }
    QOP_free(l->fcoa);
    l->fcoa = NULL;
  }

  if (l->vca) {
    if (l->vca->r) {
      for (int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->vca->r[i]);
      QOP_free(l->vca->r);
    }

    if (l->vca->p) {
      for (int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->vca->p[i]);
      QOP_free(l->vca->p);
    }
    if (l->vca->Ap) {
      for (int i = 0 ; i < l->nv ; i++) QDP_FN_destroy_V(l->vca->Ap[i]);
      QOP_free(l->vca->Ap);
    }
    QOP_free(l->vca);
    l->vca = NULL;
  }

  if (l->gcrc) {
    QOP_F_gcrFree(l->gcrc);
    l->gcrc = NULL;
  }

  if (l->gcrf) {
    QOP_F_gcrFree(l->gcrf);
    l->gcrf = NULL;
  }
  if (l->mgargs) {
    QOP_F_mgFreeArgs(l->mgargs);
    l->mgargs = NULL;
  }

  if (l->dargs) {
    QOP_F_mgFreeDslash(l->dargs);
    l->dargs = NULL;
  }

  if (l->mgblock) {
    QOP_mgFreeBlock(l->mgblock);    // also destroys l->lattice
    l->mgblock = NULL;
  }
  // sic! l->lattice is destroyed by QOP_mgFreeBlock
//  if (l->lattice) {
//    QDP_destroy_lattice(l->lattice); 
//    l->lattice = NULL;
//  }
}

static void
setNumLevels(QOP_WilsonMg *wmg, int nlevels)
{
  //int ret = wmg->nlevels;
  QDP_Lattice *lat = wmg->qdp_lattice;

  if(nlevels>0) {
    for(int i=nlevels; i<wmg->nlevels; i++) {
      free_level(&wmg->mg[i]);
    }
    wmg->mg = realloc(wmg->mg, nlevels*sizeof(QOP_WilMgLevel));
    for(int i=wmg->nlevels; i<nlevels; i++) {
      init_level(wmg->mg, i, lat);
    }
    wmg->nlevels = nlevels;
  }
  //return ret;
}

void
QOP_wilsonMgFree(QOP_WilsonMg *wmg)
{
  if(wmg->mg) {
    for(int i=0; i<wmg->nlevels; i++) {
      free_level(&wmg->mg[i]);
    }
    QOP_free(wmg->mg);
  }
  if(wmg->wilF_priv) {
    QOP_F_wilson_destroy_L(wmg->wilF_priv);
  }
  //QOP_gcrFree(gcro);
  if (wmg->gcrF) {
    QOP_F_gcrFree(wmg->gcrF);
    wmg->gcrF = NULL;
  }
  if (wmg->gcrD) {
    QOP_D_gcrFree(wmg->gcrD);
    wmg->gcrD = NULL;
  }
  QOP_free(wmg);
}

#define seti(t) if(!strcmp(s,#t)) wmg->mg[l].t = (int) val;
#define setd(t) if(!strcmp(s,#t)) wmg->mg[l].t = val;
#define setalli(t) if(!strcmp(s,#t)) for(int i=0; i<wmg->nlevels; i++) wmg->mg[i].t = (int) val;
void
QOP_wilsonMgSet(QOP_WilsonMg *wmg, int l, char *s, double val)
{
  if(l==-2) { // apply to all
    setalli(verbose);
  } else if(l==-1) { // global options
    if(!strcmp(s,"nlevels")) setNumLevels(wmg, (int)val);
    if(!strcmp(s,"kappa")) wmg->kappa = val;
    if(!strcmp(s,"kappanv")) wmg->kappanv = val;
    if(!strcmp(s,"verbose")) wmg->verbose = (int) val;
    if(!strcmp(s,"profile")) wmg->profile = (int) val;
    if(!strcmp(s,"itmax")) wmg->itmax = (int) val;
    if(!strcmp(s,"ngcr")) wmg->ngcr = (int) val;
    if(!strcmp(s,"nc")) wmg->nc = (int) val;
  } else if(l>=0 || l<wmg->nlevels) {
    seti(verbose);
    seti(nvecs);
    seti(npre);
    seti(npost);
    setd(scale);
    setd(cres);
    seti(itmax);
    seti(ngcr);
    setd(setup_res);
    setd(setup_change_fac);
    seti(setup_maxit);
    seti(setup_nvecs);
  }
}

void
QOP_wilsonMgSetArray(QOP_WilsonMg *wmg, int l, char *s, double *vals, int nval)
{
  if(l==-1) { // global options
  } else if(l>=0 || l<wmg->nlevels) {
    if(!strcmp(s,"lattice")) {
      wmg->mg[l].ndim = nval;
      QOP_free(wmg->mg[l].lattice_size);  // free default
      QOP_malloc(wmg->mg[l].lattice_size, int, nval);
      for(int i=0; i<nval; i++) wmg->mg[l].lattice_size[i] = (int) vals[i];
    }
  }
}

void
QOP_wilsonMgSetup(QOP_WilsonMg *wmg)
{
  for(int i=0; i<wmg->nlevels; i++) {
    if(!wmg->mg[i].created) {
      create_level(wmg, i);
      wmg->mg[i].created = 1;
    }
  }
}

void
QOP_wilsonMgSetLinks(QOP_WilsonMg *wmg, QOP_FermionLinksWilson *wil)
{
  wmg->wilF = wil;
  wmg->vcwaF.wil = wil;
  wmg->nvwaF.wil = wil;
}

#if 0
void
QOP_wilsonMgSetVecs(QOP_WilsonMg *wmg, int nv, QDP_DiracFermion *vv[nv])
{
  QDP_Lattice *lat = QDP_get_lattice_D(vv[0]);
  QDP_Subset all = QDP_all_L(lat);
  QOP_malloc(wmg->v, QDP_DiracFermion *, nv);
  wmg->nv = nv;
  for(int i=0; i<nv; i++) {
    wmg->v[i] = QDP_create_D_L(lat);
    QDP_D_eq_D(wmg->v[i], vv[i], all);
  }
  //mgorth(wmg, vv);
}
#endif

#else // double precision

void
QOP_wilsonMgSetLinks(QOP_WilsonMg *wmg, QOP_FermionLinksWilson *wild)
{
  QOP_F_FermionLinksWilson *wil = QOP_FD_wilson_create_L_from_L(wild);
  if (wmg->wilF_priv)
    QOP_F_wilson_destroy_L(wmg->wilF_priv);
  wmg->wilF_priv = wil;
  wmg->wilF = wil;
  wmg->vcwaF.wil = wil;
  wmg->nvwaF.wil = wil;
}

#endif // single/double precision

#if 0
void
QOP_wilsonMgProj(QOP_WilsonMg *wmg, QDP_DiracFermion *x, QDP_DiracFermion *in)
{
}

void
QOP_wilsonMgApply(QOP_WilsonMg *wmg, QDP_DiracFermion *x, QDP_DiracFermion *in)
{
}
#endif

#if QDP_Precision == 2 || QDP_Precision == 'D'
typedef struct {
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  QDP_FN_ColorVector **in, **out;
  int nv;
  QDP_Subset sub;
} d2fArgs;

// scales input to unit norm, calls op, scales output back
static void
d2f(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  d2fArgs *a = (d2fArgs *)args;
  QLA_D_Real n2in=0;
  for(int i=0; i<a->nv; i++) {
    QLA_D_Real r=0;
    QDP_DN_r_eq_norm2_V(&r, in[i], a->sub);
    n2in += r;
  }
  QLA_D_Real so = sqrt(n2in);
  QLA_D_Real si = 1/so;
  for(int i=0; i<a->nv; i++) {
    QDP_DN_V_eq_r_times_V(out[i], &si, in[i], a->sub);
    QDP_FDN_V_eq_V((a->in)[i], out[i], a->sub);
  }
  a->op(a->out, a->in, sign, a->opargs);
  for(int i=0; i<a->nv; i++) {
    QDP_DFN_V_eq_V(out[i], (a->out)[i], a->sub);
    QDP_DN_V_eq_r_times_V(out[i], &so, out[i], a->sub);
  }
}
#endif

void
QOP_wilsonMgSolve(QOP_info_t *info, QOP_WilsonMg *wmg,
		  QOP_FermionLinksWilson *flw, QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg, QLA_Real kappa,
		  QDP_DiracFermion *out, QDP_DiracFermion *in)
{
#define NC QDP_get_nc(out)
  QDP_Lattice *lat = QDP_get_lattice_D(in);
  QDP_Subset all = QDP_all_L(lat);
  QOP_WilArgs wa;
  wa.wil = flw;
  wa.kappa = kappa;
  wmg->vcwaF.kappa = kappa;
  for(int i=0; i<wmg->nlevels; i++) {
    wmg->mg[i].vca->verbose = wmg->mg[i].verbose;
    QOP_F_gcrSet(wmg->mg[i].vca->gcr, "verbose", wmg->mg[i].verbose);
    wmg->mg[i].vca->npre = wmg->mg[i].npre;
    wmg->mg[i].vca->npost = wmg->mg[i].npost;
    wmg->mg[i].vca->s = wmg->mg[i].scale;
    QOP_F_gcrSet(wmg->mg[i].sa->gcr, "verbose", wmg->mg[i].verbose);
    wmg->mg[i].sa->res = wmg->mg[i].cres;
    wmg->mg[i].sa->itmax = wmg->mg[i].itmax;
    wmg->mg[i].sa->ngcr = wmg->mg[i].ngcr;
    wmg->mg[i].vca->tpre = 0;
    wmg->mg[i].vca->tcoarse = 0;
    wmg->mg[i].vca->tpost = 0;
    wmg->mg[i].vca->count = 0;
    wmg->mg[i].mgblock->trestrict = 0;
    wmg->mg[i].mgblock->tprolong = 0;
    wmg->mg[i].mgblock->rcount = 0;
    wmg->mg[i].mgblock->pcount = 0;
  }

  int nv = wmg->mg[0].nv;
  int fnc = wmg->mg[0].fnc;
  int itmax = wmg->itmax;
  int ngcr = wmg->ngcr;
  QDP_Lattice *lat0 = QDP_get_lattice_D(in);
#ifdef EO_PREC
  QDP_Subset ssub = QDP_even_L(lat);
#else
  QDP_Subset ssub = QDP_all_L(lat);
#endif

  QDP_DiracFermion *r = QDP_create_D_L(lat);
  QDP_DiracFermion *out2 = QDP_create_D_L(lat);
  QDP_DiracFermion *eovec = QDP_create_D_L(lat);
  QDPN(ColorVector) *v2in[nv], *v2out[nv];
  for(int i=0; i<nv; i++) {
    v2in[i] = QDPN(create_V_L)(fnc,lat);
    v2out[i] = QDPN(create_V_L)(fnc,lat);
  }

  QOPP(MgOp) *pop;
  void *popargs;
#if QDP_Precision == 1 || QDP_Precision == 'F'
  pop = QOP_F_mgVcycle;
  popargs = wmg->mg[0].vca;
#else
  QDP_FN_ColorVector *v2inf[nv], *v2outf[nv];
  for(int i=0; i<nv; i++) {
    v2inf[i] = QDP_FN_create_V_L(fnc,lat);
    v2outf[i] = QDP_FN_create_V_L(fnc,lat);
  }
  d2fArgs da;
  da.op = QOP_F_mgVcycle;
  da.opargs = wmg->mg[0].vca;
  da.in = v2inf;
  da.out = v2outf;
  da.nv = nv;
  da.sub = ssub;
  pop = d2f;
  popargs = &da;
#endif

#if QDP_Precision == 'F' || QDP_Precision == 1
  if(wmg->gcrF==NULL) {
    wmg->gcrF = QOP_gcrInit(lat0, nv, fnc, ngcr);
  }
  QOP_Gcr *gcro = wmg->gcrF;
#else
  if(wmg->gcrD==NULL) {
    wmg->gcrD = QOP_gcrInit(lat0, nv, fnc, ngcr);
  }
  QOP_Gcr *gcro = wmg->gcrD;
#endif
  QOP_gcrSet(gcro, "verbose", wmg->verbose);
  QOP_gcrSet(gcro, "indent", 0);
  //QOP_gcrSet(gcro, "nsmooth", 3);
  //QOP_gcrSet(gcro, "reuse", 1);

  QLA_Real norm2in, rsqstop, rsq;
  QDP_r_eq_norm2_D(&norm2in, in, all);
  rsqstop = res_arg->rsqmin*norm2in;
  double gcrres = -sqrt(res_arg->rsqmin*norm2in);

  double secs=0;
  int its = 0;
  while(1) {
    QOP_wilsonDslash(r, out, flw, kappa, 1, QOP_EVENODD, QOP_EVENODD);
    QDP_D_eq_D_minus_D(r, in, r, all);
    QDP_r_eq_norm2_D(&rsq, r, all);
    if(its>0) QOP_printf0("iters = %4i  secs = %8f  rsq = %g\n", its, secs, rsq/norm2in);
    if(rsq<rsqstop || its>=itmax) break;

#ifdef EO_PREC
    QOP_wilEoProjectD(eovec, r, &wa);
#else
    QDP_D_eq_D(eovec, r, all);
#endif
#ifdef SPLIT_CHIRALITIES
    QOP_V2eqD(v2in, eovec, ssub);
    //QOP_V2eqD(v2out, out);
#else
    QOP_V1eqD(v2in, eovec, ssub);
#endif
    { QDP_Subset sub=ssub; V_eq_zero(v2out); }

    secs -= QDP_time();
#ifdef EO_PREC
#ifdef SPLIT_CHIRALITIES
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilEoV2, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_even_L(lat));
#else
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilEoV1, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_even_L(lat));
#endif
#else
#ifdef SPLIT_CHIRALITIES
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilDV2, &wa, pop, popargs, gcrres, itmax, ngcr, all);
#else
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilDV1, &wa, pop, popargs, gcrres, itmax, ngcr, all);
#endif
#endif
    secs += QDP_time();

#ifdef SPLIT_CHIRALITIES
    QOP_DeqV2(eovec, v2out, ssub);
#else
    QOP_DeqV1(eovec, v2out, ssub);
#endif
#ifdef EO_PREC
    //QOP_wilEoReconstructD(out, eovec, in, &wa);
    QOP_wilEoReconstructD(out2, eovec, r, &wa);
#else
    QDP_D_eq_D(out2, eovec, all);
#endif
    QDP_D_peq_D(out, out2, all);
  }

  //QOP_gcrFree(gcro);

  QDP_destroy_D(r);
  QDP_destroy_D(out2);
  QDP_destroy_D(eovec);
  for(int i=0; i<nv; i++) {
    QDPN(destroy_V)(v2in[i]);
    QDPN(destroy_V)(v2out[i]);
  }
#if QDP_Precision == 2 || QDP_Precision == 'D'
  for(int i=0; i<nv; i++) {
    QDP_FN_destroy_V(v2inf[i]);
    QDP_FN_destroy_V(v2outf[i]);
  }
#endif

  if(wmg->profile) {
    for(int i=0; i<wmg->nlevels; i++) {
      QOP_printf0("level %i count %i tpre %.2f tcoarse %.2f tpost %.2f trest %.2f tprol %.2f\n",
	      i, wmg->mg[i].vca->count, wmg->mg[i].vca->tpre,
	      wmg->mg[i].vca->tcoarse, wmg->mg[i].vca->tpost,
	      wmg->mg[i].mgblock->trestrict, wmg->mg[i].mgblock->tprolong);
      //QOP_printf0(" restrict count %i time %.2f\n",
      //wmg->mg[i].mgblock->rcount, wmg->mg[i].mgblock->trestrict);
      //QOP_printf0("  prolong count %i time %.2f\n",
      //wmg->mg[i].mgblock->pcount, wmg->mg[i].mgblock->tprolong);
    }
  }

  res_arg->final_iter = its;
  info->final_sec = secs;
#undef NC
}

#endif // USE_MG
