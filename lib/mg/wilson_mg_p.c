#include <qop_internal.h>

#ifdef USE_MG

#include "solvers.h"

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QCDPC(x) QOP_F3_##x
#define QDPN(x) QDP_FN_##x
#else
#define QOPP(x) QOP_D_##x
#define QCDPC(x) QOP_D3_##x
#define QDPN(x) QDP_DN_##x
#endif

#define NV_BICGSTAB

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

QOP_WilsonMg *
QOP_wilsonMgNew(void)
{
  QOP_WilsonMg *QOP_malloc(wmg, QOP_WilsonMg, 1);
  wmg->wilD = NULL;
  wmg->wilF = NULL;
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
  return wmg;
}

static void
init_level(QOP_WilMgLevel l[], int n)
{
  QOP_printf0("init_level %i\n", n);
  l[n].ndim = 4;
  QOP_malloc(l[n].lattice_size, int, l[n].ndim);
  if(n==0) {
    l[n].lattice_size[0] = QDP_coord_size(0)/3;
    l[n].lattice_size[1] = QDP_coord_size(1)/3;
    l[n].lattice_size[2] = QDP_coord_size(2)/3;
    l[n].lattice_size[3] = QDP_coord_size(3)/8;
    l[n].nvecs = 24;
    l[n].npre = 0;
    l[n].npost = 5;
    l[n].scale = 1;
#ifdef SPLIT_CHIRALITIES
    l[n].fnc = 6;
#else
    l[n].fnc = 12;
#endif
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
  for(int i=0; i<l[n].ndim; i++) if(l[n].lattice_size[i]<1) l[n].lattice_size[i]=1;
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
}

static QLA_Real
norm2V(QDPN(ColorVector) **v, int nv)
{
  QDP_Lattice *lat = QDPN(get_lattice_V)(v[0]);
  QDP_Subset allf = QDP_all_L(lat);
  QLA_Real nrm2 = 0;
  for(int i=0; i<nv; i++) {
    QLA_Real tt;
    QDPN(r_eq_norm2_V)(&tt, v[i], allf);
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
  int i,t;

  t = coords[0];
  for(i=1; i<QDP_ndim(); i++) {
    t = t*QDP_coord_size(i) + coords[i];
  }
  *li = t;
}

static void
QOP_randSeed(QDP_RandomState *rs, int seed)
{
  QDP_Int *li;

  li = QDP_create_I();

  QDP_I_eq_func(li, lex_int, QDP_all);
  QDP_S_eq_seed_i_I(rs, seed, li, QDP_all);

  QDP_destroy_I(li);
}

#if 0
static void
get_nvev(QOP_F_MgArgs *mgargs, QOPP(MgOp) *op, void *opargs, QDP_Subset sub,
	 QOPP(MgOp) *eop, void *eopargs, QDP_Subset esub, int nvecs, int maxits, double cfac)
{
  int nv = mgargs->nv;
  int fnc = mgargs->fnc;
  int cnc = mgargs->cnc;
  QDP_Lattice *fine = mgargs->mgb->fine;
  QDP_Subset allf = QDP_all_L(fine);
  QDP_RandomState *rs = QDP_create_S(); // hack: lives on fine lattice
  QOP_randSeed(rs, 987654321); // needs to be fixed to use any lattice
  double dt = -QDP_time();

  QDP_FN_ColorVector *pv[nvecs][nv];
  QDP_FN_ColorVector *mpv[nv][cnc];
  for(int i=0; i<nvecs; i++) {
    for(int j=0; j<nv; j++) {
      if(i<cnc) {
	pv[i][j] = mgargs->pv[j][i];
	mpv[j][i] = pv[i][j];
      } else {
	pv[i][j] = QDP_FN_create_V_L(fnc, fine);
      }
    }
  }
  QDP_FN_ColorVector *tv[nv];
  for(int j=0; j<nv; j++) tv[j] = QDP_FN_create_V_L(fnc, fine);
  for(int j=0; j<nv; j++) QDP_FN_V_eq_gaussian_S(pv[0][j], rs, allf);
  //double nrm2 = norm2V(tv, nv);
  //QOP_printf0("nrm2 = %g\n", nrm2);

  double sv[nvecs];
  int k = 0;
  while(1) {
    V_eq_V(tv, pv[k]);

    QOP_F_eignhVN(nvecs-k, nv, &pv[k], fnc, tv, op, opargs, allf, maxits);

    for(int i=k; i<nvecs; i++) {
      QLA_Real nrm = norm2V(pv[i], nv);
      op(tv, pv[i], 1, opargs);
      QLA_Real nrm2 = norm2V(tv, nv);
      sv[i] = nrm2/nrm;
      QOP_printf0("%3i %g\n", i, sv[i]);
      QLA_Real scale = 1/sqrt(nrm);
      for(int j=0; j<nv; j++) { 
	QDP_FN_V_eq_r_times_V(pv[i][j], &scale, pv[i][j], allf);
      }
    }

    for(int i=0; i<cnc; i++) {
      double svmin = sv[i];
      int imin = i;
      for(int l=i+1; l<nvecs; l++) if(sv[l]<svmin) { svmin = sv[l]; imin = l; }
      if(i!=imin) {
	sv[i] = svmin;
	for(int j=0; j<nv; j++) { 
	  QDP_FN_V_eq_V(tv[j], pv[i][j], allf);
	  QDP_FN_V_eq_V(pv[i][j], pv[imin][j], allf);
	  QDP_FN_V_eq_V(pv[imin][j], tv[j], allf);
	}
      }
      for(int l=i+1; l<nvecs; l++) {
        QLA_Complex z = dotV(pv[i], pv[l], nv);
        for(int j=0; j<nv; j++) {
	  QDP_FN_V_meq_c_times_V(pv[l][j], &z, pv[i][j], allf);
	}
	QLA_Real nrm = norm2V(pv[l], nv);
	op(tv, pv[l], 1, opargs);
	QLA_Real nrm2 = norm2V(tv, nv);
	sv[l] = nrm2/nrm;
	//QOP_printf0("%3i %g\n", l, sv[l]);
	QLA_Real scale = 1/sqrt(nrm);
	for(int j=0; j<nv; j++) { 
	  QDP_FN_V_eq_r_times_V(pv[l][j], &scale, pv[l][j], allf);
	}
      }
      QOP_printf0("%3i %g\n", i, sv[i]);
    }
    k = 0;
    while(k<nvecs && sv[0]>=cfac*sv[k]) k++;
    if(k>=cnc) break;
  }
  for(int i=cnc; i<nvecs; i++) {
    for(int j=0; j<nv; j++) {
      QDP_FN_destroy_V(pv[i][j]);
    }
  }

  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  QOP_F_mgOrthoVn(mgargs->pv, nv, 0, cnc, allf);
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int j=0; j<nv; j++) {
    QOP_F_mgOrtho(mpv[j], cnc, mgargs->mgb);
  }
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int i=0; i<cnc; i++) {
    for(int j=0; j<nv; j++) {
      QDP_FN_V_eq_V(mgargs->rv[j][i], mgargs->pv[j][i], allf);
    }
  }
  dt += QDP_time();
  QOP_printf0("%s: %g secs\n", __func__, dt);

  for(int i=0; i<nv; i++) {
    QDP_FN_V_eq_gaussian_S(tv[i], rs, allf);
  }
  QOP_F_mgTestCoarse(mgargs, tv);
  for(int j=0; j<nv; j++) {
    QDP_FN_destroy_V(tv[j]);
  }
  QDP_destroy_S(rs);
}
#endif

static void
get_nullvecs(QOP_F_MgArgs *mgargs, QOPP(MgOp) *op, void *opargs, QLA_Real res,
	     QLA_Real change_fac, int maxit, int nvecs, int verbose,
	     int (smoother)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in,
			     int sign, void *args), void *sargs)
{
  int nv = mgargs->nv;
  int fnc = mgargs->fnc;
  int cnc = mgargs->cnc;
  QDP_Lattice *fine = mgargs->mgb->fine;
  QDP_Subset allf = QDP_all_L(fine);
  QDP_RandomState *rs = QDP_create_S(); // hack: lives on fine lattice
  QOP_randSeed(rs, 987654321); // needs to be fixed to use any lattice
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

  QDP_FN_ColorVector *tv[nv], *tv2[nv];
  for(int j=0; j<nv; j++) tv[j] = QDP_FN_create_V_L(fnc, fine);
  for(int j=0; j<nv; j++) tv2[j] = QDP_FN_create_V_L(fnc, fine);

#if 0
#ifdef NV_BICGSTAB
  QOP_Bicgstab *bcg = QOP_bicgstabInit(fine, nv, fnc);
  QOP_bicgstabSet(bcg, "verbose", verbose);
  QOP_bicgstabSet(bcg, "indent", 2);
#else
  QOP_Cgls *cg = QOP_cglsInit(fine, nv, fnc);
  QOP_cglsSet(cg, "verbose", verbose);
  QOP_cglsSet(cg, "indent", 2);
#endif
#endif

#if 1
  int tnits=0;
  int maxhits=20;
  double svhist[maxhits];
  for(int i=0; i<cnc; i++) {
    if(i==0) {
      for(int j=0; j<nv; j++) QDP_FN_V_eq_gaussian_S(pv[i][j], rs, allf);
    }
    //if(i==1) maxit /= 2;
    QLA_Real sv2=FLT_MAX, sv2o;
    int nhits=0;
    int nits=0;
    do {
      nhits++;
#if 1
      //QOP_mgOrthoVn(mpv, nv, 0, i+1, allf);
      for(int k=0; k<i; k++) {
        QLA_Complex z = dotV(pv[k], pv[i], nv);
        for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(pv[i][j], &z, pv[k][j], allf); }
      }
#else
      for(int j=0; j<nv; j++) {
        QOP_F_mgOrthoVec(mpv[j], i, mgargs->mgb, 0);
      }
#endif
#if 0
#ifdef NV_BICGSTAB
      nits += QOP_bicgstabSolveS(bcg, tv, pv[i], op, opargs, NULL, NULL, 1, res, maxit, allf);
#else
      nits += QOP_cglsSolve(cg, tv, pv[i], NULL, op, opargs, NULL, NULL, 0, res, 0, maxit, allf, allf);
      //QOP_cglsSolve(cg, tv, NULL, pv[i], op, opargs, NULL, NULL, 0, 0, res, maxit, allf, allf);
#endif
#endif
      nits += smoother(tv, pv[i], 1, sargs);
      QLA_Real nrm = norm2V(tv, nv);
      op(pv[i], tv, 1, opargs);
      QLA_Real nrm2 = norm2V(pv[i], nv);
      sv2o = sv2;
      sv2 = nrm2/nrm;
      svhist[nhits-1] = sv2;
      QLA_Real scale = 1/sqrt(nrm);
      if(i+1<cnc) for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(pv[i+1][j], pv[i][j], allf); }
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_r_times_V(pv[i][j], &scale, tv[j], allf); }
      if(verbose>0) QOP_printf0("%i %g\n", i, sv2);
      if(i>0 && nhits==4) break;
    } while(sv2<change_fac*sv2o);
    tnits += nits;
    QOP_printf0("%-3i %-3i %-5i  %-10g  %-10g", i, nhits, nits, sqrt(sv2), sqrt(sv2o));
    for(int i=nhits-3; i>0; i--) QOP_printf0("  %-10g", sqrt(svhist[i]));
    QOP_printf0("\n");
#if 0
    for(int j=0; j<nv; j++) {
      QOP_F_mgOrthoVec(mpv[j], i, mgargs->mgb, 1);
    }
#endif
  }
  QOP_printf0("setup its = %i\n", tnits);
#endif

#if 0
  int tnits=0;
  for(int i=0; i<cnc; i++) {
    int nhits=0;
    int nits=0;
    if(i==0) {
      for(int j=0; j<nv; j++) QDP_FN_V_eq_gaussian_S(pv[i][j], rs, allf);
      nhits = 10;
    } else {
      QLA_Complex z = dotV(pv[i-1], pv[i], nv);
      for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(pv[i][j], &z, pv[i-1][j], allf); }
      nhits = 1;
    }
    QLA_Real sv2=FLT_MAX, sv2o;
    {
      QLA_Real nrm = norm2V(pv[i], nv);
      op(tv2, pv[i], 1, opargs);
      QLA_Real nrm2 = norm2V(tv2, nv);
      sv2o = sv2;
      sv2 = sqrt(nrm2/nrm);
      QOP_printf0("warmup %-3i %-5i %-10g\n", 0, nits, sv2);
    }
    for(int h=0; h<nhits; h++) {
      nits += QOP_bicgstabSolveS(bcg, tv, pv[i], op, opargs, NULL, NULL, 1, res, maxit, allf);
      for(int k=0; k<i; k++) {
	QLA_Complex z = dotV(pv[k], tv, nv);
	for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(tv[j], &z, pv[k][j], allf); }
      }
      QLA_Real nrm = norm2V(tv, nv);
      op(tv2, tv, 1, opargs);
      QLA_Real nrm2 = norm2V(tv2, nv);
      sv2o = sv2;
      sv2 = sqrt(nrm2/nrm);
      QLA_Real scale = 1/sqrt(nrm);
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_r_times_V(pv[i][j], &scale, tv[j], allf); }
      QOP_printf0("warmup %-3i %-5i %-10g\n", h+1, nits, sv2);
      if(sv2>change_fac*sv2o) break;
    }
    if(i+1<cnc) for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(pv[i+1][j], pv[i][j], allf); }
    {
      // Ax
      op(tv, pv[i], 1, opargs);
      // relax
#ifdef NV_BICGSTAB
      nits += QOP_bicgstabSolveS(bcg, tv2, tv, op, opargs, NULL, NULL, 1, res, maxit/2, allf);
#else
      nits += QOP_cglsSolve(cg, tv2, tv, NULL, op, opargs, NULL, NULL, 0, res, 0, maxit, allf, allf);
      //QOP_cglsSolve(cg, tv, NULL, pv[i], op, opargs, NULL, NULL, 0, 0, res, maxit, allf, allf);
#endif
      // e
      for(int j=0; j<nv; j++) { QDP_FN_V_meq_V(tv2[j], pv[i][j], allf); }
      // ortho
      for(int k=0; k<i; k++) {
	QLA_Complex z = dotV(pv[k], tv2, nv);
	for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(tv2[j], &z, pv[k][j], allf); }
      }
      // RQ
      // norm
      QLA_Real nrm = norm2V(tv2, nv);
      op(tv, tv2, 1, opargs);
      QLA_Real nrm2 = norm2V(tv, nv);
      sv2 = sqrt(nrm2/nrm);
      QLA_Real scale = 1/sqrt(nrm);
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_r_times_V(pv[i][j], &scale, tv2[j], allf); }
      QOP_printf0("%-3i %-5i  %-10g\n", i, nits, sv2);
    }
    tnits += nits;
  }
  QOP_printf0("setup its = %i\n", tnits);
#endif

#if 0
  QLA_F_Complex evs[cnc];
  double min[cnc], ave[cnc], max[cnc];

  int tnits=0, ncycles=0, ncyclesmax=5;
  int imin=0;
  //double tol=1e-4;
  double tol=1e-1;
  for(int j=0; j<nv; j++) QDP_FN_V_eq_gaussian_S(pv[imin][j], rs, allf);
  do {
    int nits=0;
    ncycles++;
    for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(tv[j], pv[imin][j], allf); }
    //for(int i=cnc-1; i>=imin; i--) {
    for(int i=imin; i<cnc; i++) {
      if(ncycles!=1) {
	for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(tv[j], pv[i][j], allf); }
      }
      nits += QOP_bicgstabSolveS(bcg, pv[i], tv, op, opargs, NULL, NULL, 1, res, maxit, allf);
      //for(int k=0; k<cnc; k++) {
      //if(k>=imin && k<=i) continue;
      for(int k=0; k<i; k++) {
	QLA_Complex z = dotV(pv[k], pv[i], nv);
	for(int j=0; j<nv; j++) { QDP_FN_V_meq_c_times_V(pv[i][j], &z, pv[k][j], allf); }
      }
      QLA_Real nrm = norm2V(pv[i], nv);
      op(tv, pv[i], 1, opargs);
      QLA_Real nrm2 = norm2V(tv, nv);
      QLA_Real scale = 1/sqrt(nrm);
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_r_times_V(pv[i][j], &scale, pv[i][j], allf); }
      for(int j=0; j<nv; j++) { QDP_FN_V_eq_V(tv[j], pv[i][j], allf); }
      QOP_printf0("%-3i  %-10g\n", i, sqrt(nrm2/nrm));
    }

    //QOP_F_rRitzHarm(pv[imin], cnc-imin, nv, op, opargs, evs+imin, QOP_SORT_ASCEND, QOP_NORM_2, allf);
    QOP_F_rRitzHarm(pv[0], cnc, nv, op, opargs, evs, QOP_SORT_ASCEND, QOP_NORM_2, allf);
    int imin0 = imin;
    //double rstop = tol*fabs(QLA_real(evs[0]));
    double rstop = tol*QLA_norm2_c(evs[0]);
    for(int i=0; i<cnc; i++) {
      double err;
      //QOP_printf0("ev%i = %-12g %-12g  %-12g\n", i, QLA_real(evs[i]), QLA_imag(evs[i]), sqrt(QLA_norm2_c(evs[i])));
      {
	QLA_Real nrm = norm2V(pv[i], nv);
	op(tv, pv[i], 1, opargs);
	QLA_Complex z = dotV(pv[i], tv, nv);
	QLA_c_eq_r_times_c(z, 1/nrm, z);
	//QOP_printf0("ev%i = %-12g %-12g  %-12g\n", i, QLA_real(z), QLA_imag(z), sqrt(QLA_norm2_c(z)));
	QLA_c_meq_c(z, evs[i]);
	err = sqrt(QLA_norm2_c(z));
	QOP_printf0("ev%i = %-12g %-12g  %-12g  %-12g\n", i, QLA_real(evs[i]), QLA_imag(evs[i]), sqrt(QLA_norm2_c(evs[i])), err);
      }
      //if(imin0<i && fabs(QLA_real(evs[i]))>rstop) imin0 = i;
      //if(imin0<i && QLA_norm2_c(evs[i])>rstop) imin0 = i;
      if(imin0==i-1 && err<rstop) imin0 = i;
    }
    if(imin0>cnc0) imin0=cnc0;
#if 0
    for(int j=0; j<nv; j++) {
      QOP_F_mgOrthoSort(mpv[j], imin, cnc, mgargs->mgb, min, ave, max);
      for(int i=0; i<cnc; i++) {
	QOP_printf0("min,ave,max %i %-2i = %-12g %-12g %-12g\n", j, i, min[i], ave[i], max[i]);
	if(imin0<i && max[i]>tol*max[0]) imin0 = i;
      }
    }
#endif
    imin = imin0;
    tnits += nits;
    QOP_printf0("cycle %i  its %i  imin = %i\n", ncycles, nits, imin);
  } while(ncycles<ncyclesmax);
  //} while(imin<cnc0-1 && ncycles<ncyclesmax);
  QOP_printf0("setup its = %i\n", tnits);
#endif

#if 0
  QOP_F_rRitzHarm(*pv, cnc, nv, op, opargs, evs, QOP_SORT_ASCEND, QOP_NORM_E, allf);
  for(int i=0; i<cnc; i++) {
    QLA_Real nrm = norm2V(pv[i], nv);
    op(tv, pv[i], 1, opargs);
    QLA_Complex z = dotV(pv[i], tv, nv);
    QLA_c_eq_r_times_c(z, 1/nrm, z);
    //QOP_printf0("ev%i = %-12g %-12g  %-12g\n", i, QLA_real(z), QLA_imag(z), sqrt(QLA_norm2_c(z)));
    QLA_c_meq_c(z, evs[i]);
    double err = sqrt(QLA_norm2_c(z));
    QOP_printf0("ev%i = %-12g %-12g  %-12g  %-12g\n", i, QLA_real(evs[i]), QLA_imag(evs[i]), sqrt(QLA_norm2_c(evs[i])), err);
  }
  for(int j=0; j<nv; j++) {
    QOP_F_mgOrthoSort(mpv[j], 0, cnc, mgargs->mgb, min, ave, max);
    for(int i=0; i<cnc; i++) {
      QOP_printf0("min,ave,max %i %-2i = %-12g %-12g %-12g\n", j, i, min[i], ave[i], max[i]);
    }
  }
#else
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  QOP_F_mgOrthoVn(mgargs->pv, nv, 0, cnc0, allf);
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int j=0; j<nv; j++) {
    QOP_F_mgOrtho(mpv[j], cnc, mgargs->mgb);
  }
#endif
  print_norms(*mpv, nv, fnc, cnc, op, opargs, 1);
  for(int i=0; i<cnc0; i++) {
    for(int j=0; j<nv; j++) {
#ifndef SPLIT_CHIRALITIES
      if(fnc==12) { // hack -- top level without splitting chiralities
	QDP_F3_D_eq_gamma_times_D((QDP_F3_DiracFermion*)mgargs->rv[j][i], (QDP_F3_DiracFermion*)mgargs->pv[j][i], 15, allf);
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
#if 0
#ifdef NV_BICGSTAB
  QOP_bicgstabFree(bcg);
#else
  QOP_cglsFree(cg);
#endif
#endif
}

static void
create_level(QOP_WilsonMg *wmg, int n)
{
  QOP_printf0("creating MG level %i\n", n);

  QOP_WilMgLevel *l = wmg->mg;
  QDP_Lattice *lat0;
  int nv = l[n].nv;

  wmg->nvwaF.kappa = wmg->kappanv;
  wmg->vcwaF.kappa = wmg->kappa;

  void (*smoothop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  if(n==0) {
    lat0 = QDP_get_default_lattice();
#ifdef SPLIT_CHIRALITIES
#ifdef EO_PREC
    l[n].vcop = QOP_F3_wilEoV2;
    l[n].nvop = QOP_F3_wilPV2;
    l[n].cfop = QOP_F3_wilPV2;
#else
    l[n].vcop = QOP_F3_wilDV2;
    l[n].nvop = QOP_F3_wilDV2;
    l[n].cfop = QOP_F3_wilDV2;
#endif
    smoothop = QOP_F3_wilEoV2;
#else
#ifdef EO_PREC
    l[n].vcop = QOP_F3_wilEoV1;
    l[n].nvop = QOP_F3_wilPV1;
    l[n].cfop = QOP_F3_wilPV1;
#else
    l[n].vcop = QOP_F3_wilDV1;
    l[n].nvop = QOP_F3_wilDV1;
    l[n].cfop = QOP_F3_wilDV1;
#endif
    smoothop = QOP_F3_wilEoV1;
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
  QOP_printf0("creating lattice:");
  for(int i=0; i<l[n].ndim; i++) QOP_printf0(" %i", l[n].lattice_size[i]);
  QOP_printf0("\n");
  l[n].lattice = QDP_create_lattice(QDP_get_default_layout(), NULL,
				    l[n].ndim, l[n].lattice_size);

  QOP_printf0("\ncreating block\n");
  l[n].mgblock = QOP_mgCreateBlockFromLattice(lat0, l[n].lattice);

  QOP_printf0("creating args\n");
  l[n].mgargs = QOP_mgCreateArgs(l[n].mgblock, l[n].nvecs, l[n].fnc, l[n].nv, NULL, NULL);

  // create some vectors for the vcycle and also used in setup
  QDP_FN_ColorVector **QOP_malloc(fr, QDP_FN_ColorVector *, nv);
  QDP_FN_ColorVector **QOP_malloc(fp, QDP_FN_ColorVector *, nv);
  QDP_FN_ColorVector **QOP_malloc(fAp, QDP_FN_ColorVector *, nv);
  for(int i=0; i<nv; i++) {
    fr[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
    fp[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
    fAp[i] = QDP_FN_create_V_L(l[n].fnc, lat0);
  }

  // get nullvecs
  QOP_printf0("creating null vecs\n");
  {
    QOP_Bicgstab *bcg = QOP_bicgstabInit(lat0, l[n].nv, l[n].fnc);
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
      nvsa->project = QOP_F3_wilEoProjectV2;
      nvsa->reconstruct = QOP_F3_wilEoReconstructPV2;
#else
      nvsa->project = QOP_F3_wilEoProjectV1;
      nvsa->reconstruct = QOP_F3_wilEoReconstructPV1;
#endif
#else
#ifdef SPLIT_CHIRALITIES
      nvsa->project = QOP_F3_wilEoProjectV2;
      nvsa->reconstruct = QOP_F3_wilEoReconstructV2;
#else
      nvsa->project = QOP_F3_wilEoProjectV1;
      nvsa->reconstruct = QOP_F3_wilEoReconstructV1;
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
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose);
    //#ifdef EO_PREC
    //get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveEo, nvsa);
    get_nullvecs(l[n].mgargs, l[n].nvop, l[n].nvopargs, l[n].setup_res, l[n].setup_change_fac, l[n].setup_maxit, l[n].setup_nvecs, l[n].verbose, QOP_F_bicgstabSolveEo, nvsa);
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

  // create coarse Dslash
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

  // create true coarse op
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
  QOP_F_mgCloneOp(QOP_F_mgC2fOp, cfoa, dargs);
#endif
  l[n].dargs = dargs;

  // coarse solver
  QOP_F_Gcr *gcrc = QOP_F_gcrInit(l[n].lattice, nv, l[n].nvecs, l[n].ngcr);
  QOP_F_gcrSet(gcrc, "verbose", l[n].verbose);
  QOP_F_gcrSet(gcrc, "indent", 2*(n+2));
  l[n].gcrc = gcrc;

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
  sa->opargs = dargs;
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

  // fine to coarse solver
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
  mgoa2->opargs = sa;
  mgoa2->cin = cin;
  mgoa2->cout = cout;

  // smoother and vcycle
#define MAX(a,b) ((a)<(b)?(b):(a));
  int ngcrf = MAX(l[n].npre, l[n].npost);
  QOP_F_Gcr *gcrf = QOP_F_gcrInit(lat0, nv, l[n].fnc, ngcrf);
  QOP_F_gcrSet(gcrf, "verbose", l[n].verbose);
  QOP_F_gcrSet(gcrf, "indent", 2*(n+1));
  l[n].gcrf = gcrf;
  QOP_F_MgVcycleArgs *QOP_malloc(vc, QOP_F_MgVcycleArgs, 1);
  l[n].vca = vc;
  vc->op = l[n].vcop;
  vc->opargs = l[n].vcopargs;
  vc->cop = QOP_F_mgF2cOp;
  vc->copargs = mgoa2;
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

  if(n>0) {
    // need to set vc as preconditioner for prev solver    
    l[n-1].sa->pop = QOP_F_mgVcycle;
    l[n-1].sa->popargs = vc;
  }

}

static void
free_level(QOP_WilMgLevel *l)
{
  QOP_printf0("*** freeing level %p ***\n", l);
  QOP_free(l->lattice_size);
  // FIXME
}

static void
setNumLevels(QOP_WilsonMg *wmg, int nlevels)
{
  //int ret = wmg->nlevels;
  if(nlevels>0) {
    for(int i=nlevels; i<wmg->nlevels; i++) {
      free_level(&wmg->mg[i]);
    }
    wmg->mg = realloc(wmg->mg, nlevels*sizeof(QOP_WilMgLevel));
    for(int i=wmg->nlevels; i<nlevels; i++) {
      init_level(wmg->mg, i);
    }
    wmg->nlevels = nlevels;
  }
  //return ret;
}

void
QOP_wilsonMgFree(QOP_WilsonMg *wmg)
{
  for(int i=0; i<wmg->nlevels; i++) {
    free_level(&wmg->mg[i]);
  }
  if(wmg->mg) {
    QOP_free(wmg->mg);
  }
  //QOP_gcrFree(gcro);
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

#if 0
void
QOP_wilsonMgSetVecs(QOP_WilsonMg *wmg, int nv, QDP_DiracFermion *vv[nv])
{
  QOP_malloc(wmg->v, QDP_DiracFermion *, nv);
  wmg->nv = nv;
  for(int i=0; i<nv; i++) {
    wmg->v[i] = QDP_create_D();
    QDP_D_eq_D(wmg->v[i], vv[i], QDP_all);
  }
  //mgorth(wmg, vv);
}
#endif

void
QOP_wilsonMgSetLinks(QOP_WilsonMg *wmg, QOP_FermionLinksWilson *wil)
{
  wmg->wilF = wil;
  wmg->vcwaF.wil = wil;
  wmg->nvwaF.wil = wil;
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

#endif // single precision

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
  QDP_Subset ssub = QDP_even;
#else
  QDP_Subset ssub = QDP_all;
#endif

  QDP_DiracFermion *r = QDP_create_D();
  QDP_DiracFermion *out2 = QDP_create_D();
  QDP_DiracFermion *eovec = QDP_create_D();
  QDPN(ColorVector) *v2in[nv], *v2out[nv];
  for(int i=0; i<nv; i++) {
    v2in[i] = QDPN(create_V)(fnc);
    v2out[i] = QDPN(create_V)(fnc);
  }

  QOPP(MgOp) *pop;
  void *popargs;
#if QDP_Precision == 1 || QDP_Precision == 'F'
  pop = QOP_F_mgVcycle;
  popargs = wmg->mg[0].vca;
#else
  QDP_FN_ColorVector *v2inf[nv], *v2outf[nv];
  for(int i=0; i<nv; i++) {
    v2inf[i] = QDP_FN_create_V(fnc);
    v2outf[i] = QDP_FN_create_V(fnc);
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
  QDP_r_eq_norm2_D(&norm2in, in, QDP_all);
  rsqstop = res_arg->rsqmin*norm2in;
  double gcrres = -sqrt(res_arg->rsqmin*norm2in);

  double secs=0;
  int its = 0;
  while(1) {
    QOP_wilsonDslash(r, out, flw, kappa, 1, QOP_EVENODD, QOP_EVENODD);
    QDP_D_eq_D_minus_D(r, in, r, QDP_all);
    QDP_r_eq_norm2_D(&rsq, r, QDP_all);
    if(its>0) QOP_printf0("iters = %4i  secs = %8f  rsq = %g\n", its, secs, rsq/norm2in);
    if(rsq<rsqstop || its>=itmax) break;

#ifdef EO_PREC
    QOP_wilEoProjectD(eovec, r, &wa);
#else
    QDP_D_eq_D(eovec, r, QDP_all);
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
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilEoV2, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_even);
#else
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilEoV1, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_even);
#endif
#else
#ifdef SPLIT_CHIRALITIES
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilDV2, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_all);
#else
    its += QOP_gcrSolve(gcro, v2out, v2in, QOP_wilDV1, &wa, pop, popargs, gcrres, itmax, ngcr, QDP_all);
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
    QDP_D_eq_D(out2, eovec, QDP_all);
#endif
    QDP_D_peq_D(out, out2, QDP_all);
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
}

#endif // USE_MG
