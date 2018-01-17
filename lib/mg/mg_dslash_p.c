#include <qop_internal.h>

#include <math.h>
#include <limits.h>
#include <qla_fn.h>
#include <qla_dn.h>
#include <qdp_fn.h>
#include <qdp_dn.h>
#include <qdp_dfn.h>
#include "solvers.h"

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QDPN(x) QDP_FN_##x
#define QLAN(x) QLA_FN_##x
#else
#define QOPP(x) QOP_D_##x
#define QDPN(x) QDP_DN_##x
#define QLAN(x) QLA_DN_##x
#endif

#define printf0 QOP_printf0

static int
addpath(int k, int nd, int *path, int *ls, int *pidx)
{
  int x[nd];
  for(int i=0; i<nd; i++) {
    int j = i;
    if(pidx) j = pidx[i];
    x[j] = (k - path[j] + abs(path[j])*ls[j]) % ls[j];
    k /= ls[j];
  }
  k = 0;
  for(int i=nd-1; i>=0; i--) {
    int j = i;
    if(pidx) j = pidx[i];
    k = k*ls[j] + x[j];
  }
  return k;
}

static int
checkfit(int *used, int k, int ns, int npaths, int nd, int *paths, int *ls, int *pidx)
{
  int fits = 1;
  for(int ip=0; ip<npaths; ip++) {
    int j = addpath(k, nd, &paths[nd*ip], ls, pidx);
    if(used[j]==ns) {
      for(int i=0; i<ip; i++) {
	int jj = addpath(k, nd, &paths[nd*i], ls, pidx);
	used[jj] = -1;
      }
      fits = 0;
      break;
    }
    used[j] = ns;
  }
  return fits;
}

QOPP(MgDslashArgs) *
QOPP(mgCreateDslash)(int nvout, int nvin, int nc, int npaths,
		     int *paths, QDP_Lattice *lat)
{
  QOPP(MgDslashArgs) *QOP_malloc(da, QOPP(MgDslashArgs), 1);
  da->nvin = nvin;
  da->nvout = nvout;
  da->nc = nc;
  da->lat = lat;
  da->all = QDP_all_L(lat);
  da->even = QDP_even_L(lat);
  da->odd = QDP_odd_L(lat);
  int nd = QDP_ndim_L(lat);
  int ls[nd];
  QDP_latsize_L(lat, ls);

  int dests[npaths];
  for(int i=0; i<npaths; i++) {
    int k = addpath(0, nd, &paths[nd*i], ls, NULL);
    dests[i] = k;
  }
  int havelocal = 0;
  int ns = 0;
  int pathidx[npaths];
  for(int i=0; i<npaths; i++) {
    if(dests[i]==0) { havelocal = 1; continue; }
    pathidx[ns] = i;
    ns++;
    for(int j=0; j<i; j++) {
      if(dests[j]==dests[i]) { ns--; break; }
    }
  }
  da->nshifts = ns;
  int np = ns;
  if(havelocal) np++;

  da->npaths = np;
  QOP_malloc(da->paths, int, np*nd);
  if(ns) {
    QOP_malloc(da->shifts, QDP_Shift, ns);
    QOP_malloc(da->fb, QDP_ShiftDir, ns);
    //QOP_malloc(da->shiftsadj, QDP_Shift, ns);
    //QOP_malloc(da->fbadj, QDP_ShiftDir, ns);
    for(int i=0; i<ns; i++) {
      printf0(" shift %2i :", i);
      for(int j=0; j<nd; j++) {
	da->paths[nd*i+j] = paths[nd*pathidx[i]+j];
	printf0(" %2i", da->paths[nd*i+j]);
      }
      printf0("\n");
      da->shifts[i] = QDP_create_shift_L(lat, &paths[nd*pathidx[i]]);
      da->fb[i] = QDP_forward;
      //da->shiftsadj[i] = da->shifts[i]
      //da->fbadj[ns] = QDP_backward;
    }
  } else {
    da->shifts = NULL;
    da->fb = NULL;
  }
  if(havelocal) {
    int i = ns;
    printf0(" shift %2i :", i);
    for(int j=0; j<nd; j++) {
      da->paths[nd*i+j] = 0;
      printf0(" %2i", da->paths[nd*i+j]);
    }
    printf0("\n");
  }

  int nl = 2*ns + 2;
  QOP_malloc(da->temp, QDPN(ColorVector) **, nl);
  for(int i=0; i<nl; i++) {
    QOP_malloc(da->temp[i], QDPN(ColorVector) *, nvin);
    for(int j=0; j<nvin; j++) {
      da->temp[i][j] = QDPN(create_V_L)(nc, lat);
    }
  }
  QOP_malloc(da->links, QDPN(ColorMatrix) ***, nl);
  for(int i=0; i<nl; i++) {
    QOP_malloc(da->links[i], QDPN(ColorMatrix) **, nvout);
    for(int j=0; j<nvout; j++) {
      QOP_malloc(da->links[i][j], QDPN(ColorMatrix) *, nvin);
      for(int k=0; k<nvin; k++) {
	da->links[i][j][k] = QDPN(create_M_L)(nc, lat);
      }
    }
  }
  return da;
}

void
QOPP(mgFreeDslash)(QOPP(MgDslashArgs) *da)
{
  int nl = 2*da->nshifts + 2;
  int nvout = da->nvout;
  int nvin = da->nvin;
  for(int i=0; i<nl; i++) {
    for(int j=0; j<nvin; j++) {
      QDPN(destroy_V)(da->temp[i][j]);
    }
    QOP_free(da->temp[i]);
  }
  QOP_free(da->temp);
  for(int i=0; i<nl; i++) {
    for(int j=0; j<nvout; j++) {
      for(int k=0; k<nvin; k++) {
	QDPN(destroy_M)(da->links[i][j][k]);
      }
      QOP_free(da->links[i][j]);
    }
    QOP_free(da->links[i]);
  }
  QOP_free(da->links);
  if(da->nshifts) {
    for (int i = 0 ; i < da->nshifts ; i++)
        QDP_destroy_shift(da->shifts[i]);
    QOP_free(da->shifts);
  }
  if (da->fb)
    QOP_free(da->fb);
  QOP_free(da->paths);
  QOP_free(da);
}

static int
getfac(int n)
{
  if(n==1) return 1;
  int f = 2;
  while(n%f) f++;
  return f;
}

static void
getper(int *period, int n, int *latsize, int nd)
{
  int r[nd];
  for(int i=0; i<nd; i++) r[i] = latsize[i];
  for(int i=0; i<nd; i++) period[i] = 1;
  int nf = 1;
  int k = 0;
  while(n>0) {
    if(n&1) {
      for(int i=0; i<nf; i++) {
	int f = getfac(r[k]);
	if(f==1) { period[0] = 0; return; }
	r[k] /= f;
	period[k] *= f;
      }
    }
    n >>= 1;
    if(++k>=nd) { k = 0; nf *= 2; }
  }
}

static void
getperm(int *pidx, int *period, int p, int nd)
{
  for(int i=0; i<nd; i++) pidx[i] = i;
  for(int i=0; i<nd; i++) {
    int k = i + (p %(nd-i));
    int t = pidx[k];
    for(int j=k-1; j>=i; j--) {
      if(pidx[j]==t) { pidx[0] = 0; return; }
      pidx[j+1] = pidx[j];
    }
    pidx[i] = t;
    p /= nd-i;
  }
}

static int
getns(int *sp, int *used, int *paths, int npaths, int *per, int *pidx, int nd)
{
  int vol = 1;
  for(int i=0; i<nd; i++) vol *= per[i];
  for(int i=0; i<vol; i++) { sp[i] = -1; used[i] = -1; }
  int ns = 0;
  while(1) {
    int n = 0;
    for(int k=0; k<vol; k++) {
      if(sp[k]<0) {
	if(checkfit(used, k, ns, npaths, nd, paths, per, pidx)) {
	  sp[k] = ns;
	  n++;
	}
      }
    }
    if(n==0) break;
    ns++;
  }
  return ns;
}

static void
get_srcpat(int **srcpat, int *ns, int *period, int *paths, int npaths,
	   int *latsize, int nd)
{
  int vol = 1;
  for(int i=0; i<nd; i++) vol *= latsize[i];
  int *sp = malloc(vol*sizeof(int));
  int *used = malloc(vol*sizeof(int));
  int nmax = 255; // assume 4d
  int nperm = 24; // ditto
  int bestn = 0;
  int bestp = 0;
  int bestns = INT_MAX;
  int pidx[nd];
  for(int n=0; n<=nmax; n++) {
    getper(period, n, latsize, nd);
    if(period[0]>0) {
      for(int p=0; p<nperm; p++) {
	getperm(pidx, period, p, nd);
	if(pidx[0]>0) {
	  int s = getns(sp, used, paths, npaths, period, pidx, nd);
	  //printf0("%i %i %i %i %i %i %i\n", n, p, s, pidx[0],pidx[1],pidx[2],pidx[3]);
	  if(s>0) {
	    if(s<bestns) {
	      bestns = s;
	      bestn = n;
	      bestp = p;
	    }
	  }
	}
      }
    }
  }

  getper(period, bestn, latsize, nd);
  getperm(pidx, period, bestp, nd);
  int pvol = 1;
  for(int i=0; i<nd; i++) pvol *= period[i];
  QOP_malloc(*srcpat, int, pvol);
  getns(sp, used, paths, npaths, period, pidx, nd);
  //printf0("%i %i %i\n", bestn, bestp, bestns);
  //printf0("%i %i %i %i\n", period[0], period[1], period[2], period[3]);
  //printf0("%i %i %i %i\n", pidx[0], pidx[1], pidx[2], pidx[3]);
  for(int i=0; i<pvol; i++) {
    int k = i;
    int x[nd];
    for(int j=0; j<nd; j++) {
      x[pidx[j]] = k % period[pidx[j]];
      k /= period[pidx[j]];
    }
    k = 0;
    for(int j=nd-1; j>=0; j--) {
      k = k*period[j] + x[j];
    }
    (*srcpat)[k] = sp[i];
    //printf(" %i -> %i\n", i, k);
  }
  *ns = bestns;
  free(sp);
  free(used);
}

typedef struct {
  int *table;
  int *path;
  int *ls;
  int nd;
} patfunc_args_t;

static int
patfunc(QDP_Lattice *lat, int x[], void *args)
{
  patfunc_args_t *pfa = (patfunc_args_t *)args;
  int k = 0;
  for(int i=pfa->nd-1; i>=0; i--) {
    k = k*pfa->ls[i] + (x[i] + pfa->path[i] + abs(pfa->path[i])*pfa->ls[i])%pfa->ls[i];
  }
  return pfa->table[k];
}

void
QOPP(mgCloneOp)(void op(QDPN(ColorVector) *Ax[],
			QDPN(ColorVector) *x[], int sign, void *args),
		void *args, QOPP(MgDslashArgs) *da)
{
  int nvout = da->nvout;
  int nvin = da->nvin;
  int nc = da->nc;
  int npaths = da->npaths;
  int *paths = da->paths;
  QDP_Lattice *lat = da->lat;
  QDP_Subset all = da->all;
  double dtime = -QDP_time();
  int nd = QDP_ndim_L(lat);
  int latsize[nd];
  QDP_latsize_L(lat, latsize);

  // get pattern
  int period[nd];
  int *srcpat;
  int ns;
  get_srcpat(&srcpat, &ns, period, paths, npaths, latsize, nd);
  if(0) {
    int pvol = 1;
    for(int i=0; i<nd; i++) pvol *= period[i];
    printf0("ns = %i\n", ns);
    for(int i=0; i<pvol; i++) {
      printf0(" %i %i\n", i, srcpat[i]);
    }
  }

  patfunc_args_t pfa;
  pfa.table = srcpat;
  pfa.ls = period;
  int pfapath[nd];
  pfa.path = pfapath;
  for(int i=0; i<nd; i++) pfa.path[i] = 0;
  pfa.nd = nd;
  QDP_Subset *sub = QDP_create_subset_L(lat, patfunc, &pfa, sizeof(pfa), ns);
  QDP_Subset *nsub[npaths];
  for(int i=0; i<npaths; i++) {
    pfa.path = &paths[nd*i];
    nsub[i] = QDP_create_subset_L(lat, patfunc, &pfa, sizeof(pfa), ns);
  }

  QDPN(ColorVector) *in[nvin], *out[nvout];
  for(int iv=0; iv<nvin; iv++) in[iv] = QDPN(create_V_L)(nc, lat);
  for(int iv=0; iv<nvout; iv++) out[iv] = QDPN(create_V_L)(nc, lat);

  printf0("%i:", ns);
  for(int s=0; s<ns; s++) {
    if(s%5==0) printf0(" %i", s);
    for(int iv=0; iv<nvin; iv++) {
      for(int ic=0; ic<nc; ic++) {
	//printf0("%i %i %i\n", s, iv, ic);
	QLAN(ColorVector)(nc, cv);
	QLAN(V_eq_zero)(nc, &cv);
	QLA_c_eq_r(QLA_elem_V(cv,ic), 1);
	for(int jv=0; jv<nvin; jv++) {
	  QDPN(V_eq_zero)(in[jv], all);
	}
	QDPN(V_eq_v)(in[iv], &cv, sub[s]);
	op(out, in, 1, args);
	for(int i=0; i<npaths; i++) {
	  for(int jv=0; jv<nvout; jv++) {
	    QDPN(M_eq_colorvec_V)(da->links[i][jv][iv], out[jv], ic, nsub[i][s]);
	  }
	}
      }
    }
  }
  printf0("\n");

  QDP_destroy_subset(sub);
  for(int i=0; i<npaths; i++) {
    QDP_destroy_subset(nsub[i]);
  }
  QOP_free(srcpat);
  for(int iv=0; iv<nvin; iv++) QDPN(destroy_V)(in[iv]);
  for(int iv=0; iv<nvout; iv++) QDPN(destroy_V)(out[iv]);

  if(nvin==nvout) QOPP(mgSetDiaginv)(da);
  QOPP(mgTestDslash)(op, args, da);

  dtime += QDP_time();
  printf0("%s: %g secs\n", __func__, dtime);
}

static void
links2zmat(int n, QLAN(ColorMatrix)(n, (*m)),
	   int nc, int nv, QLAN(ColorMatrix)(nc, (*q[nv][nv])), int s)
{
  for(int i=0; i<n; i++) {
    int ic = i%nc;
    int iv = i/nc;
    for(int j=0; j<n; j++) {
      int jc = j%nc;
      int jv = j/nc;
      QLA_c_eq_c(QLA_elem_M(*m,i,j), QLA_elem_M(q[iv][jv][s],ic,jc));
    }
  }
}

static void
zmat2links(int nc, int nv, QLAN(ColorMatrix)(nc, (*q[nv][nv])), int s,
	   int n, QLAN(ColorMatrix)(n, (*m)))
{
  for(int i=0; i<n; i++) {
    int ic = i%nc;
    int iv = i/nc;
    for(int j=0; j<n; j++) {
      int jc = j%nc;
      int jv = j/nc;
      QLA_c_eq_c(QLA_elem_M(q[iv][jv][s], ic, jc), QLA_elem_M(*m,i,j));
    }
  }
}

void
QOPP(mgSetDiaginv)(QOPP(MgDslashArgs) *da)
{
  int nv = da->nvin; // assume nvin==nvout
  int nc = da->nc;
  int nl = da->nshifts;
  QDP_Lattice *lat = da->lat;
  QLAN(ColorMatrix)(nc, (*qs[nv][nv]));
  QLAN(ColorMatrix)(nc, (*qd[nv][nv]));

  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      qs[i][j] = QDPN(expose_M)(da->links[nl][i][j]);
      qd[i][j] = QDPN(expose_M)(da->links[nl+1][i][j]);
    }
  }

  int ns = QDP_sites_on_node_L(lat);
  int n = nv*nc;

#pragma omp parallel for
  for(int s=0; s<ns; s++) {
    QLAN(ColorMatrix)(n, m1);
    QLAN(ColorMatrix)(n, m2);
    links2zmat(n, &m1, nc, nv, qs, s);
    QLAN(M_eq_inverse_M)(n, &m2, &m1);
    zmat2links(nc, nv, qd, s, n, &m2);
  }

  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      QDPN(reset_M)(da->links[nl][i][j]);
      QDPN(reset_M)(da->links[nl+1][i][j]);
    }
  }

  QDPN(ColorMatrix) *qcm[nv][nv];
  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      qcm[i][j] = QDPN(create_M_L)(nc, lat);
    }
  }
  QLAN(ColorMatrix)(nc, (*qs2[nv][nv]));
  for(int l=0; l<nl; l++) {
    for(int i=0; i<nv; i++) {
      for(int j=0; j<nv; j++) {
	QDPN(M_eq_sM)(qcm[i][j], da->links[nl+1][i][j], da->shifts[l], da->fb[l], da->all);
      }
    }

    for(int i=0; i<nv; i++) {
      for(int j=0; j<nv; j++) {
	qs[i][j] = QDPN(expose_M)(da->links[l][i][j]);
	qs2[i][j] = QDPN(expose_M)(qcm[i][j]);
	qd[i][j] = QDPN(expose_M)(da->links[nl+2+l][i][j]);
      }
    }

#pragma omp parallel for
    for(int s=0; s<ns; s++) {
      QLAN(ColorMatrix)(n, m1);
      QLAN(ColorMatrix)(n, m2);
      QLAN(ColorMatrix)(n, m3);
      links2zmat(n, &m1, nc, nv, qs, s);
      links2zmat(n, &m2, nc, nv, qs2, s);
      QLAN(M_eq_M_times_M)(n, &m3, &m1, &m2);
      zmat2links(nc, nv, qd, s, n, &m3);
    }

    for(int i=0; i<nv; i++) {
      for(int j=0; j<nv; j++) {
	QDPN(reset_M)(da->links[l][i][j]);
	QDPN(reset_M)(qcm[i][j]);
	QDPN(reset_M)(da->links[nl+2+l][i][j]);
      }
    }
  }
  for(int i=0; i<nv; i++) {
    for(int j=0; j<nv; j++) {
      QDPN(destroy_M)(qcm[i][j]);
    }
  }
}

static void
vpfunc(int nc, QLAN(ColorVector)(nc, (*dest)), int index)
{
  if(index==0 && QDP_this_node==0) {
    QLA_c_eq_r(QLA_elem_V(*dest, 0), 1);
  }
}

static void
vcomp(QDPN(ColorVector) *v1[], QDPN(ColorVector) *v2[], int nv, int nc)
{
  QDP_Lattice *lat = QDPN(get_lattice_V)(v1[0]);
  QDPN(ColorVector) *d = QDPN(create_V_L)(nc, lat);
  QDPN(ColorVector) *s = QDPN(create_V_L)(nc, lat);
  for(int i=0; i<nv; i++) {
    QLA_Real dn, sn, v1n, v2n;
    QDPN(V_eq_V_minus_V)(d, v1[i], v2[i], QDP_all_L(lat));
    QDPN(V_eq_V_plus_V)(s, v1[i], v2[i], QDP_all_L(lat));
    QDPN(r_eq_norm2_V)(&dn, d, QDP_all_L(lat)); 
    QDPN(r_eq_norm2_V)(&sn, s, QDP_all_L(lat));
    QDPN(r_eq_norm2_V)(&v1n, v1[i], QDP_all_L(lat)); 
    QDPN(r_eq_norm2_V)(&v2n, v2[i], QDP_all_L(lat));
    printf0("vcomp %i  ( %g %g )  %g / %g = %g\n", i, v1n, v2n, dn, sn, dn/sn);
  }
  QDPN(destroy_V)(d);
  QDPN(destroy_V)(s);
}

static void
print_nonzeros(int nc, QLAN(ColorVector)(nc, (*dest)), int index)
{
  double eps = 1e-10;
  for(int i=0; i<nc; i++) {
    QLA_Real t = QLA_norm2_c(QLA_elem_V(*dest,i));
    if(t>eps) printf("%i %i %i %g\n", QDP_this_node, index, i, sqrt(t));
  }
}

void
QOPP(mgTestDslash)(void op(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[],
			   int sign, void *args),
		   void *args, QOPP(MgDslashArgs) *da)
{
  int nvout = da->nvout;
  int nvin = da->nvin;
  int nc = da->nc;
  QDP_Lattice *lat = da->lat;
  QDP_Subset all = QDP_all_L(lat);
  QDPN(ColorVector) *in[nvin], *out1[nvout], *out2[nvout], *out3[nvout];
  for(int i=0; i<nvin; i++) {
    in[i] = QDPN(create_V_L)(nc, lat);
    QDPN(V_eq_zero)(in[i], all);
  }
  for(int i=0; i<nvout; i++) {
    out1[i] = QDPN(create_V_L)(nc, lat);
    out2[i] = QDPN(create_V_L)(nc, lat);
    out3[i] = QDPN(create_V_L)(nc, lat);
  }

  QDPN(V_eq_funcit)(in[0], vpfunc, all);

  op(out1, in, 1, args);
  QOPP(mgDslash)(out2, in, 1, da, QOP_EVENODD, QOP_EVENODD);
  vcomp(out1, out2, nvout, nc);
  for(int i=0; i<nvout; i++) {
    //QDPN(V_eq_funci)(out1[i], print_nonzeros, all);
    //QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
    QDPN(V_meq_V)(out2[i], out1[i], all);
    QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
  }

  if(nvin==nvout) {
    QOPP(mgDiag)(out2, in, 1, da, QOP_EVEN);
    QOPP(mgDiaginv)(out1, out2, 1, da, QOP_EVEN);
    QOPP(mgDiag)(out2, in, 1, da, QOP_ODD);
    QOPP(mgDiaginv)(out1, out2, 1, da, QOP_ODD);
    vcomp(in, out1, nvin, nc);
    for(int i=0; i<nvin; i++) {
      QDPN(V_meq_V)(out1[i], in[i], all);
      QDPN(V_eq_funci)(out1[i], print_nonzeros, all);
    }
  }

#if 0
  for(int i=0; i<nv; i++) { QDPN(V_eq_zero)(out2[i], all); }
  QOPP(mgDiaginv)(out1, in, 1, da, QOP_EVEN);
  QOPP(mgDslash)(out2, out1, 1, da, QOP_ODD, QOP_EVEN);
  for(int i=0; i<nv; i++) { QDPN(V_eq_zero)(out1[i], all); }
  QOPP(mgDslashP)(out1, in, 1, da, QOP_ODD, QOP_EVEN);
  vcomp(out1, out2, nv, nc);
  for(int i=0; i<nv; i++) {
    QDPN(V_meq_V)(out2[i], out1[i], all);
    QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
  }
#endif

#if 0
  for(int i=0; i<nv; i++) { QDPN(V_eq_zero)(out2[i], all); }
  QOPP(mgDiaginv)(out1, in, 1, da, QOP_ODD);
  QOPP(mgDslash)(out2, out1, 1, da, QOP_EVEN, QOP_ODD);
  for(int i=0; i<nv; i++) { QDPN(V_eq_zero)(out1[i], all); }
  QOPP(mgDslashP)(out1, in, 1, da, QOP_EVEN, QOP_ODD);
  vcomp(out1, out2, nv, nc);
  for(int i=0; i<nv; i++) {
    QDPN(V_meq_V)(out2[i], out1[i], all);
    QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
  }
#endif

  if(nvin==nvout) {
    QOPP(mgDiaginv)(out1, in, 1, da, QOP_EVEN);
    QOPP(mgDiaginv)(out1, in, 1, da, QOP_ODD);
    QOPP(mgDslash)(out2, out1, 1, da, QOP_EVENODD, QOP_EVENODD);
    QOPP(mgDslashP)(out1, in, 1, da, QOP_EVENODD, QOP_EVENODD);
    vcomp(out1, out2, nvout, nc);
    for(int i=0; i<nvout; i++) {
      QDPN(V_meq_V)(out2[i], out1[i], all);
      QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
    }
  }

  if(nvin==nvout) {
    for(int i=0; i<nvout; i++) { QDPN(V_eq_zero)(out2[i], all); }
    QOPP(mgDslash)(out1, in, 1, da, QOP_EVENODD, QOP_EVENODD);
    QOPP(mgDslashEoProject)(out2, out1, da);
    for(int i=0; i<nvout; i++) { QDPN(V_eq_zero)(out1[i], all); }
    QOPP(mgDiag)(out3, in, 1, da, QOP_EVEN);
    QOPP(mgDslashEo)(out1, out3, 1, da);
    vcomp(out1, out2, nvout, nc);
#if 0
    for(int i=0; i<nvout; i++) {
      //QDPN(V_eq_funci)(out1[i], print_nonzeros, all);
      //QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
      QDPN(V_meq_V)(out2[i], out1[i], all);
      QDPN(V_eq_funci)(out2[i], print_nonzeros, all);
    }
#endif
  }

  for(int i=0; i<nvin; i++) {
    QDPN(destroy_V)(in[i]);
  }
  for(int i=0; i<nvout; i++) {
    QDPN(destroy_V)(out1[i]);
    QDPN(destroy_V)(out2[i]);
    QDPN(destroy_V)(out3[i]);
  }
}

#define qop2sub(q, e, o, a) (((q)==QOP_EVEN)?(e):(((q)==QOP_ODD)?(o):(a)))
#define PAROPP(p) (QOP_EVEN + QOP_ODD - p)
static void
QOPP(mgDslash2)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in0[], int sign,
	     int prec, void *args, QOP_evenodd_t parout, QOP_evenodd_t parin)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nvout = da->nvout;
  int nvin = da->nvin;
  int ns = da->nshifts;
  QDPN(ColorMatrix) ****links = da->links;
  QDPN(ColorVector) ***temp = da->temp;
  QDPN(ColorVector) **in = temp[ns];
  QDP_Shift *shifts = da->shifts;
  QDP_ShiftDir *fb = da->fb;
  QDP_Subset insub = qop2sub(parin, da->even, da->odd, da->all);
  QDP_Subset outsub = qop2sub(parout, da->even, da->odd, da->all);
  if(parout==QOP_ODD) temp = temp + ns + 2;
  if(prec) links = da->links + ns + 2;

  /* Start gathers from all directions */
  for(int j=0; j<nvin; j++) {
    if(sign>0 || j<nvin/2) { // assume gamma_5 Hermiticity in nvin
      QDPN(V_eq_V)(in[j], in0[j], insub);
    } else {
      QDPN(V_eqm_V)(in[j], in0[j], insub);
    }
    for(int i=0; i<ns; i++) {
      QDPN(V_eq_sV)(temp[i][j], in[j], shifts[i], fb[i], outsub);
    }
  }

  if(parin==PAROPP(parout)) {
    for(int j=0; j<nvout; j++) {
      QDPN(V_eq_zero)(out[j], outsub);
    }
  } else {
    if(prec) { // preconditioned operator
      for(int j=0; j<nvout; j++) {
	QDPN(V_eq_V)(out[j], in[j], outsub);
      }
    } else {
      for(int j=0; j<nvout; j++) {
	QDPN(V_eq_M_times_V)(out[j], links[ns][j][0], in[0], outsub);
	for(int k=1; k<nvin; k++) {
	  QDPN(V_peq_M_times_V)(out[j], links[ns][j][k], in[k], outsub);
	}
      }
    }
  }

  /* Multiply by matrix for all directions and accumulate */
  for(int j=0; j<nvout; j++) {
    //QDPN(V_eq_zero)(out[j], outsub);
    for(int k=0; k<nvin; k++) {
      for(int i=0; i<ns; i++) {
	QDPN(V_peq_M_times_V)(out[j], links[i][j][k], temp[i][k], outsub);
      }
    }
    if(sign<0 && j>=nvout/2) {
      QDPN(V_eqm_V)(out[j], out[j], outsub);
    }
  }
  for(int i=0; i<ns; i++) {
    for(int j=0; j<nvin; j++) {
      QDPN(discard_V)(temp[i][j]);
    }
  }
}

void
QOPP(mgDslash)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign,
	    void *args, QOP_evenodd_t parout, QOP_evenodd_t parin)
{
  if(QDP_volume_L(QDPN(get_lattice_V)(out[0]))==1) {
    if(parout==QOP_EVENODD) parout = QOP_EVEN;
    if(parin==QOP_EVENODD) parin = QOP_EVEN;
    if(parin==QOP_ODD || parout==QOP_ODD) {
      QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
      int nv = da->nvout;
      QDP_Subset sub = qop2sub(parout, da->even, da->odd, da->all);
      V_eq_zero(out);
      return;
    }
  }
  if(parout!=QOP_ODD) { // do even output sites
    if(parin==QOP_EVEN) {
      QOPP(mgDiag)(out, in, sign, args, QOP_EVEN);
    } else {
      QOPP(mgDslash2)(out, in, sign, 0, args, QOP_EVEN, parin);
    }
  }
  if(parout!=QOP_EVEN) { // do odd output sites
    if(parin==QOP_ODD) {
      QOPP(mgDiag)(out, in, sign, args, QOP_ODD);
    } else {
      QOPP(mgDslash2)(out, in, sign, 0, args, QOP_ODD, parin);
    }
  }
}

void
QOPP(mgDslashAll)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args)
{
  QOPP(mgDslash)(out, in, sign, args, QOP_EVENODD, QOP_EVENODD);
}

void
QOPP(mgDiag)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args,
	  QOP_evenodd_t par)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nvout = da->nvout;
  int nvin = da->nvin;
  int ns = da->nshifts;
  QDPN(ColorMatrix) ****links = da->links;
  QDP_Subset sub = qop2sub(par, da->even, da->odd, da->all);

  for(int j=0; j<nvout; j++) {
    if(sign>0 || j<nvout/2) { // assume gamma_5 Hermiticity in nv
      QDPN(V_eq_M_times_V)(out[j], links[ns][j][0], in[0], sub);
    } else {
      QDPN(V_eqm_M_times_V)(out[j], links[ns][j][0], in[0], sub);
    }
    for(int k=1; k<nvin; k++) {
      if(sign>0 || (j<nvout/2 && k<nvin/2) || (j>=nvout/2 && k>=nvin/2)) {
	QDPN(V_peq_M_times_V)(out[j], links[ns][j][k], in[k], sub);
      } else {
	QDPN(V_meq_M_times_V)(out[j], links[ns][j][k], in[k], sub);
      }
    }
  }
}

void
QOPP(mgDiaginv)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args,
	     QOP_evenodd_t par)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nvout = da->nvout;
  int nvin = da->nvin;
  int ns = da->nshifts + 1;
  QDPN(ColorMatrix) ****links = da->links;
  QDP_Subset sub = qop2sub(par, da->even, da->odd, da->all);

  for(int j=0; j<nvout; j++) {
    if(sign>0 || j<nvout/2) { // assume gamma_5 Hermiticity in nv
      QDPN(V_eq_M_times_V)(out[j], links[ns][j][0], in[0], sub);
    } else {
      QDPN(V_eqm_M_times_V)(out[j], links[ns][j][0], in[0], sub);
    }
    for(int k=1; k<nvin; k++) {
      if(sign>0 || (j<nvout/2 && k<nvin/2) || (j>=nvout/2 && k>=nvin/2)) {
	QDPN(V_peq_M_times_V)(out[j], links[ns][j][k], in[k], sub);
      } else {
	QDPN(V_meq_M_times_V)(out[j], links[ns][j][k], in[k], sub);
      }
    }
  }
}

void
QOPP(mgDslashP)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args,
		QOP_evenodd_t parout, QOP_evenodd_t parin)
{
  if(QDP_volume_L(QDPN(get_lattice_V)(out[0]))==1) {
    if(parout==QOP_EVENODD) parout = QOP_EVEN;
    if(parin==QOP_EVENODD) parin = QOP_EVEN;
    if(parin==QOP_ODD || parout==QOP_ODD) {
      QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
      int nv = da->nvout;
      QDP_Subset sub = qop2sub(parout, da->even, da->odd, da->all);
      V_eq_zero(out);
      return;
    }
  }
  if(parout!=QOP_ODD) { // do even output sites
    if(parin==QOP_EVEN) {
      QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
      int nv = da->nvout;
      QDP_Subset sub = da->even;
      V_eq_V(out, in);
    } else {
      QOPP(mgDslash2)(out, in, sign, 1, args, QOP_EVEN, parin);
    }
  }
  if(parout!=QOP_EVEN) { // do odd output sites
    if(parin==QOP_ODD) {
      QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
      int nv = da->nvout;
      QDP_Subset sub = da->odd;
      V_eq_V(out, in);
    } else {
      QOPP(mgDslash2)(out, in, sign, 1, args, QOP_ODD, parin);
    }
  }
}

void
QOPP(mgDslashPAll)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args)
{
  QOPP(mgDslashP)(out, in, sign, args, QOP_EVENODD, QOP_EVENODD);
}

void
QOPP(mgDslashEo)(QDPN(ColorVector) *out[], QDPN(ColorVector) *in[], int sign, void *args)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nv = da->nvout;
  QDP_Subset sub = da->even;
  QDPN(ColorVector) **t = da->temp[da->nshifts+1];
  QOPP(mgDslashP)(t, in, sign, args, QOP_ODD, QOP_EVEN);
  QOPP(mgDslashP)(out, t, sign, args, QOP_EVEN, QOP_ODD);
  V_eq_V_minus_V(out, in, out);
}

void
QOPP(mgDslashEoProject)(QDPN(ColorVector) *ineo[], QDPN(ColorVector) *in[], void *args)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nv = da->nvin;
  QDP_Subset sub = da->even;
  //{ QLA_Real nrm2; sub = da->odd; r_eq_norm2_V(&nrm2, in); printf0("nrm2 in = %g\n", nrm2); sub = da->even; }
  QOPP(mgDslashP)(ineo, in, 1, da, QOP_EVEN, QOP_ODD);
  //{ QLA_Real nrm2; r_eq_norm2_V(&nrm2, in); printf0("nrm2 ineo = %g\n", nrm2); }
  V_eq_V_minus_V(ineo, in, ineo);
  //{ QLA_Real nrm2; r_eq_norm2_V(&nrm2, in); printf0("nrm2 ineo = %g\n", nrm2); }
}

void
QOPP(mgDslashEoReconstruct)(QDPN(ColorVector) *out[], QDPN(ColorVector) *outeo[],
			    QDPN(ColorVector) *in[], void *args)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nv = da->nvin;
  QDP_Subset sub = da->odd;
  QDPN(ColorVector) **t = da->temp[da->nshifts+1];
  QOPP(mgDiaginv)(out, outeo, 1, da, QOP_EVEN);
  QOPP(mgDslash)(t, out, 1, da, QOP_ODD, QOP_EVEN);
  V_eq_V_minus_V(t, in, t);
  QOPP(mgDiaginv)(out, t, 1, da, QOP_ODD);
}

void
QOPP(mgDslashEoReconstructP)(QDPN(ColorVector) *out[], QDPN(ColorVector) *outeo[],
			     QDPN(ColorVector) *in[], void *args)
{
  QOPP(MgDslashArgs) *da = (QOPP(MgDslashArgs) *) args;
  int nv = da->nvin;
  QDP_Subset sub = da->even;
  V_eq_V(out, outeo);
  QOPP(mgDslashP)(out, out, 1, da, QOP_ODD, QOP_EVEN);
  sub = da->odd;
  V_eq_V_minus_V(out, in, out);
}
