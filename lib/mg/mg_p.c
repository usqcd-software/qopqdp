#include <qop_internal.h>

#include <math.h>
#include <float.h>
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

static int g_ic;
static QLA_Complex *g_z;

static void
V_eq_elem_c(int nc, QLAN(ColorVector)(nc,(*v)), int index)
{
  //QOP_printf0("g_nc = %i  i = %i\n", g_nc, i);
  //QOP_printf0("v = %p\n", v);
  //fflush(stdout);
  //QOP_printf0("g_z = %g %g\n", QLA_real(*g_z), QLA_imag(*g_z));
  //fflush(stdout);
  QLA_c_eq_c(QLA_elem_V(*v, g_ic), *g_z);
}

static void
QDPN(V_eq_elem_c_multi)(QDPN(ColorVector) *v, int nc, int ic,
			QLA_Complex z[], QDP_Subset subs[], int ns)
{
  g_ic = ic;
  for(int s=0; s<ns; s++) {
    g_z = &z[s];
    //QOP_printf0("g_z = %p\n", g_z);
    //QOP_printf0("z[s] = %g %g\n", QLA_real(z[s]), QLA_imag(z[s]));
    QDPN(V_eq_funcit)(v, V_eq_elem_c, subs[s]);
  }
}

#if 0
static void
QDPN(V_eq_r_times_V_multi)(QDPN(ColorVector) *r, QLA_Real z[],
			   QDPN(ColorVector) *v, QDP_Subset subs[], int ns)
{
  for(int s=0; s<ns; s++) {
    QDPN(V_eq_r_times_V)(r, &z[s], v, subs[s]);
  }
}
#endif

static void
QDPN(V_peq_c_times_V_multi)(QDPN(ColorVector) *r, QLA_Complex z[],
			    QDPN(ColorVector) *v, QDP_Subset subs[], int ns)
{
  for(int s=0; s<ns; s++) {
    QDPN(V_peq_c_times_V)(r, &z[s], v, subs[s]);
  }
}

#if 0
static void
QDPN(V_meq_c_times_V_multi)(QDPN(ColorVector) *r, QLA_Complex z[],
			    QDPN(ColorVector) *v, QDP_Subset subs[], int ns)
{
  for(int s=0; s<ns; s++) {
    QDPN(V_meq_c_times_V)(r, &z[s], v, subs[s]);
  }
}
#endif

static void
inv_sqrt_re(int nc, QLAN(ColorVector)(nc,(*v)), int index)
{
  for(int c=0; c<nc; c++) {
    QLA_D_Complex z;
    QLA_c_eq_c(z, QLA_elem_V(*v,c));
    QLA_c_eq_r(z, 1/sqrt(QLA_real(z)));
    QLA_c_eq_c(QLA_elem_V(*v,c), z);
  }
}

#if 0
static void
inv(int nc, QLAN(ColorVector)(nc,(*v)), int index)
{
  for(int c=0; c<nc; c++) {
    QLA_D_Complex z, z1, z2;
    QLA_c_eq_c(z, QLA_elem_V(*v,c));
    QLA_c_eq_r(z1, 1);
    QLA_c_eq_c_div_c(z2, z1, z);
    QLA_c_eq_c(QLA_elem_V(*v,c), z2);
  }
}
#endif

// globally orthogonalize all vectors against previous ones
void
QOPP(mgOrthoVn)(QDPN(ColorVector) ***v, int nv, int i0, int n, QDP_Subset sub)
{
  QLA_D_Complex z;
  QLA_D_Real r;
  QDPN(ColorVector) *cv1[nv], *cv2[nv];

  for(int i=0; i<n; i++) {
    for(int k=0; k<nv; k++) cv1[k] = v[k][i0+i];
    for(int j=0; j<i; j++) {
      for(int k=0; k<nv; k++) cv2[k] = v[k][i0+j];
      c_eq_V_dot_V(&z, cv2, cv1);
      V_meq_c_times_V(cv1, &z, cv2);
    }
    r_eq_norm2_V(&r, cv1);
    r = 1/sqrt(r);
    V_eq_r_times_V(cv1, &r, cv1);
  }
}

#if 0
static inline dcmplx dcmplxRI(double r, double i) { union U { double d[2]; dcmplx z; }; return (union U){{r,i}}.z; }
#define q2c(x) dcmplxRI(QLA_real(x),QLA_imag(x))
#define c2q(x,y) QLA_c_eq_r_plus_ir(x, creal(y), cimag(y))

void
QOPP(rRitzHarm)(QDPN(ColorVector) *vv[], int n, int nv, QOP_MgOp *op, void *opargs, QLA_Complex evs[], int sort, int norm, QDP_Subset sub)
{
  QDPN(ColorVector) *Av[n][nv], *(*v)[nv] = (QDPN(ColorVector) *(*)[nv]) vv;
  QDP_Lattice *lat = QDPN(get_lattice_V)(v[0][0]);
  int nc = QDP_get_nc(v[0][0]);
  zmat m;
  zmat_alloc(&m, n, n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<nv; j++) {
      Av[i][j] = QDPN(create_V_L)(nc, lat);
    }
    op(Av[i], v[i], 1, opargs);
    for(int k=0; k<i; k++) {
      QLA_Complex z;
      c_eq_V_dot_V(&z, Av[k], Av[i]);
      V_meq_c_times_V(Av[i], &z, Av[k]);
      V_meq_c_times_V(v[i], &z, v[k]);
    }
    QLA_Real nA;
    r_eq_norm2_V(&nA, Av[i]);
    QOP_printf0("|Av[%i]|^2 = %g\n", i, nA);
    nA = 1/sqrt(nA);
    V_eq_r_times_V(v[i], &nA, v[i]);
    V_eq_r_times_V(Av[i], &nA, Av[i]);
    //r_eq_norm2_V(&nA, Av[i]);
    //QOP_printf0("|Av[%i]|^2 = %g\n", i, nA);
    for(int k=0; k<=i; k++) {
      QLA_Complex z;
      c_eq_V_dot_V(&z, Av[i], v[k]);
      zmat_set(&m, i, k, q2c(z));
    }
    for(int k=0; k<i; k++) {
      QLA_Complex z;
      c_eq_V_dot_V(&z, Av[k], v[i]);
      zmat_set(&m, k, i, q2c(z));
    }
  }

  zmat vl, vr;
  zvec e;
  zmat_alloc(&vl, n, n);
  zmat_alloc(&vr, n, n);
  zvec_alloc(&e, n);
  zgeigsv(&m, &e, &vl, &vr, n);

  // sort ascending
  QLA_Real ae[n];
  int idx[n];
  for(int i=0; i<n; i++) {
    int ii = n-1-i;
    dcmplx z = 1/zvec_get(&e,ii);
    QLA_Real ai = cabs(z);
    //QLA_Real ai = sqrt(DBL_EPSILON+creal(z)*creal(z));
    int k = i;
    while(k>0 && ae[k-1]>ai) {
      idx[k] = idx[k-1];
      evs[k] = evs[k-1];
      ae[k] = ae[k-1];
      k--;
    }
    idx[k] = ii;
    c2q(evs[k], z);
    ae[k] = ai;
  }

  // E norm for now
  for(int i=0; i<n; i++) {
    //QLA_Real s;
    //s = 1/sqrt(ae[i]);
    //V_eq_zero(v[i]);
    V_eq_zero(Av[i]);
    int ii = idx[i];
    for(int k=0; k<n; k++) {
      QLA_Complex z;
      //c2q(z, s*zmat_get(&vr,k,ii));
      c2q(z, zmat_get(&vr,k,ii));
      //V_peq_c_times_V(v[i], &z, Av[k]);
      V_peq_c_times_V(Av[i], &z, v[k]);
    }
  }
  for(int i=0; i<n; i++) {
    QLA_Real s, nA;
    r_eq_norm2_V(&nA, Av[i]);
    switch(norm) {
    case QOP_NORM_E: s = 1/sqrt(nA*ae[i]); break;
    case QOP_NORM_2:
    default: s = 1/sqrt(nA);
    }
    V_eq_r_times_V(v[i], &s, Av[i]);
  }
  zmat_free(&m);
  zmat_free(&vr);
  zmat_free(&vl);
  zvec_free(&e);

  for(int i=0; i<n; i++) {
    QLA_Real nrm2;
    r_eq_norm2_V(&nrm2, v[i]);
    QOP_printf0("|v[%i]|^2 = %g\n", i, nrm2);
    for(int j=0; j<nv; j++) {
      QDPN(destroy_V)(Av[i][j]);
    }
  }
}
#endif

// locally orthogonalize all vectors against previous ones
void
QOPP(mgOrtho)(QDPN(ColorVector) *cv[], int nv, QOP_MgBlock *mgb)
{
#if 1
  int fnc = QDP_get_nc(cv[0]);
  QDPN(ColorVector) *t = QDPN(create_V_L)(fnc, mgb->fine);
  QDPN(ColorVector) *z = QDPN(create_V_L)(1, mgb->coarse);
  QDPN(V_eq_zero)(t, QDP_all_L(mgb->fine));
  QDPN(V_eq_zero)(z, QDP_all_L(mgb->coarse));
  for(int i=0; i<nv; i++) {
    for(int j=0; j<i; j++) {
      //QDPN(c_eq_V_dot_V_multi)(z, cv[j], cv[i], mgb->fs, ns);
      QOPP(mgRestrict)(&z, &cv[i], 1, &cv[j], 1, fnc, mgb, QOP_EVENODD);
      //QDPN(V_meq_c_times_V_multi)(cv[i], z, cv[j], mgb->fs, ns);
      QOPP(mgProlong)(&t, &z, 1, &cv[j], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_meq_V)(cv[i], t, QDP_all_L(mgb->fine));
    }
    //QDPN(r_eq_norm2_V_multi)(r, cv[i], mgb->fs, ns);
    QOPP(mgRestrict)(&z, &cv[i], 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
    //for(int j=0; j<ns; j++) r[j] = 1/sqrt(r[j]);
    QDPN(V_eq_funcit)(z, inv_sqrt_re, QDP_all_L(mgb->coarse));
    //QDPN(V_eq_r_times_V_multi)(cv[i], r, cv[i], mgb->fs, ns);
    QOPP(mgProlong)(&t, &z, 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
    QDPN(V_eq_V)(cv[i], t, QDP_all_L(mgb->fine));
  }
  QDPN(destroy_V)(t);
  QDPN(destroy_V)(z);
#else
  QLA_Real *r;
  QLA_Complex *z;
  int ns = mgb->ns;
  QOP_malloc(r, QLA_Real, ns);
  QOP_malloc(z, QLA_Complex, ns);
  for(int i=0; i<nv; i++) {
    for(int j=0; j<i; j++) {
      QDPN(c_eq_V_dot_V_multi)(z, cv[j], cv[i], mgb->fs, ns);
      QDPN(V_meq_c_times_V_multi)(cv[i], z, cv[j], mgb->fs, ns);
    }
    QDPN(r_eq_norm2_V_multi)(r, cv[i], mgb->fs, ns);
    for(int j=0; j<ns; j++) {
      r[j] = 1/sqrt(r[j]);
    }
    QDPN(V_eq_r_times_V_multi)(cv[i], r, cv[i], mgb->fs, ns);
  }
  QOP_free(r);
  QOP_free(z);
#endif
}

static void
ones_func(int nc, QLAN(ColorVector)(nc,(*v)), int index)
{
  for(int i=0; i<nc; i++) {
    QLA_c_eq_r(QLA_elem_V(*v,i), 1);
  }
}

static void
get_index_larger(QDPN(ColorVector) *zk, QDPN(ColorVector) *zn, QDPN(ColorVector) *z, int j)
{
  int n = QDP_subset_len(QDP_all_L(QDPN(get_lattice_V)(zk)));
  QLAN(ColorVector)(1, (*qzk)) = QDPN(expose_V)(zk);
  QLAN(ColorVector)(1, (*qzn)) = QDPN(expose_V)(zn);
  QLAN(ColorVector)(1, (*qz)) = QDPN(expose_V)(z);
#pragma omp parallel for
  for(int i=0; i<n; i++) {
    QLA_Real rzn = QLA_real(QLA_elem_V(qzn[i],0));
    QLA_Real rz = QLA_real(QLA_elem_V(qz[i],0));
    if(rz>=rzn) {
      QLA_c_eq_r(QLA_elem_V(qzk[i],0), j);
      QLA_c_eq_c(QLA_elem_V(qzn[i],0), QLA_elem_V(qz[i],0));
    }
  }
  QDPN(reset_V)(zk);
  QDPN(reset_V)(zn);
  QDPN(reset_V)(z);
#if 0
  int i;
  QDP_loop_sites(i, QDP_all_L(QDPN(get_lattice_V)(zk)), {
    QLAN(ColorVector)(1, (*qzn));
    qzn = QDPN(site_ptr_readonly_V)(zn,i);
    QLAN(ColorVector)(1, (*qz)) = QDPN(site_ptr_readonly_V)(z,i);
    QLA_Real rzn = QLA_real(QLA_elem_V(*qzn,0));
    QLA_Real rz = QLA_real(QLA_elem_V(*qz,0));
    if(rz>=rzn) {
      QLAN(ColorVector)(1, (*qzk)) = QDPN(site_ptr_readwrite_V)(zk,i);
      QLA_c_eq_r(QLA_elem_V(*qzk,0), j);
      }
  });
#endif
}

static void
swap_sites(QDPN(ColorVector) *cv[], int nv, int nc, QDPN(ColorVector) *tk, int j)
{
  int n = QDP_subset_len(QDP_all_L(QDPN(get_lattice_V)(tk)));
  QLAN(ColorVector)(1, (*qtk)) = QDPN(expose_V)(tk);
  QLAN(ColorVector)(nc, (*qcv[nv]));
  for(int i=0; i<nv; i++) qcv[i] = QDPN(expose_V)(cv[i]);
#pragma omp parallel for
  for(int i=0; i<n; i++) {
    int k = round(QLA_real(QLA_elem_V(qtk[i],0)));
    if(k!=j) {
      QLAN(ColorVector)(nc,t);
      QLAN(V_eq_V)(nc, &t, &qcv[j][i]);
      QLAN(V_eq_V)(nc, &qcv[j][i], &qcv[k][i]);
      QLAN(V_eq_V)(nc, &qcv[k][i], &t);
    }
  }
  QDPN(reset_V)(tk);
  for(int i=0; i<nv; i++) QDPN(reset_V)(cv[i]);
#if 0
  int i;
  QDP_loop_sites(i, QDP_all_L(QDPN(get_lattice_V)(tk)), {
    QLAN(ColorVector)(1, (*qtk)) = QDPN(site_ptr_readonly_V)(tk,i);
    int k = round(QLA_real(QLA_elem_V(*qtk,0)));
    if(k!=j) {
      QLAN(ColorVector)(nc,t);
      QLAN(V_eq_V)(nc, &t, QDPN(site_ptr_readonly_V)(cv[j],i));
      QLAN(V_eq_V)(nc, QDPN(site_ptr_readwrite_V)(cv[j],i), QDPN(site_ptr_readonly_V)(cv[k],i));
      QLAN(V_eq_V)(nc, QDPN(site_ptr_readwrite_V)(cv[k],i), &t);
    }
  });
#endif
}

static double g_min, g_sum, g_max, g_n;
static void
stats_func(int nc, QLAN(ColorVector)(nc,(*v)), int index)
{
  QLA_Real r;
  QLAN(r_eq_norm2_V)(nc, &r, v);
  r = sqrt(r);
  g_n += 1;
  g_sum += r;
  if(r<g_min) g_min = r;
  if(r>g_max) g_max = r;
}

static void
get_stats(QDPN(ColorVector) *v, double *min, double *ave, double *max)
{
  g_min = DBL_MAX;
  g_sum = g_max = g_n = 0;
  QDPN(V_eq_funci)(v, stats_func, QDP_all_L(QDPN(get_lattice_V)(v)));
  QMP_min_double(&g_min);
  QMP_sum_double(&g_n);
  QMP_sum_double(&g_sum);
  QMP_max_double(&g_max);
  *min = g_min;
  *ave = g_sum/g_n;
  *max = g_max;
}

// locally orthogonalize all vectors against previous ones
//  and sort
void
QOPP(mgOrthoSort)(QDPN(ColorVector) *cv[], int imin, int n, QOP_MgBlock *mgb,
		  double min[], double ave[], double max[])
{
  int fnc = QDP_get_nc(cv[0]);
  QDPN(ColorVector) *t = QDPN(create_V_L)(fnc, mgb->fine);
  QDPN(ColorVector) *tk = QDPN(create_V_L)(1, mgb->fine);
  QDPN(ColorVector) *ones = QDPN(create_V_L)(1, mgb->fine);
  QDPN(ColorVector) *z = QDPN(create_V_L)(1, mgb->coarse);
  QDPN(ColorVector) *zk = QDPN(create_V_L)(1, mgb->coarse);
  QDPN(ColorVector) *zn = QDPN(create_V_L)(1, mgb->coarse);
  QDPN(V_eq_zero)(t, QDP_all_L(mgb->fine));
  QDPN(V_eq_zero)(z, QDP_all_L(mgb->coarse));
  QDPN(V_eq_funcit)(ones, ones_func, QDP_all_L(mgb->fine));
  for(int i=0; i<imin; i++) {
    for(int j=imin; j<n; j++) {
      QOPP(mgRestrict)(&z, &cv[j], 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QOPP(mgProlong)(&t, &z, 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_meq_V)(cv[j], t, QDP_all_L(mgb->fine));
    }
  }
  for(int i=imin; i<n; i++) {
    // find remaining local vector with largest norm
    QDPN(V_eq_zero)(zk, QDP_all_L(mgb->coarse));
    QDPN(V_eq_zero)(zn, QDP_all_L(mgb->coarse));
    for(int j=i; j<n; j++) {
      QOPP(mgRestrict)(&z, &cv[j], 1, &cv[j], 1, fnc, mgb, QOP_EVENODD);
      get_index_larger(zk, zn, z, j);
    }
    get_stats(zn, &min[i], &ave[i], &max[i]);
    QOPP(mgProlong)(&tk, &zk, 1, &ones, 1, 1, mgb, QOP_EVENODD); // really should be 1->X shift
    swap_sites(cv, n, fnc, tk, i);
    {
      QLA_Real nrm2;
      QDPN(r_eq_norm2_V)(&nrm2, cv[i], QDP_all_L(mgb->fine));
      QOP_printf0("|cv[%i]|^2 = %g\n", i, nrm2);
      QOPP(mgRestrict)(&z, &cv[i], 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_meq_V)(z, zn, QDP_all_L(mgb->coarse));
      QDPN(r_eq_norm2_V)(&nrm2, z, QDP_all_L(mgb->coarse));
      QOP_printf0("|z-zn|^2 = %g\n", nrm2);
    }

    // normalize
    QDPN(V_eq_funcit)(zn, inv_sqrt_re, QDP_all_L(mgb->coarse));
    QOPP(mgProlong)(&t, &zn, 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
    QDPN(V_eq_V)(cv[i], t, QDP_all_L(mgb->fine));
    {
      QLA_Real nrm2;
      QDPN(r_eq_norm2_V)(&nrm2, cv[i], QDP_all_L(mgb->fine));
      QOP_printf0("|cv[%i]|^2 = %g\n", i, nrm2);
    }

    // project out from remaining vectors
    for(int j=i+1; j<n; j++) {
      QOPP(mgRestrict)(&z, &cv[j], 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QOPP(mgProlong)(&t, &z, 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_meq_V)(cv[j], t, QDP_all_L(mgb->fine));
    }
  }
  QDPN(destroy_V)(t);
  QDPN(destroy_V)(tk);
  QDPN(destroy_V)(ones);
  QDPN(destroy_V)(z);
  QDPN(destroy_V)(zk);
  QDPN(destroy_V)(zn);
}


// locally orthogonalize the nv'th vector against previous ones
// assumes previous ones are already locally orthonormalized
void
QOPP(mgOrthoVec)(QDPN(ColorVector) *cv[], int nv, QOP_MgBlock *mgb, int norm)
{
  if(nv>0 || norm) {
    int fnc = QDP_get_nc(cv[0]);
    QDPN(ColorVector) *t = QDPN(create_V_L)(fnc, mgb->fine);
    QDPN(ColorVector) *z = QDPN(create_V_L)(1, mgb->coarse);
    QDPN(V_eq_zero)(t, QDP_all_L(mgb->fine));
    QDPN(V_eq_zero)(z, QDP_all_L(mgb->coarse));
    int i=nv;
    for(int j=0; j<i; j++) {
      //QDPN(c_eq_V_dot_V_multi)(z, cv[j], cv[i], mgb->fs, ns);
      QOPP(mgRestrict)(&z, &cv[i], 1, &cv[j], 1, fnc, mgb, QOP_EVENODD);
      //QDPN(V_meq_c_times_V_multi)(cv[i], z, cv[j], mgb->fs, ns);
      QOPP(mgProlong)(&t, &z, 1, &cv[j], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_meq_V)(cv[i], t, QDP_all_L(mgb->fine));
    }
    if(norm) {
      //QDPN(r_eq_norm2_V_multi)(r, cv[i], mgb->fs, ns);
      QOPP(mgRestrict)(&z, &cv[i], 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      //for(int j=0; j<ns; j++) r[j] = 1/sqrt(r[j]);
      QDPN(V_eq_funcit)(z, inv_sqrt_re, QDP_all_L(mgb->coarse));
      //QDPN(V_eq_r_times_V_multi)(cv[i], r, cv[i], mgb->fs, ns);
      QOPP(mgProlong)(&t, &z, 1, &cv[i], 1, fnc, mgb, QOP_EVENODD);
      QDPN(V_eq_V)(cv[i], t, QDP_all_L(mgb->fine));
    }
    QDPN(destroy_V)(t);
    QDPN(destroy_V)(z);
  }
}

#if 0
static void
QDPN(Minv)(QDPN(ColorMatrix) *m, int nc)
{
  QLAN(ColorMatrix)(nc, (*qm));
  qm = QDPN(expose_M)(m);
  int ns = QDP_sites_on_node_L(QDPN(get_lattice_M)(m));

  zmat t;
  zmat_alloc(&t, nc, nc);

  for(int s=0; s<ns; s++) {
    for(int i=0; i<nc; i++) {
      for(int j=0; j<nc; j++) {
	cmplx *z = (cmplx *)&QLA_elem_M(qm[s], i, j);
	zmat_set(&t, i, j, (dcmplx)*z);
      }
    }
    zinv(&t);
    for(int i=0; i<nc; i++) {
      for(int j=0; j<nc; j++) {
	dcmplx z = zmat_get(&t, i, j);
	c2q(QLA_elem_M(qm[s], i, j), z);
      }
    }
  }

  QDPN(reset_M)(m);
  zmat_free(&t);
}

// orthonormalizes cv2 against cv1 on MG blocks
void
QOPP(mgOrtho2)(QDPN(ColorVector) *cv1[], QDPN(ColorVector) *cv2[], int nv, QOP_MgBlock *mgb)
{
  int fnc = QDP_get_nc(cv1[0]);
  QDPN(ColorMatrix) *m = QDPN(create_M_L)(nv, mgb->coarse);
  QDPN(ColorVector) *z = QDPN(create_V_L)(nv, mgb->coarse);
  QDPN(ColorVector) *t[nv];

  QDPN(V_eq_zero)(z, QDP_all_L(mgb->coarse));
  for(int i=0; i<nv; i++) {
    t[i] = QDPN(create_V_L)(fnc, mgb->fine);
    QDPN(V_eq_zero)(t[i], QDP_all_L(mgb->fine));
  }

  for(int i=0; i<nv; i++) {
    QOPP(mgRestrict)(&z, &cv2[i], 1, cv1, nv, fnc, mgb, QOP_EVENODD);
    QDPN(M_eq_colorvec_V)(m, z, i, QDP_all_L(mgb->coarse));
  }
  QDPN(Minv)(m, nv);
  for(int i=0; i<nv; i++) {
    QDPN(V_eq_colorvec_M)(z, m, i, QDP_all_L(mgb->coarse));
    QOPP(mgProlong)(&t[i], &z, 1, cv2, nv, fnc, mgb, QOP_EVENODD);
  }
  for(int i=0; i<nv; i++) {
    QDPN(V_eq_V)(cv2[i], t[i], QDP_all_L(mgb->fine));
    QDPN(destroy_V)(t[i]);
  }
  QDPN(destroy_V)(z);
  QDPN(destroy_M)(m);
}
#endif

// cv[i][block][color] = pv[color][block]' * fv[i][block]
// i=0..nv-1, color=0..cnc-1
void
QOPP(mgRestrict)(QDPN(ColorVector) *cv[], QDPN(ColorVector) *fv[], int nv,
		 QDPN(ColorVector) *pv[], int cnc, int fnc, QOP_MgBlock *mgb,
		 QOP_evenodd_t par)
{
  double t0 = QDP_time();
#define CHKNC(nc,v) { int _nc = QDP_get_nc(v); if(nc!=_nc) { QOP_printf0("%s mismatch in %s %s %i %i\n", __func__, #nc, #v, nc, _nc); QDP_abort(1); } }
  //CHKNC(cnc, cv[0]);
  //CHKNC(fnc, fv[0]);
  //CHKNC(fnc, pv[0]);
  if(mgb->local) {
    QLAN(ColorVector)(cnc, (*qcv[nv]));
    QLAN(ColorVector)(fnc, (*qfv[nv]));
    QLAN(ColorVector)(fnc, (*qpv[cnc]));
    for(int i=0; i<nv; i++) {
      qcv[i] = QDPN(expose_V)(cv[i]);
      qfv[i] = QDPN(expose_V)(fv[i]);
    }
    for(int cc=0; cc<cnc; cc++) {
      if(pv!=fv) {
	qpv[cc] = QDPN(expose_V)(pv[cc]);
      } else {
	qpv[cc] = qfv[cc];
      }
    }
    //printf("size = %i\n", sizeof(qcv[i][b]));
    for(int b=0; b<mgb->nlb; b++) {
      QOP_MgLocalBlock *lb = &(mgb->lb[b]);
      int *sites = par==QOP_ODD ? lb->sites[1] : lb->sites[0];
      int nsites = (par==QOP_ODD ? 0 : lb->nsites[0]) + (par==QOP_EVEN ? 0 : lb->nsites[1]);
      for(int i=0; i<nv; i++) {
	for(int cc=0; cc<cnc; cc++) {
	  QLA_Complex z;
	  QLAN(c_xeq_V_dot_V)(fnc, &z, qpv[cc], qfv[i], sites, nsites);
	  //QOP_printf0("%i %i %p\n", i, b, &qcv[i][b]);
	  QLA_c_eq_c(QLA_elem_V(qcv[i][b],cc), z);
	}
      }
    }
    for(int i=0; i<nv; i++) {
      QDPN(reset_V)(cv[i]);
      QDPN(reset_V)(fv[i]);
    }
    if(pv!=fv) {
      for(int cc=0; cc<cnc; cc++) {
	QDPN(reset_V)(pv[cc]);
      }
    }
  } else {
    QLA_Complex *z;
    int ns = mgb->ns;
    //QLA_Complex z[ns];
    QOP_malloc(z, QLA_Complex, ns);
    //QOP_printf0("ns = %i\n", ns);
    for(int i=0; i<nv; i++) {
      for(int cc=0; cc<cnc; cc++) {
	//QOP_printf0("z[0] = %g %g\n", QLA_real(z[0]), QLA_imag(z[0])); fflush(stdout);
	QDPN(c_eq_V_dot_V_multi)(z, pv[cc], fv[i], mgb->fs, ns);
	if(0) for(int s=0; s<ns; s++) {
	    QLA_Complex z2, z3, z4;
	    QDPN(c_eq_V_dot_V)(&z2, pv[cc], fv[i], mgb->fs[s]);
	    QLA_c_eq_c(z3, z[s]);
	    QLA_c_eq_c_minus_c(z4, z2, z3);
	    if(QLA_norm2_c(z4)>1e-6) {
#define PZ(z) QLA_real(z), QLA_imag(z)
	      QOP_printf0("* %i %i %i : %g %g  %g %g\n", s, cc, i, PZ(z2), PZ(z3));
	    }
	  }
	//QOP_printf0("z[0] = %g %g\n", QLA_real(z[0]), QLA_imag(z[0])); fflush(stdout);
	QDPN(V_eq_elem_c_multi)(cv[i], cnc, cc, z, mgb->cs, ns);
      }
    }
    QOP_free(z);
  }
  if(0) {
    int ns = mgb->ns;
    for(int s=0; s<ns; s++) {
      for(int i=0; i<nv; i++) {
	QLAN(ColorVector)(cnc, zv);
	QDPN(v_eq_sum_V)(&zv, cv[i], mgb->cs[s]);
	for(int cc=0; cc<cnc; cc++) {
	  QLA_Complex z2, z3, z4;
	  QDPN(c_eq_V_dot_V)(&z2, pv[cc], fv[i], mgb->fs[s]);
	  QLA_c_eq_c(z3, QLA_elem_V(zv,cc));
	  QLA_c_eq_c_minus_c(z4, z2, z3);
	  if(QLA_norm2_c(z4)>1e-6) {
#define PZ(z) QLA_real(z), QLA_imag(z)
	    QOP_printf0("** %i %i %i : %g %g  %g %g\n", s, cc, i, PZ(z2), PZ(z3));
	  }
	}
      }
    }
  }
  mgb->trestrict += QDP_time() - t0;
  mgb->rcount++;
}

// fv[i][block] = pv[color][block] * cv[i][block][color]
// i=0..nv-1, color=0..cnc-1
void
QOPP(mgProlong)(QDPN(ColorVector) *fv[], QDPN(ColorVector) *cv[], int nv,
		QDPN(ColorVector) *pv[], int cnc, int fnc, QOP_MgBlock *mgb,
		QOP_evenodd_t par)
{
  double t0 = QDP_time();
  if(mgb->local) {
    QLAN(ColorVector)(cnc, (*qcv[nv]));
    QLAN(ColorVector)(fnc, (*qfv[nv]));
    QLAN(ColorVector)(fnc, (*qpv[cnc]));
    for(int i=0; i<nv; i++) {
      qcv[i] = QDPN(expose_V)(cv[i]);
      qfv[i] = QDPN(expose_V)(fv[i]);
    }
    for(int cc=0; cc<cnc; cc++) {
      if(pv!=fv) {
	qpv[cc] = QDPN(expose_V)(pv[cc]);
      } else {
	qpv[cc] = qfv[cc];
      }
    }
    for(int b=0; b<mgb->nlb; b++) {
      QOP_MgLocalBlock *lb = &(mgb->lb[b]);
      int *sites = par==QOP_ODD ? lb->sites[1] : lb->sites[0];
      int nsites = (par==QOP_ODD ? 0 : lb->nsites[0]) + (par==QOP_EVEN ? 0 : lb->nsites[1]);
      for(int i=0; i<nv; i++) {
	QLAN(V_xeq_zero)(fnc, qfv[i], sites, nsites);
	for(int cc=0; cc<cnc; cc++) {
	  QLA_Complex z;
	  QLA_c_eq_c(z, QLA_elem_V(qcv[i][b],cc));
	  if(QLA_norm2_c(z)>0) {
	    QLAN(V_xpeq_c_times_V)(fnc, qfv[i], &z, qpv[cc], sites, nsites);
	  }
	}
      }
    }
    for(int i=0; i<nv; i++) {
      QDPN(reset_V)(cv[i]);
      QDPN(reset_V)(fv[i]);
    }
    if(pv!=fv) {
      for(int cc=0; cc<cnc; cc++) {
	QDPN(reset_V)(pv[cc]);
      }
    }
  } else {
    QLA_Complex *z;
    int ns = mgb->ns;
    //QLA_Complex z[ns];
    QOP_malloc(z, QLA_Complex, ns);
    QDP_Complex *ct = QDP_create_C_L(mgb->coarse);
    for(int i=0; i<nv; i++) {
      QDPN(V_eq_zero)(fv[i], QDP_all_L(mgb->fine));
      for(int cc=0; cc<cnc; cc++) {
	QDPN(C_eq_elem_V)(ct, cv[i], cc, QDP_all_L(mgb->coarse));
	QDP_c_eq_sum_C_multi(z, ct, mgb->cs, ns);
	QDPN(V_peq_c_times_V_multi)(fv[i], z, pv[cc], mgb->fs, ns);
      }
    }
    QOP_free(z);
  }
  mgb->tprolong += QDP_time() - t0;
  mgb->pcount++;
}

QOPP(MgArgs) *
QOPP(mgCreateArgs)(QOP_MgBlock *mgb, int cnc, int fnc, int nv,
		   QDPN(ColorVector) ***rv, QDPN(ColorVector) ***pv)
{
  QOPP(MgArgs) *QOP_malloc(mga, QOPP(MgArgs), 1);

  mga->cnc = cnc;
  mga->fnc = fnc;
  mga->nv = nv;
  mga->maxnv = nv;
  mga->mgb = mgb;
  if(rv) {
    mga->rv = rv;
    mga->we_malloced_rv = 0;
  } else {
    QOP_malloc(rv, QDPN(ColorVector) **, nv);
    for(int i=0; i<nv; i++) {
      QOP_malloc(rv[i], QDPN(ColorVector) *, cnc);
      for(int j=0; j<cnc; j++) {
	rv[i][j] = QDPN(create_V_L)(fnc, mgb->fine);
      }
    }
    mga->rv = rv;
    mga->we_malloced_rv = 1;
  }
  if(pv) {
    mga->pv = pv;
    mga->we_malloced_pv = 0;
  } else {
    QOP_malloc(pv, QDPN(ColorVector) **, nv);
    for(int i=0; i<nv; i++) {
      QOP_malloc(pv[i], QDPN(ColorVector) *, cnc);
      for(int j=0; j<cnc; j++) {
	pv[i][j] = QDPN(create_V_L)(fnc, mgb->fine);
      }
    }
    mga->pv = pv;
    mga->we_malloced_pv = 1;
  }

  return mga;
}

void
QOPP(mgFreeArgs)(QOPP(MgArgs) *mga)
{
  if(mga->we_malloced_rv) {
    for(int i=0; i<mga->maxnv; i++) {
      for(int j=0; j<mga->cnc; j++) {
	QDPN(destroy_V)((mga->rv)[i][j]);
      }
      QOP_free((mga->rv)[i]);
    }
    QOP_free(mga->rv);
  }
  if(mga->we_malloced_pv) {
    for(int i=0; i<mga->maxnv; i++) {
      for(int j=0; j<mga->cnc; j++) {
	QDPN(destroy_V)((mga->pv)[i][j]);
      }
      QOP_free((mga->pv)[i]);
    }
    QOP_free(mga->pv);
  }
  QOP_free(mga);
}

void
QOPP(mgF2c)(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
	    QOPP(MgArgs) *w, QOP_evenodd_t par)
{
  QOP_MgBlock *mgb = w->mgb;
  int cnc = w->cnc;
  int fnc = w->fnc;
  int nv = w->nv;
  QDPN(ColorVector) ***rv = w->rv;

  for(int i=0; i<nv; i++) {
    QOPP(mgRestrict)(&out[i], &in[i], 1, rv[i], cnc, fnc, mgb, par);
  }
}

void
QOPP(mgC2f)(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
	    QOPP(MgArgs) *w, QOP_evenodd_t par)
{
  QOP_MgBlock *mgb = w->mgb;
  int cnc = w->cnc;
  int fnc = w->fnc;
  int nv = w->nv;
  QDPN(ColorVector) ***pv = w->pv;

  for(int i=0; i<nv; i++) {
    QOPP(mgProlong)(&out[i], &in[i], 1, pv[i], cnc, fnc, mgb, par);
  }
}

// takes fine vectors, goes to coarse, applies op, then brings back to fine
void
QOPP(mgF2cOp)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  QOPP(MgF2cOpArgs) *mgoa = (QOPP(MgF2cOpArgs) *) args;
  QOPP(MgArgs) *mga = mgoa->mga;
  QDPN(ColorVector) **cin = mgoa->cin;
  QDPN(ColorVector) **cout = mgoa->cout;

  QOPP(mgF2c)(cin, in, mga, mgoa->fpar);
  mgoa->op(cout, cin, sign, mgoa->opargs);
  QOPP(mgC2f)(out, cout, mga, mgoa->fpar);
}

// takes coarse vectors, goes to fine, applies op, then brings back to coarse
void
QOPP(mgC2fOp)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  QOPP(MgC2fOpArgs) *mgoa = (QOPP(MgC2fOpArgs) *) args;
  QOPP(MgArgs) *mga = mgoa->mga;
  QDPN(ColorVector) **fin = mgoa->fin;
  QDPN(ColorVector) **fout = mgoa->fout;

  QOPP(mgC2f)(fin, in, mga, mgoa->fpar);
  { QDP_Subset sub = QDP_all_L(mga->mgb->fine); int nv = mga->nv; V_eq_zero(fout); }
  mgoa->op(fout, fin, sign, mgoa->opargs);
  QOPP(mgF2c)(out, fout, mga, mgoa->fpar);
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
    QOP_printf0("vcomp %i  ( %g %g )  %g / %g = %g\n", i, v1n, v2n, dn, sn, dn/sn);
  }
  QDPN(destroy_V)(d);
  QDPN(destroy_V)(s);
}

void
QOPP(mgTestCoarse)(QOPP(MgArgs) *mga, QDPN(ColorVector) **vf)
{
  QOP_MgBlock *mgb = mga->mgb;
  int fnc = mga->fnc;
  int cnc = mga->cnc;
  int nv = mga->nv;
  QDP_Lattice *flat = mgb->fine;
  QDP_Lattice *clat = mgb->coarse;

  QDPN(ColorVector) *vc[nv];
  for(int i=0; i<nv; i++) {
    vc[i] = QDPN(create_V_L)(cnc, clat);
  }

  QOPP(mgF2c)(vc, vf, mga, QOP_EVENODD);
  QOPP(mgC2f)(vf, vc, mga, QOP_EVENODD);

  QDPN(ColorVector) *vc2[nv];
  for(int i=0; i<nv; i++) {
    vc2[i] = QDPN(create_V_L)(cnc, clat);
  }

  QOPP(mgF2c)(vc2, vf, mga, QOP_EVENODD);

  vcomp(vc, vc2, nv, cnc);

  QDPN(ColorVector) *vf2[nv];
  for(int i=0; i<nv; i++) {
    vf2[i] = QDPN(create_V_L)(fnc, flat);
  }

  QOPP(mgC2f)(vf2, vc2, mga, QOP_EVENODD);

  vcomp(vf, vf2, nv, fnc);

  for(int i=0; i<nv; i++) {
    QDPN(destroy_V)(vc[i]);
    QDPN(destroy_V)(vc2[i]);
    QDPN(destroy_V)(vf2[i]);
  }
}
