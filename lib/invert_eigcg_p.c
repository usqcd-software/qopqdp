/* eigCG of Stathopoulos and Orginos (arXiv:0707.0131) */
#ifndef HAVE_LAPACK

QOP_status_t
QOPPCV(invert_eigcg)(QOPPCV(linop_t) *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     Vector *out,
		     Vector *in,
		     Vector *p,
		     QDP_Subset subset,
		     QOPPCV(eigcg_t) *eigcg
		     vIndexDef)
{
  QOP_error("need to compile QOP with --enable-lapack for eigcg!\n");
  return QOP_FAIL;
}

#else

#include <complex.h>
#include <float.h>
#include "linalg.h"

#define printf0 QOP_printf0
//#define TRACE printf0("%s %s %i\n", __FILE__, __func__, __LINE__);
#define TRACE

#if 1
#define BEGIN_TIMER {double _dt=-QOP_time();
#define END_TIMER _dt+=QOP_time(); VERB(MED, "%i: %i\n", __LINE__, (int)(1e6*_dt));}
#else
#define BEGIN_TIMER
#define END_TIMER
#endif

#define drotate_vecs(v,r,s) drotate_vecs_func(v,&(r),s)
#define zrotate_vecs(v,r,s) zrotate_vecs_func(v,&(r),s)

static int ngood;

typedef struct {
  dvec td, td1, td2;
  dmat t, t1, ta, ta1, ta2, tb, tb0, tb1, tb2, tb3, y, y0, y1, h;
} eigcg_temps_t;

static eigcg_temps_t *
create_temps(int nev, int m)
{
  // t  : R  m x m  symm
  // t1 : t(0,0 : m-1,m-1)
  // ta : R  m x m
  // ta1: ta(0,0 : m-1,m-1)
  // ta2: ta(0,0 : 2*nev,2*nev)
  // tb : R  m x m
  // tb0: tb(0,0 : m,nev)
  // tb1: tb(0,0 : m-1,m-1)
  // tb2: tb(m-1,0 : 1,nev)
  // tb3: tb(0,0 : m,2*nev)
  // td : R  m
  // td1: td(0 : m-1)
  // td2: td(0 : 2*nev)
  // y  : R  m x 2*nev
  // y0 : y(0,0 : m,nev)
  // y1 : y(0,nev : m,nev)
  // h  : R  2*nev x 2*nev
  eigcg_temps_t *et;
  et = malloc(sizeof(eigcg_temps_t));
  dvec_alloc(et->td, m);
  dmat_alloc(et->t, m, m);
  dmat_alloc(et->ta, m, m);
  dmat_alloc(et->tb, m, m);
  dmat_alloc(et->y, m, 2*nev);
  dmat_alloc(et->h, 2*nev, 2*nev);
  //#if 0
  dsubvec(et->td1, et->td, 0, m-1);
  dsubvec(et->td2, et->td, 0, 2*nev);
  dsubmat(et->t1, et->t, 0, 0, m-1, m-1);
  dsubmat(et->ta1, et->ta, 0, 0, m-1, m-1);
  dsubmat(et->ta2, et->ta, 0, 0, 2*nev, 2*nev);
  dsubmat(et->tb0, et->tb, 0, 0, m, nev);
  dsubmat(et->tb1, et->tb, 0, 0, m-1, m-1);
  dsubmat(et->tb2, et->tb, m-1, 0, 1, nev);
  dsubmat(et->tb3, et->tb, 0, 0, m, 2*nev);
  dsubmat(et->y0, et->y, 0, 0, m, nev);
  dsubmat(et->y1, et->y, 0, nev, m, nev);
  //#endif
  // === end eigCG struct
  return et;
}

static void
free_temps(eigcg_temps_t *et)
{
  dvec_free(et->td);
  dmat_free(et->t);
  dmat_free(et->ta);
  dmat_free(et->tb);
  dmat_free(et->y);
  dmat_free(et->h);
  free(et);
}

static void
drotate_vecs_func(Vector *v[], dmat *r, QDP_Subset subset)
{
  int i, j;
  int n1 = r->size1;
  int n2 = r->size2;
  Vector **tv;

#if 0
  printf("%i %i\n", n1, n2);
  dmat_print(*r);
  for(j=0; j<n1; j++) {
    QLA_Real nrm2;
    r_eq_norm2_V(&nrm2, v[j], subset);
    printf(" %i %g\n", j, nrm2);
  }
#endif
  tv = malloc(n2*sizeof(Vector *));
  for(j=0; j<n2; j++) {
    TRACE;
    create_V(tv[j]);
    TRACE;
    V_eq_zero(tv[j], subset);
    for(i=0; i<n1; i++) {
      QLA_Real s = dmat_get(*r, i, j);
      TRACE;
      V_peq_r_times_V(tv[j], &s, v[i], subset);
      TRACE;
    }
  }
#if 0
  for(j=0; j<n2; j++) {
    QLA_Real nrm2;
    r_eq_norm2_V(&nrm2, tv[j], subset);
    printf(" %i %g\n", j, nrm2);
  }
#endif
  for(j=0; j<n2; j++) {
    V_eq_V(v[j], tv[j], subset);
    destroy_V(tv[j]);
  }
}

static void
zrotate_vecs_func(Vector *v[], zmat *r, QDP_Subset subset)
{
  QLA_Complex s;
  int i, j;
  int n1 = r->size1;
  int n2 = r->size2;
  Vector **tv;

#if 0
  printf("%i %i\n", n1, n2);
  zmat_print(*r);
  for(j=0; j<n1; j++) {
    QLA_Real nrm2;
    r_eq_norm2_V(&nrm2, v[j], subset);
    printf(" %i %g\n", j, nrm2);
  }
#endif
  tv = malloc(n2*sizeof(Vector *));
  for(j=0; j<n2; j++) {
    create_V(tv[j]);
    V_eq_zero(tv[j], subset);
    for(i=0; i<n1; i++) {
      _Complex double cs = zmat_get(*r, i, j);
      QLA_c_eq_r_plus_ir(s, creal(cs), cimag(cs));
      V_peq_c_times_V(tv[j], &s, v[i], subset);
    }
  }
#if 0
  for(j=0; j<n2; j++) {
    QLA_Real nrm2;
    r_eq_norm2_V(&nrm2, tv[j], subset);
    printf(" %i %g\n", j, nrm2);
  }
#endif
  for(j=0; j<n2; j++) {
    V_eq_V(v[j], tv[j], subset);
    destroy_V(tv[j]);
  }
}

static void
normalize(Vector *src[], int n, QDP_Subset subset)
{
  QLA_Real r;
  int i;
  for(i=0; i<n; i++) {
    r_eq_norm2_V(&r, src[i], subset);
    r = 1.0/sqrt(r);
    V_eq_r_times_V(src[i], &r, src[i], subset);
  }
}

static void
orthonormalize(Vector *v[], int i0, int n, QDP_Subset subset)
{
  QLA_Complex z;
  int i, j, nn=i0+n;
  for(i=i0; i<nn; i++) {
    for(j=0; j<i; j++) {
      c_eq_V_dot_V(&z, v[j], v[i], subset);
      V_meq_c_times_V(v[i], &z, v[j], subset);
    }
    normalize(v+i, 1, subset);
  }
}

static void
rayleighRitz(Vector *u[], QLA_Real l[], int nu, int nv, QLA_Real *rsq,
	     QOPPCV(linop_t) *linop, Vector *p, Vector *Mp, QDP_Subset subset)
{
  int i, j, n=nu+nv;
  dvec vv;
  zmat mm, tt;

  if(n>0) {
    VERB(MED, "eigCG: rr: nu %i nv %i\n", nu, nv);
    //printf0("CG: rr: nu %i nv %i\n", nu, nv);
    dvec_alloc(vv, n);
    zmat_alloc(mm, n, n);
    zmat_alloc(tt, n, n);

    BEGIN_TIMER;
    orthonormalize(u, nu, nv, subset);
    END_TIMER;

    BEGIN_TIMER;
    for(i=0; i<nu; i++) {
      _Complex double zz;
      zz = 0.;
      for(j=0; j<nu; j++) {
	zmat_set(mm, i, j, zz);
      }
      zz = l[i];
      zmat_set(mm, i, i, zz);
    }
    //printf("test\n");

    for(j=nu; j<n; j++) {
      _Complex double zz;
      QLA_Complex qz;
      //QLA_Real nrm;
      //r_eq_norm2_V(&nrm, u[j], subset);
      //printf("nrm = %g\n", nrm);
      V_eq_V(p, u[j], subset);
      linop(Mp, p, subset);
      for(i=0; i<j; i++) {
	c_eq_V_dot_V(&qz, u[i], Mp, subset);
	//printf("%i %i %g %g\n", i, j, QLA_real(qz), QLA_imag(qz));
	zz = QLA_real(qz) + I*QLA_imag(qz);
	//printf("%i %i %g %g\n", i, j, creal(zz), cimag(zz));
	zmat_set(mm, i, j, zz);
	zz = conj(zz);
	zmat_set(mm, j, i, zz);
      }
      c_eq_V_dot_V(&qz, u[j], Mp, subset);
      zz = QLA_real(qz) + I*QLA_imag(qz);
      zmat_set(mm, j, j, zz);
    }
    END_TIMER;
    //printf(" mm:\n");
    //gsl_matrix_complex_fprintf(stdout, mm.mat, "  %10g");
    BEGIN_TIMER;
    //printf("zeigs\n");
    //printf("mm\n"); zmat_print(mm);
    zeigs(mm, vv, tt, n);
    //printf("zeigs\n");
    END_TIMER;
    //printf(" tt:\n");
    //gsl_matrix_complex_fprintf(stdout, tt.mat, "  %10g");
    BEGIN_TIMER;
    //printf("zrot\n");
    zrotate_vecs(u, tt, subset);
    //printf("zrot\n");
    END_TIMER;
    if(!rsq) {
      for(i=0; i<n; i++) {
	l[i] = vv.data[i];
      }
    } else {
      Vector *r;
      create_V(r);
      for(i=0; i<n; i++) {
	QLA_Real ev, evi, n2, rn2;
	V_eq_V(p, u[i], subset);
	linop(Mp, p, subset);
	r_eq_re_V_dot_V(&ev, u[i], Mp, subset);
	r_eq_norm2_V(&n2, u[i], subset);
	ev = ev/n2;
	evi = -1/ev;
	V_eq_r_times_V_plus_V(r, &evi, Mp, u[i], subset);
	r_eq_norm2_V(&rn2, r, subset);
	rn2 /= n2;
	//printf0(" %i %-10g %-10g %-10g %-10g\n", i, vv.vec->data[i], ev, rn2, n2);
	l[i] = ev;
	rsq[i] = ev*sqrt(rn2);
      }
      destroy_V(r);
    }
    dvec_free(vv);
    zmat_free(mm);
    zmat_free(tt);
    if(n==1) l[1] = 0;
    VERB(MED, "eigCG: ev[0] %g ev[1] %g ev[%i] %g\n", l[0], l[1], n-1, l[n-1]);
    if(rsq) {
      int ni, ia[6]={0,1,0,0,0,0};
      for(i=0; i<n; i++) if(rsq[i]>1e-4) break;
      ngood = i;
      ni = 2;
      if(i-1>1 && i-1<n-2) ia[ni++] = i-1;
      if(i>1 && i<n-2) ia[ni++] = i;
      ia[ni++] = n-2;
      ia[ni++] = n-1;
      for(i=0; i<ni; i++) printf0("%-12i", ia[i]); printf0("\n");
      for(i=0; i<ni; i++) printf0("%-12g", l[ia[i]]); printf0("\n");
      for(i=0; i<ni; i++) printf0("%-12g", rsq[ia[i]]); printf0("\n");
    }
  }
}

static void
initcg(Vector *x, Vector *u[], QLA_Real l[], int n, Vector *b,
       QDP_Subset subset)
{
  if(n>0) {
    int i;
    V_eq_zero(x, subset);
    for(i=0; i<n; i++) {
      QLA_Complex z;
      c_eq_V_dot_V(&z, u[i], b, subset);
      QLA_c_eq_c_div_r(z, z, l[i]);
      V_peq_c_times_V(x, &z, u[i], subset);
    }
  }
}

static void
initcgm(Vector *x, Vector *u[], dvec *d, int n, Vector *b,
	QDP_Subset subset)
{
  if(n>0) {
    int i;
    V_eq_zero(x, subset);
    for(i=0; i<n; i++) {
      QLA_Complex z;
      QLA_Real l;
      c_eq_V_dot_V(&z, u[i], b, subset);
      l = dvec_get(*d, i);
      QLA_c_eq_c_div_r(z, z, l);
      V_peq_c_times_V(x, &z, u[i], subset);
    }
  }
}

#define gett(ta,v,m,linop,p,Mp,subset) gett_func(&ta,v,m,linop,p,Mp,subset) 
static void
gett_func(dmat *tt, Vector *v[], int n, QOPPCV(linop_t) *linop,
	  Vector *p, Vector *Mp, QDP_Subset subset)
{
  int i, j;
  Vector *sp, *sMp;
  create_V(sp);
  create_V(sMp);
  V_eq_V(sp, p, subset);
  V_eq_V(sMp, Mp, subset);
  for(j=0; j<n; j++) {
    V_eq_V(p, v[j], subset);
    linop(Mp, p, subset);
    for(i=0; i<n; i++) {
      QLA_Complex zz;
      c_eq_V_dot_V(&zz, v[i], Mp, subset);
      dmat_set(*tt, i, j, QLA_real(zz));
      //printf("%i %i %g\n", i, j, QLA_real(zz));
    }
  }
  V_eq_V(p, sp, subset);
  V_eq_V(Mp, sMp, subset);
  destroy_V(sp);
  destroy_V(sMp);
}

#define setv() { \
      QLA_Real s = 1/sqrt(rsq); \
      V_eq_r_times_V(v[nv], &s, r, subset); \
      nv++; \
    }

#if 0
static void
diag_t_old(eigcg_temps_t *et, int nev, int m, Vector *v[], int nv,
	   QDP_Subset subset)
{
  // t  : R  m x m  symm
  // t1 : t(0,0 : m-1,m-1)
  // ta : R  m x m
  // ta1: ta(0,0 : m-1,m-1)
  // ta2: ta(0,0 : 2*nev,2*nev)
  // tb : R  m x m
  // tb0: tb(0,0 : m,nev)
  // tb1: tb(0,0 : m-1,m-1)
  // tb2: tb(m-1,0 : 1,nev)
  // tb3: tb(0,0 : m,2*nev)
  // td : R  m
  // td1: td(0 : m-1)
  // td2: td(0 : 2*nev)
  // y  : R  m x 2*nev
  // y0 : y(0,0 : m,nev)
  // y1 : y(0,nev : m,nev)
  // h  : R  2*nev x 2*nev
#if 0
  dmat 
  dsubvec(et->td1, et->td, 0, m-1);
  dsubvec(et->td2, et->td, 0, 2*nev);
  dsubmat(et->t1, et->t, 0, 0, m-1, m-1);
  dsubmat(et->ta1, et->ta, 0, 0, m-1, m-1);
  dsubmat(et->ta2, et->ta, 0, 0, 2*nev, 2*nev);
  dsubmat(et->tb0, et->tb, 0, 0, m, nev);
  dsubmat(et->tb1, et->tb, 0, 0, m-1, m-1);
  dsubmat(et->tb2, et->tb, m-1, 0, 1, nev);
  dsubmat(et->tb3, et->tb, 0, 0, m, 2*nev);
  dsubmat(et->y0, et->y, 0, 0, m, nev);
  dsubmat(et->y1, et->y, 0, nev, m, nev);
#endif
#if 0
  printf("cmp1:\n");
  gett(ta, v, m, linop, p, Mp, subset);
  compmr(ta, t);
  deigs(ta, td, tb, nev);
  printf("t:\n");
  printvr(td);
#endif
  dmat_copy(et->ta, et->t);
  deigs(et->ta, et->td, et->tb, nev);
  //printf("tm:\n");
  //printvr(td);

  {
    dmat tb2;
    dsubmat(tb2, et->tb, 0, 0, m, 2*nev);
    dmat_copy(et->y, tb2);
  }
#if 0
  dmat_copy(et->y0, et->tb0);
  dmat_copy(et->ta1, et->t1);
  deigs(et->ta1, et->td1, et->tb1, nev);
  dvec_zero(et->tb2);
  dmat_copy(et->y1, et->tb0);
  //mmulanr(h, y, y);
  //printf("y^T y:\n");
  //printmr(h);
  dorthocols(et->y1, et->y0);
#endif
  //mmulanr(h, y, y);
  //printf("y^T y:\n");
  //printmr(h);
  TRACE;
  dmulmnn(et->tb3, et->t, et->y);
  TRACE;
  dmulman(et->h, et->y, et->tb3);
  TRACE;
  deigs(et->h, et->td2, et->ta2, 2*nev);
  //printf("td2:\n");
  //printvr(td2);
  TRACE;
  dmulmnn(et->tb3, et->y, et->ta2);
#if 0
  {
    //mmulnnr(y, t, tb3);
    //mmulanr(h, tb3, y);
    dmulman(h, tb3, tb3);
    printf("q^T y^T t y q:\n");
    printmr(h);
  }
#endif
  TRACE;
  drotate_vecs(v, et->tb3, subset);
  TRACE;
}
#endif

static void
diag_t(eigcg_temps_t *et, int nev, int nv)
{
  dmat t0, ta0, tb0;
  dvec td0;

  dsubmat(t0, et->t, 0, 0, nv, nv);
  dsubmat(ta0, et->ta, 0, 0, nv, nv);
  dsubvec(td0, et->td, 0, nv);
  dsubmat(tb0, et->tb, 0, 0, nv, nv);
  dmat_copy(ta0, t0);
  //printf("ta0\n"); dmat_print(ta0);
  deigs(ta0, td0, tb0, nev);
  //dvec_print(td0);
  //printf("tb0\n"); dmat_print(tb0);
}

static void
set_y(eigcg_temps_t *et, int noff, int nev, int nv, int n)
{
  dmat y0, tb0;
  int nz = n-nv;

  dsubmat(y0, et->y, 0, noff, nv, nev);
  dsubmat(tb0, et->tb, 0, 0, nv, nev);
  dmat_copy(y0, tb0);
  //printf("tb0\n"); dmat_print(tb0);
  //printf("y0\n"); dmat_print(y0);
  if(nz) {
    dmat y1;
    dsubmat(y1, et->y, nv, noff, nz, nev);
    dmat_zero(y1);
  }
  if(noff) {
    dmat y2, y3;
    dsubmat(y2, et->y, 0, 0, n, noff);
    dsubmat(y3, et->y, 0, noff, n, nev);
    //printf("y2\n"); dmat_print(y2);
    //printf("y3\n"); dmat_print(y3);
    dorthocols(y3, y2);
    //printf("y2\n"); dmat_print(y2);
    //printf("y3\n"); dmat_print(y3);

    dmat tb1, t0, y4;
    dsubmat(tb1, et->tb, 0, 0, n, noff+nev);
    dsubmat(t0, et->t, 0, 0, n, n);
    dsubmat(y4, et->y, 0, 0, n, noff+nev);
    dmulmnn(tb1, t0, y4);
    //printf("y4\n"); dmat_print(y4);

    dmat h0;
    dsubmat(h0, et->h, 0, 0, noff+nev, noff+nev);
    dmulman(h0, y4, tb1);

    dmat ta0;
    dvec td0;
    dsubmat(ta0, et->ta, 0, 0, noff+nev, noff+nev);
    dsubvec(td0, et->td, 0, noff+nev);
    deigs(h0, td0, ta0, noff+nev);

    dmulmnn(tb1, y4, ta0);
  }
}

static void
rot_t(eigcg_temps_t *et, int n, Vector *v[], int nv, QDP_Subset subset)
{
  dmat tb0;
  if(n>nv) n = nv;
  dsubmat(tb0, et->tb, 0, 0, nv, n);
  //printf("tb0\n"); dmat_print(tb0);
  drotate_vecs(v, tb0, subset);
}

static void
reset_t(eigcg_temps_t *et, int nev, int m, Vector *v[], int *nvp,
	QLA_Real a0, QLA_Real a, QLA_Real b, Vector *r, QLA_Real rsq,
	QOPPCV(linop_t) *linop, Vector *p, Vector *Mp, QDP_Subset subset)
{
  int i, nv = *nvp;
  Vector *sp, *sMp;
  create_V(sp);
  create_V(sMp);
  V_eq_V(sp, p, subset);
  V_eq_V(sMp, Mp, subset);

  nv = 2*nev;
  setv();
  TRACE;
  dmat_zero(et->t);
  V_eq_V(p, v[nv-1], subset);
  linop(Mp, p, subset);
  for(i=0; i<nv-1; i++) {
    QLA_Real rr;
    r_eq_re_V_dot_V(&rr, v[i], Mp, subset);
    dmat_set(et->t, i, nv-1, rr);
    dmat_set(et->t, nv-1, i, rr);
    dmat_set(et->t, i, i, dvec_get(et->td2,i));
  }
  dmat_set(et->t, nv-1, nv-1, 1./a + b/a0);
  TRACE;
#if 0
  {
    dvec tv1;
    dmat tt1, tt2, tt3;
    dsubvec(tv1, et->td, 0, 2*nev+1);
    dsubmat(tt1, et->t, 0, 0, 2*nev+1, 2*nev+1);
    dsubmat(tt2, et->ta, 0, 0, 2*nev+1, 2*nev+1);
    dsubmat(tt3, et->tb, 0, 0, 2*nev+1, 2*nev+1);
    dmat_copy(tt2, tt1);
    deigs(tt2, tv1, tt3, 2*nev+1);
    printf("t:2nev+1: theory\n");
    printvr(tv1);

    gett(tt2, v, 2*nev+1, linop, p, Mp, subset);
    printf("t:2nev+1: exact\n");
    compmr(tt2, tt1);
    deigs(tt2, tv1, tt3, 2*nev+1);
    printf("t:2nev+1: exact\n");
    printvr(tv1);
  }
#endif
  V_eq_V(p, sp, subset);
  V_eq_V(Mp, sMp, subset);
  destroy_V(sp);
  destroy_V(sMp);
  *nvp = nv;
}


QOP_status_t
QOPPCV(invert_eigcg)(QOPPCV(linop_t) *linop,
		     QOP_invert_arg_t *inv_arg,
		     QOP_resid_arg_t *res_arg,
		     Vector *out,
		     Vector *in,
		     Vector *p,
		     QDP_Subset subset,
		     QOPPCV(eigcg_t) *eigcg
		     vIndexDef)
{
  QLA_Real a, b, a0;
  QLA_Real rsq, relnorm2, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  Vector *r, *Mp;
  int i, iteration=0, total_iterations=0, nrestart=-1;
  int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;

  // === begin eigCG struct
  TRACE;
  Vector **v = NULL;
  Vector **u = eigcg->u;
  QLA_Real *l = eigcg->l;
  QLA_Real *evrsq=NULL;
  //QLA_Real new_low=FLT_MAX;
  int nev = eigcg->nev;
  int m = eigcg->m;
  int numax = eigcg->numax;
  int nu = eigcg->nu;
  int nv = eigcg->nv;
  int addvecs = eigcg->addvecs;
  int nn = eigcg->nn;
  eigcg_temps_t *et;

  TRACE;
  if(u==NULL) {
    eigcg->u = u = malloc(numax*sizeof(Vector *));
    eigcg->l = l = malloc(numax*sizeof(QLA_Real));
    for(i=0; i<numax; i++) {
      create_V(u[i]);
    }
    nu = 0;
    nv = 0;
    nn = 0;
    addvecs = 1;
  }
  TRACE;
  VERB(MED, "eigCG: numax %i m %i nev %i nu %i nn %i nv %i\n", numax, m, nev,nu,nn,nv);

  et = create_temps(nev, m);
  if(QOP_common.verbosity>=QOP_VERB_MED) {
    evrsq = (QLA_Real *)malloc(numax*sizeof(QLA_Real));
  }
  TRACE;

  create_V(r);
  create_V(Mp);

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "eigCG: rsqstop = %g\n", rsqstop);
  rsq = 0;
  relnorm2 = 1.;

  /* Default output values unless reassigned */
  res_arg->final_rsq = 0;
  res_arg->final_rel = 0;
  res_arg->final_iter = 0;
  res_arg->final_restart = 0;

  /* Special case of exactly zero source */
  if(insq == 0.){
    VERB(LOW, "eigCG: exiting because of zero source\n");
    destroy_V(r);
    destroy_V(Mp);
    V_eq_zero(out, subset);
    
    return QOP_SUCCESS;
  }

  //nur = numax - 2*m - nev*max_restarts;
  //if(nur<nev) nur = nev;
  //nn = 0;
  //if(nu+nv<nur+m) { nn = nv; nv = 0; }
#if 0
  if(nu+nn+nv>numax-m) {
    if(nn+nv!=0) {
      BEGIN_TIMER;
      //if(nv>nev) nv=nev;
      nu += nn + nv;
      rayleighRitz(u, l, 0, nu, evrsq, linop, p, Mp, subset);
      //rayleighRitz(u, l, nu, nv, evrsq, linop, p, Mp, subset);
      //if(nv>nev) nv=nev;
      //nu += nv; nv = 0;
      nn = 0;
      nv = 0;
      END_TIMER;
      VERB(MED, "eigCG: ngood %i rsq[0] %g rsq[1] %g rsq[%i] %g\n", ngood,
	   evrsq[0], evrsq[1], nu-1, evrsq[nu-1]);
    }
    addvecs = 0;
  }
#endif
  //nv = 0;
  //if(nu>ngood+nev) nu = ngood+nev;
  //if( nu>=numax-m && ngood>(numax-m)/2 ) { addvecs = 0; }

#if 0
  if(nn>=m) {
    BEGIN_TIMER;
    nu += nn;
    rayleighRitz(u, l, 0, nu, evrsq, linop, p, Mp, subset);
    //rayleighRitz(u, l, nu, nv, evrsq, linop, p, Mp, subset);
    nn = 0;
    END_TIMER;
    VERB(MED, "eigCG: ngood %i rsq[0] %g rsq[1] %g rsq[%i] %g\n", ngood,
	 evrsq[0], evrsq[1], nu-1, evrsq[nu-1]);
    if(nu>numax-2*m) 
  }

  v = u + nu + nn;
#endif

  do {

    BEGIN_TIMER;
    V_eq_V(p, out, subset);
    linop(Mp, p, subset);
    iteration = 1;
    total_iterations++;
    V_eq_V_minus_V(r, in, Mp, subset);
    r_eq_norm2_V(&rsq, r, subset);
    END_TIMER;
    VERB(LOW, "eigCG: (re)start: iter %i rsq = %g rel = %g\n", 
	 total_iterations, rsq, relnorm2);
    if( ((rsqstop <= 0 || rsq<rsqstop) &&
	 (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ||
	(total_iterations>=max_iterations) ) break;
    
    if(addvecs && nn>=m) {
      double oldl, delta;
      int nr;
      nr = numax - 2*m - 1;
      if(nr>=nu) nr = nu-2;
      if(nr<0) nr = 0;
      oldl = l[nr];
      if(nu==0) oldl = 9999;
      BEGIN_TIMER;
      //rayleighRitz(u, l, 0, nu+nn, evrsq, linop, p, Mp, subset);
      rayleighRitz(u, l, nu, nn, evrsq, linop, p, Mp, subset);
      nu += nn;
      nn = 0;
      END_TIMER;
      VERB(MED, "eigCG: ev[0] %g ev[1] %g ev[%i] %g ev[%i] %g\n", l[0], l[1],
	   nr, l[nr], nu-1, l[nu-1]);
      VERB(MED, "eigCG: ngood %i rsq[0] %g rsq[1] %g rsq[%i] %g rsq[%i] %g\n", ngood,
	   evrsq[0], evrsq[1], nr, evrsq[nr], nu-1, evrsq[nu-1]);
      delta = fabs( 2*(oldl-l[nr])/(oldl+l[nr]) );
      VERB(MED, "eigCG: oldl = %g l[%i] = %g  delta = %g\n", oldl, nr, l[nr], delta);
      if(delta<1e-4) addvecs = 0;
      else if(nu>numax-2*m) nu = numax - 2*m;
    }

    if(nu>0) {
      BEGIN_TIMER;
      initcg(p, u, l, nu, r, subset);
      V_peq_V(out, p, subset);
      END_TIMER;
      BEGIN_TIMER;
      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration = 1;
      total_iterations++;
      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);
      
      END_TIMER;
      VERB(LOW, "eigCG: proj u: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);
      if( ((rsqstop <= 0 || rsq<rsqstop) &&
	   (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ||
	  (total_iterations>=max_iterations) ) break;
      //if(nv==0 && nn==0 && nu>nur) nu = nur;
    }

#if 0
    if(nv>0) {
      rayleighRitz(u+nu+nn, l+nu+nn, 0, nv, NULL, linop, p, Mp, subset);

      BEGIN_TIMER;
      initcg(p, u+nu+nn, l+nu+nn, nv, r, subset);
      V_peq_V(out, p, subset);
      END_TIMER;
      BEGIN_TIMER;
      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration = 1;
      total_iterations++;
      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);
      
      END_TIMER;
      VERB(LOW, "eigCG: proj v: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);
      if( (rsq<rsqstop) ||
	  (total_iterations>=max_iterations) ) break;
      if(nv<nev) nn += nv; else nn += nev;
      nv = 0;
    }
    //if(nu+nn>numax-m) nn = numax - nu - m;
#endif

#if 0
    if(nv<2) nv = 0;
    if(nv) {
      //BEGIN_TIMER;
      //if(nv>nev) nv=nev;
      //rayleighRitz(u, l, 0, nu+nv, evrsq, linop, p, Mp, subset);
      //rayleighRitz(u, l, nu, nv, evrsq, linop, p, Mp, subset);
      //if(nv>nev) nv=nev;
      //nu += nv; nv = 0;
      //END_TIMER;
      //VERB(MED, "eigCG: ngood %i rsq[0] %g rsq[1] %g rsq[%i] %g\n", ngood,
      //evrsq[0], evrsq[1], nu-1, evrsq[nu-1]);

      int nnv = nv/2;
      BEGIN_TIMER;
      //gett(et->t, v, nv, linop, p, Mp, subset);
      if(nnv>nev) nnv = nev;
      diag_t(et, nnv, nv);
      set_y(et, 0, nnv, nv, nv);
      diag_t(et, nnv, nv-1);
      set_y(et, nnv, nnv, nv-1, nv);
      END_TIMER;
      BEGIN_TIMER;
      rot_t(et, 2*nnv, v, nv, subset);
      END_TIMER;
      nv = 2*nnv;

      BEGIN_TIMER;
      initcgm(p, v, &et->td, nv, r, subset);
      V_peq_V(out, p, subset);
      END_TIMER;
      BEGIN_TIMER;
      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration = 1;
      total_iterations++;
      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);

      END_TIMER;
      VERB(LOW, "eigCG: proj v: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);
      if( ((rsqstop <= 0 || rsq<rsqstop) &&
	   (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ||
	  (total_iterations>=max_iterations) ) break;
      nn += nv/2;
      nv = 0;
      if(nu+nn>numax-m) nn = numax - nu - m;
    }
#endif

    if(addvecs) {
      v = u + nu + nn;
      dmat_zero(et->t);
      TRACE;
    }

    V_eq_V(p, r, subset);

    while(1) {

      linop(Mp, p, subset);
      iteration++;
      total_iterations++;

      r_eq_re_V_dot_V(&pkp, p, Mp, subset);

      a0 = a;
      a = rsq / pkp;

      VERB(HI, "eigCG: a0 = %g  a = %g  b = %g\n", a0, a, b);
      if(addvecs) {
	if(nv==0) {
	  double td = 1./a;
	  VERB(HI, "eigCG: T %g\n", td);
	  dmat_set(et->t, 0, 0, td);
	  setv();
	} else if(nv<m) {
	  double td = 1./a + b/a0;
	  double ts = -sqrt(b)/a0;
	  VERB(HI, "eigCG: T %g %g\n", td, ts);
	  dmat_set(et->t, nv, nv, td);
	  dmat_set(et->t, nv, nv-1, ts);
	  dmat_set(et->t, nv-1, nv, ts);
	  setv();
	} else {
	  BEGIN_TIMER;
#if 0
	  Vector *sp, *sMp;
	  create_V(sp);
	  create_V(sMp);
	  V_eq_V(sp, p, subset);
	  V_eq_V(sMp, Mp, subset);
	  rayleighRitz(v, l+nu+nn, 0, nv, evrsq, linop, p, Mp, subset);
	  V_eq_V(p, sp, subset);
	  V_eq_V(Mp, sMp, subset);
	  destroy_V(sp);
	  destroy_V(sMp);
	  nv = nev;
#else
	  //orthonormalize(v, 0, nv, subset);
 	  //gett(et->t, v, nv, linop, p, Mp, subset);
	  diag_t(et, nev, nv);
	  set_y(et, 0, nev, nv, nv);
	  diag_t(et, nev, nv-1);
	  set_y(et, nev, nev, nv-1, nv);
	  END_TIMER;
	  BEGIN_TIMER;
	  rot_t(et, 2*nev, v, nv, subset);
	  END_TIMER;
	  BEGIN_TIMER;
	  reset_t(et, nev, m, v, &nv, a0, a, b, r, rsq, linop, p, Mp, subset);
#endif
	  END_TIMER;
	}
	if(QOP_common.verbosity>=QOP_VERB_HI) {
	  dmat ta0, t0, tb0;
	  dvec td0;
	  dsubmat(t0, et->t, 0, 0, nv, nv);
	  dsubmat(ta0, et->ta, 0, 0, nv, nv);
	  dsubmat(tb0, et->tb, 0, 0, nv, nv);
	  dsubvec(td0, et->td, 0, nv);
	  printf0("cmp1:\n");
	  gett(ta0, v, nv, linop, p, Mp, subset);
	  if(QDP_this_node==0) dmat_comp(ta0, t0);
	  deigs(ta0, td0, tb0, nv);
	  printf0("t:\n");
	  //dvec_print(td0);
	  for(i=0; i<nv; i++) {
	    QLA_Real nrm;
	    r_eq_norm2_V(&nrm, v[i], subset);
	    printf("nrm(v[%i]) = %g\n", i, nrm);
	  }
	}
      }

      V_peq_r_times_V(out, &a, p, subset);
      V_meq_r_times_V(r, &a, Mp, subset);
      oldrsq = rsq;
      r_eq_norm2_V(&rsq, r, subset);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);
      
      VERB(HI, "eigCG: iter %i rsq = %g\n", 
	   total_iterations, rsq, relnorm2);

      if( ((rsqstop <= 0 || rsq<rsqstop) &&
	   (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ||
	  (total_iterations>=max_iterations) ) break;

      b = rsq / oldrsq;
      V_eq_r_times_V_plus_V(p, &b, p, r, subset);

    } // end of main CG loop

    if(addvecs) {
      if(nv>nev) {
#if 0
	Vector *sp, *sMp;
	create_V(sp);
	create_V(sMp);
	V_eq_V(sp, p, subset);
	V_eq_V(sMp, Mp, subset);
	rayleighRitz(v, l+nu+nn, 0, nv, evrsq, linop, p, Mp, subset);
	V_eq_V(p, sp, subset);
	V_eq_V(Mp, sMp, subset);
	destroy_V(sp);
	destroy_V(sMp);
#else
	orthonormalize(v, 0, nv, subset);
	gett(et->t, v, nv, linop, p, Mp, subset);
	diag_t(et, nev, nv);
	set_y(et, 0, nev, nv, nv);
	rot_t(et, nev, v, nv, subset);
#endif
	nn += nev;
      }
      nv = 0;
    }

    nrestart++;
  } while( (total_iterations<max_iterations) &&
	   (nrestart<max_restarts) );

  VERB(LOW, "done: nu %i nn %i nv %i\n", nu, nn, nv);
  if(addvecs) {
    rayleighRitz(u, l, 0, nu+nn+nv, evrsq, linop, p, Mp, subset);
    nu += nn + nv;
    nn = 0;
    nv = 0;
  }

#if 0
  if(nv>nev) {
    diag_t(et, nev, nv);
    set_y(et, 0, nev, nv, nv);
    rot_t(et, nev, v, nv, subset);
    nn += nev;
  }
  nv = 0;
#endif

  //if( nuold>=numax-m && ngood>(numax-m)/2 ) { nu = nuold; nv = 0; }

  eigcg->nu = nu;
  eigcg->nv = nv;
  eigcg->nn = nn;
  eigcg->addvecs = addvecs;

  free_temps(et);
  free(evrsq);

  destroy_V(r);
  destroy_V(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

  return QOP_SUCCESS;
}

#endif //HAVE_LAPACK
