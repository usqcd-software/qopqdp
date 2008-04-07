/* eigCG of Stathopoulos and Orginos (arXiv:0707.0131) */
#include <float.h>
//#include "linalg.h"

//#define TRACE printf("%s %s %i\n", __FILE__, __func__, __LINE__);
#define TRACE
#define printf0 if(QDP_this_node==0) printf

#if 0
#define BEGIN_TIMER {double _dt=-QOP_time();
#define END_TIMER _dt+=QOP_time(); printf0("%i: %i\n", __LINE__, (int)(1e6*_dt));}
#else
#define BEGIN_TIMER
#define END_TIMER
#endif

#if 0

static int ngood;

typedef struct {
  vecr td, td1, td2;
  matr t, t1, ta, ta1, ta2, tb, tb0, tb1, tb2, tb3, y, y0, y1, h;
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
  alloc_vec(et->td, m);
  alloc_mat(et->t, m, m);
  alloc_mat(et->ta, m, m);
  alloc_mat(et->tb, m, m);
  alloc_mat(et->y, m, 2*nev);
  alloc_mat(et->h, 2*nev, 2*nev);
  subvec(et->td1, et->td, 0, m-1);
  subvec(et->td2, et->td, 0, 2*nev);
  submat(et->t1, et->t, 0, 0, m-1, m-1);
  submat(et->ta1, et->ta, 0, 0, m-1, m-1);
  submat(et->ta2, et->ta, 0, 0, 2*nev, 2*nev);
  submat(et->tb0, et->tb, 0, 0, m, nev);
  submat(et->tb1, et->tb, 0, 0, m-1, m-1);
  submat(et->tb2, et->tb, m-1, 0, 1, nev);
  submat(et->tb3, et->tb, 0, 0, m, 2*nev);
  submat(et->y0, et->y, 0, 0, m, nev);
  submat(et->y1, et->y, 0, nev, m, nev);
  // === end eigCG struct
  return et;
}

static void
free_temps(eigcg_temps_t *et)
{
  free_vec(et->td);
  free_mat(et->t);
  free_mat(et->ta);
  free_mat(et->tb);
  free_mat(et->y);
  free_mat(et->h);
  free(et);
}

static void
rotate_vecs_func(Vector *v[], gsl_matrix *r, QDP_Subset subset)
{
  int i, j;
  int n1 = r->size1;
  int n2 = r->size2;
  Vector **tv;

  //printf("%i %i\n", n1, n2);
  tv = malloc(n2*sizeof(Vector *));
  for(j=0; j<n2; j++) {
    TRACE;
    create_V(tv[j]);
    TRACE;
    V_eq_zero(tv[j], subset);
    for(i=0; i<n1; i++) {
      QLA_Real s = gsl_matrix_get(r, i, j);
      TRACE;
      V_peq_r_times_V(tv[j], &s, v[i], subset);
      TRACE;
    }
  }
  for(j=0; j<n2; j++) {
    V_eq_V(v[j], tv[j], subset);
    destroy_V(tv[j]);
  }
}

static void
rotate_vecsc_func(Vector **v, gsl_matrix_complex *r, QDP_Subset subset)
{
  int i, j;
  int n1 = r->size1;
  int n2 = r->size2;
  Vector **tv;

  tv = malloc(n2*sizeof(Vector *));
  for(j=0; j<n2; j++) {
    create_V(tv[j]);
    V_eq_zero(tv[j], subset);
    for(i=0; i<n1; i++) {
      gsl_complex gs = gsl_matrix_complex_get(r, i, j);
      QLA_Complex s;
      QLA_c_eq_r_plus_ir(s, GSL_REAL(gs), GSL_IMAG(gs));
      V_peq_c_times_V(tv[j], &s, v[i], subset);
    }
  }
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
  QLA_Complex zz;
  int i, j, n=nu+nv;
  vecr vv;
  matc mm, tt;

  if(n>0) {
    VERB(MED, "eigCG: rr: nu %i nv %i\n", nu, nv);
    //printf0("CG: rr: nu %i nv %i\n", nu, nv);
    alloc_vec(vv, n);
    alloc_matc(mm, n, n);
    alloc_matc(tt, n, n);

    BEGIN_TIMER;
    orthonormalize(u, nu, nv, subset);
    END_TIMER;

    BEGIN_TIMER;
    for(i=0; i<nu; i++) {
      QLA_c_eq_r(zz, 0.);
      for(j=0; j<nu; j++) {
	setc(mm, i, j, zz);
      }
      QLA_c_eq_r(zz, l[i]);
      setc(mm, i, i, zz);
    }

    for(j=nu; j<n; j++) {
      //QLA_Real nrm;
      //r_eq_norm2_V(&nrm, u[j], subset);
      //printf("nrm = %g\n", nrm);
      V_eq_V(p, u[j], subset);
      linop(Mp, p, subset);
      for(i=0; i<j; i++) {
	c_eq_V_dot_V(&zz, u[i], Mp, subset);
	//printf("%i %i %g %g\n", i, j, QLA_real(zz), QLA_imag(zz));
	setc(mm, i, j, zz);
	QLA_c_eq_ca(zz, zz);
	setc(mm, j, i, zz);
      }
      c_eq_V_dot_V(&zz, u[j], Mp, subset);
      setc(mm, j, j, zz);
    }
    END_TIMER;
    //printf(" mm:\n");
    //gsl_matrix_complex_fprintf(stdout, mm.mat, "  %10g");
    BEGIN_TIMER;
    eigsc(mm, vv, tt, n);
    END_TIMER;
    //printf(" tt:\n");
    //gsl_matrix_complex_fprintf(stdout, tt.mat, "  %10g");
    BEGIN_TIMER;
    rotate_vecsc(u, tt, subset);
    END_TIMER;
    if(!rsq) {
      for(i=0; i<n; i++) {
	l[i] = vv.vec->data[i];
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
	rsq[i] = rn2;
      }
      destroy_V(r);
    }
    free_vec(vv);
    free_matc(mm);
    free_matc(tt);
    if(n==1) l[1] = 0;
    VERB(MED, "CG: ev[0] %g ev[1] %g ev[%i] %g\n", l[0], l[1], n-1, l[n-1]);
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
gett_func(gsl_matrix *tt, Vector *v[], int n, QOPPCV(linop_t) *linop,
	  Vector *p, Vector *Mp, QDP_Subset subset)
{
  int i, j;
  for(j=0; j<n; j++) {
    V_eq_V(p, v[j], subset);
    linop(Mp, p, subset);
    for(i=0; i<n; i++) {
      QLA_Complex zz;
      c_eq_V_dot_V(&zz, v[i], Mp, subset);
      gsl_matrix_set(tt, i, j, QLA_real(zz));
      //printf("%i %i %g\n", i, j, QLA_real(zz));
    }
  }
}

#define setv() { \
      QLA_Real s = 1/sqrt(rsq); \
      V_eq_r_times_V(v[nv], &s, r, subset); \
      nv++; \
    }

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
  matr 
  subvec(et->td1, et->td, 0, m-1);
  subvec(et->td2, et->td, 0, 2*nev);
  submat(et->t1, et->t, 0, 0, m-1, m-1);
  submat(et->ta1, et->ta, 0, 0, m-1, m-1);
  submat(et->ta2, et->ta, 0, 0, 2*nev, 2*nev);
  submat(et->tb0, et->tb, 0, 0, m, nev);
  submat(et->tb1, et->tb, 0, 0, m-1, m-1);
  submat(et->tb2, et->tb, m-1, 0, 1, nev);
  submat(et->tb3, et->tb, 0, 0, m, 2*nev);
  submat(et->y0, et->y, 0, 0, m, nev);
  submat(et->y1, et->y, 0, nev, m, nev);
#endif
#if 0
  printf("cmp1:\n");
  gett(ta, v, m, linop, p, Mp, subset);
  compmr(ta, t);
  eigsr(ta, td, tb, nev);
  printf("t:\n");
  printvr(td);
#endif
  copyr(et->ta, et->t);
  eigsr(et->ta, et->td, et->tb, nev);
  //printf("tm:\n");
  //printvr(td);

  {
    matr tb2;
    submat(tb2, et->tb, 0, 0, m, 2*nev);
    copyr(et->y, tb2);
  }
#if 0
  copyr(et->y0, et->tb0);
  copyr(et->ta1, et->t1);
  eigsr(et->ta1, et->td1, et->tb1, nev);
  zeror(et->tb2);
  copyr(et->y1, et->tb0);
  //mmultanr(h, y, y);
  //printf("y^T y:\n");
  //printmr(h);
  orthocolsr(et->y1, et->y0);
#endif
  //mmultanr(h, y, y);
  //printf("y^T y:\n");
  //printmr(h);
  TRACE;
  mmultnnr(et->tb3, et->t, et->y);
  TRACE;
  mmultanr(et->h, et->y, et->tb3);
  TRACE;
  eigsr(et->h, et->td2, et->ta2, 2*nev);
  //printf("td2:\n");
  //printvr(td2);
  TRACE;
  mmultnnr(et->tb3, et->y, et->ta2);
#if 0
  {
    //mmultnnr(y, t, tb3);
    //mmultanr(h, tb3, y);
    mmultanr(h, tb3, tb3);
    printf("q^T y^T t y q:\n");
    printmr(h);
  }
#endif
  TRACE;
  rotate_vecs(v, et->tb3, subset);
  TRACE;
}

static void
diag_t(eigcg_temps_t *et, int nev, int m, Vector *v[], int nv,
       QDP_Subset subset)
{
  matr t0, ta0, tb0;
  vecr td0;

  submat(t0, et->t, 0, 0, nv, nv);
  submat(ta0, et->ta, 0, 0, nv, nv);
  subvec(td0, et->td, 0, nv);
  submat(tb0, et->tb, 0, 0, nv, nv);
  copyr(ta0, t0);
  eigsr(ta0, td0, tb0, nv);
}

static void
rot_t(eigcg_temps_t *et, int nev, int m, Vector *v[], int nv,
      QDP_Subset subset)
{
  matr tb0;
  int n=2*nev;
  if(n>nv) n = nv;
  submat(tb0, et->tb, 0, 0, nv, n);
  rotate_vecs(v, tb0, subset);
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
  zeror(et->t);
  V_eq_V(p, v[nv-1], subset);
  linop(Mp, p, subset);
  for(i=0; i<nv-1; i++) {
    QLA_Real rr;
    r_eq_re_V_dot_V(&rr, v[i], Mp, subset);
    setr(et->t, i, nv-1, rr);
    setr(et->t, nv-1, i, rr);
    setr(et->t, i, i, getvr(et->td2,i));
  }
  setr(et->t, nv-1, nv-1, 1./a + b/a0);
  TRACE;
#if 0
  {
    vecr tv1;
    matr tt1, tt2, tt3;
    subvec(tv1, td, 0, 2*nev+1);
    submat(tt1, t, 0, 0, 2*nev+1, 2*nev+1);
    submat(tt2, ta, 0, 0, 2*nev+1, 2*nev+1);
    submat(tt3, tb, 0, 0, 2*nev+1, 2*nev+1);
    copyr(tt2, tt1);
    eigsr(tt2, tv1, tt3, 2*nev+1);
    printf("t:2nev+1: theory\n");
    printvr(tv1);

    gett(tt2, v, 2*nev+1, linop, p, Mp, subset);
    printf("t:2nev+1: exact\n");
    compmr(tt2, tt1);
    eigsr(tt2, tv1, tt3, 2*nev+1);
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

#endif

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
#if 0
  QLA_Real a, b, a0;
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  Vector *r, *Mp;
  int i, iteration=0, total_iterations=0, nrestart=-1;
  int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;
  int addvecs = 1;

  // === begin eigCG struct
  TRACE;
  Vector **v = NULL;
  Vector **u = eigcg->u;
  QLA_Real *l = eigcg->l;
  QLA_Real *evrsq;
  QLA_Real new_low=FLT_MAX;
  int nev = eigcg->nev;
  int m = eigcg->m;
  int numax = eigcg->numax;
  int nu = eigcg->nu;
  int nv = eigcg->nv;
  int nuold;
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
  }
  TRACE;
  VERB(MED, "CG: numax %i m %i nev %i nu %i nv %i\n", numax, m, nev, nu, nv);

  et = create_temps(nev, m);
  evrsq = (QLA_Real *)malloc(numax*sizeof(QLA_Real));
  TRACE;

  create_V(r);
  create_V(Mp);

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "CG: rsqstop = %g\n", rsqstop);
  rsq = 0;
  oldrsq = rsq;

  if(nv) {
    rayleighRitz(u, l, 0, nu+nv, evrsq, linop, p, Mp, subset);
    nu += nv; nv = 0;
  }
  //if(nu>ngood+nev) nu = ngood+nev;
  if( nuold>=numax-m && ngood>(numax-m)/2 ) { addvecs = 0; }
#if 0
  if(nv>nev) {
    BEGIN_TIMER;
    rayleighRitz(u+nu, l+nu, 0, nv, linop, p, Mp, subset);
    nv = nev;
    END_TIMER;
  }
  if(nv==nev) {
    BEGIN_TIMER;
    //nv = 1;
    rayleighRitz(u, l, nu, nv, linop, p, Mp, subset);
    nu += nv; nv = 0;
    END_TIMER;
  }
#endif
  if(nv) {
    BEGIN_TIMER;
    rayleighRitz(u, l, nu, nv, NULL, linop, p, Mp, subset);
    nu += nv; nv = 0;
    END_TIMER;
  }
  BEGIN_TIMER;
  V_eq_V(p, out, subset);
  linop(Mp, p, subset);
  V_eq_V_minus_V(r, in, Mp, subset);
  initcg(p, u, l, nu, r, subset);
  V_peq_V(out, p, subset);
  END_TIMER;
  if(nu>numax-m) nu = numax - m;
  nuold = nu;
  while(1) {

    if( (total_iterations==0) ||
	(iteration>=restart_iterations) ||
	(total_iterations>=max_iterations) ||
	(rsq<rsqstop) ) {  /* only way out */

      if( (total_iterations>=max_iterations) ||
	  (nrestart>=max_restarts) ) break;
      nrestart++;

      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration = 1;
      total_iterations++;

      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);
      VERB(LOW, "CG: (re)start: iter %i rsq = %g\n", total_iterations, rsq);
      if( (rsq<rsqstop) ||
	  (total_iterations>=max_iterations) ) break;

#if 0
      if(nv>nev) {
	BEGIN_TIMER;
	rayleighRitz(v, l+nu, 0, nv, NULL, linop, p, Mp, subset);
	nv = nev;
	END_TIMER;
      }
      if(nv==nev) {
	BEGIN_TIMER;
	//nv = 1;
	rayleighRitz(u, l, nu, nv, NULL, linop, p, Mp, subset);
	nu += nv; nv = 0;
	END_TIMER;
      }
#endif
      if(addvecs) {
      printf0(" %i %i %i %i\n", total_iterations, nuold, nu, nv);
      if(nv) {
	int nuadd=nv;
	//if(nuadd>nev) nuadd = nev;
	diag_t(et, nev, m, v, nv, subset);
	printf0("ev: %g %g %g %g\n", l[0], l[nuold-1], new_low, getvr(et->td,0));
	if(getvr(et->td,0)<new_low) new_low = getvr(et->td,0);
	rot_t(et, nev, m, v, nv, subset);
	nu += nuadd;
	nv = 0;
      }
      if(nu>numax-m) {
	if(new_low<l[nuold-1]) {
	  //int nl = numax - m - nev;
	  int nl = 0;
	  rayleighRitz(u+nl, l+nl, 0, nu-nl, NULL, linop, p, Mp, subset);
	  if(nuold>nl) nuold = nl;
	}
	nu = numax - m;
	new_low = FLT_MAX;
      }
      //if(nu+m>numax) addvecs = 0;
      v = u + nu;
      //nv = 0;
      zeror(et->t);
      TRACE;
      }

      V_eq_V(p, r, subset);

    } else {

      //r_eq_re_V_dot_V(&b, Mp, r, subset);
      //b = -a*b/oldrsq;
      b = rsq / oldrsq;
      V_eq_r_times_V_plus_V(p, &b, p, r, subset);

    }

    linop(Mp, p, subset);
    iteration++;
    total_iterations++;

    r_eq_re_V_dot_V(&pkp, p, Mp, subset);

    a0 = a;
    a = rsq / pkp;

    if(addvecs) {
      if(nv==0) {
	setr(et->t, 0, 0, 1./a);
	setv();
      } else if(nv<m) {
	QLA_Real ts;
	setr(et->t, nv, nv, 1./a + b/a0);
	ts = -sqrt(b)/a0;
	setr(et->t, nv, nv-1, ts);
	setr(et->t, nv-1, nv, ts);
	setv();
      } else {
	BEGIN_TIMER;
	diag_t(et, nev, m, v, nv, subset);
	END_TIMER;
	BEGIN_TIMER;
	rot_t(et, nev, m, v, nv, subset);
	END_TIMER;
	BEGIN_TIMER;
	reset_t(et, nev, m, v, &nv, a0, a, b, r, rsq, linop, p, Mp, subset);
	END_TIMER;
      }
    }
#if 0
    for(i=0; i<nv; i++) {
      QLA_Real nrm;
      r_eq_norm2_V(&nrm, v[i], subset);
      printf("nrm(v[%i]) = %g\n", i, nrm);
    }
#endif

    V_peq_r_times_V(out, &a, p, subset);
    V_meq_r_times_V(r, &a, Mp, subset);
    oldrsq = rsq;
    r_eq_norm2_V(&rsq, r, subset);
    VERB(HI, "CG: iter %i rsq = %g\n", total_iterations, rsq);
  }
  //rayleighRitz(u, l, nu, nv, linop, p, Mp, subset);
  VERB(LOW, "done: nu %i nv %i\n", nu, nv);

  //if( nuold>=numax-m && ngood>(numax-m)/2 ) { nu = nuold; nv = 0; }

  eigcg->nu = nu;
  eigcg->nv = nv;

  free_temps(et);
  free(evrsq);

  destroy_V(r);
  destroy_V(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

#endif
  return QOP_SUCCESS;
}
