/* eigCG of Stathopoulos and Orginos (arXiv:0707.0131) */
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

//#define TRACE printf("%s %s %i\n", __FILE__, __func__, __LINE__);
#define TRACE

typedef struct {
  gsl_vector *vec;
  gsl_vector_view view;
} vecr;

typedef struct {
  gsl_matrix *mat;
  gsl_matrix_view view;
} matr;

#define alloc_vec(v,n) { (v).vec = gsl_vector_alloc(n); }
#define alloc_mat(m,r,c) { (m).mat = gsl_matrix_alloc(r,c); }
#define free_vec(v) gsl_vector_free((v).vec)
#define free_mat(m) gsl_matrix_free((m).mat)
#define subvec(v1,v2,k,n) { (v1).view=gsl_vector_subvector((v2).vec,k,n); (v1).vec=&(v1).view.vector; }
#define submat(m1,m2,k1,k2,n1,n2) { (m1).view=gsl_matrix_submatrix((m2).mat,k1,k2,n1,n2); (m1).mat=&(m1).view.matrix; }
#define copyr(m1,m2) gsl_matrix_memcpy((m1).mat,(m2).mat)
#define zeror(m) gsl_matrix_set_zero((m).mat)
#define getvr(v,k) gsl_vector_get((v).vec,k)
#define setr(m,k1,k2,r) gsl_matrix_set((m).mat,k1,k2,r)
#define eigsr(m,d,v,n) {gsl_eigen_symmv((m).mat,(d).vec,(v).mat,ewrk((m).mat->size1)); gsl_eigen_symmv_sort((d).vec,(v).mat,GSL_EIGEN_SORT_VAL_ASC);}
#define orthocolsr(y1,y0) orthocolsr_func(y1.mat,y0.mat)
#define mmultnnr(a,b,c) gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,(b).mat,(c).mat,0.0,(a).mat)
#define mmultanr(a,b,c) gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,(b).mat,(c).mat,0.0,(a).mat)
#define rotate_vecs(v,r,s) rotate_vecs_func(v,(r).mat,s)
#define printvr(v) gsl_vector_fprintf(stdout, v.vec, "  %-12g")
#define printmr(m) printmr_func(m.mat)
#define gett(ta,v,m,linop,p,Mp,subset) gett_func(ta.mat,v,m,linop,p,Mp,subset) 
#define compmr(a,b) compmr_func(a.mat,b.mat)

typedef struct {
  gsl_vector_complex *vec;
  gsl_vector_complex_view view;
} vecc;

typedef struct {
  gsl_matrix_complex *mat;
  gsl_matrix_complex_view view;
} matc;

#define alloc_vecc(v,n) { (v).vec = gsl_vector_complex_alloc(n); }
#define alloc_matc(m,r,c) { (m).mat = gsl_matrix_complex_alloc(r,c); }
#define free_vecc(v) gsl_vector_complex_free((v).vec)
#define free_matc(m) gsl_matrix_complex_free((m).mat)
#define setc(m,k1,k2,r) {gsl_complex zt; GSL_SET_COMPLEX(&zt,QLA_real(r),QLA_imag(r));gsl_matrix_complex_set((m).mat,k1,k2,zt);}
#define eigsc(m,d,v,n) {gsl_eigen_hermv((m).mat,(d).vec,(v).mat,ewrkc((m).mat->size1)); gsl_eigen_hermv_sort((d).vec,(v).mat,GSL_EIGEN_SORT_VAL_ASC);}
#define rotate_vecsc(v,r,s) rotate_vecsc_func(v,(r).mat,s)

void
printmr_func(gsl_matrix *m)
{
  int i, j, n1, n2;
  n1 = m->size1;
  n2 = m->size2;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {
      double mv = gsl_matrix_get(m, i, j);
      printf("%12g ", mv);
    }
    printf("\n");
  }
}

void
compmr_func(gsl_matrix *a, gsl_matrix *b)
{
  int i, j, n1, n2;
  n1 = a->size1;
  n2 = a->size2;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {
      double av = gsl_matrix_get(a, i, j);
      double bv = gsl_matrix_get(b, i, j);
      printf("%12g %12g %12g\n", av, bv, 0.5*fabs(av-bv)/(av+bv));
    }
  }
}

void
orthocolsr_func(gsl_matrix *v, gsl_matrix *u)
{
  int i, j, ni, nj;
  ni = v->size2;
  nj = u->size2;
  for(i=0; i<ni; i++) {
    gsl_vector_view vv = gsl_matrix_column(v, i);
    for(j=0; j<nj; j++) {
      double s;
      gsl_vector_view uu = gsl_matrix_column(u, j);
      gsl_blas_ddot(&uu.vector, &vv.vector, &s);
      s = -s;
      gsl_blas_daxpy(s, &uu.vector, &vv.vector);
    }
    for(j=0; j<i; j++) {
      double s;
      gsl_vector_view uu = gsl_matrix_column(v, j);
      gsl_blas_ddot(&uu.vector, &vv.vector, &s);
      s = -s;
      gsl_blas_daxpy(s, &uu.vector, &vv.vector);
    }
    double f;
    f = 1./gsl_blas_dnrm2(&vv.vector);
    gsl_blas_dscal(f, &vv.vector);
  }
}

gsl_eigen_symmv_workspace *
ewrk(int n)
{
  static gsl_eigen_symmv_workspace *w=NULL;
  static int s=0;
  if(n>s) {
    if(w) gsl_eigen_symmv_free(w);
    w = gsl_eigen_symmv_alloc(n);
  }
  return w;
}

gsl_eigen_hermv_workspace *
ewrkc(int n)
{
  static gsl_eigen_hermv_workspace *w=NULL;
  static int s=0;
  if(n>s) {
    if(w) gsl_eigen_hermv_free(w);
    w = gsl_eigen_hermv_alloc(n);
  }
  return w;
}

void
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

void
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
rayleighRitz(Vector *u[], QLA_Real l[], int nu, int nv, QOPPCV(linop_t) *linop,
	     Vector *p, Vector *Mp, QDP_Subset subset)
{
  QLA_Complex zz;
  int i, j, n=nu+nv;
  vecr vv;
  matc mm, tt;

  if(n>0) {
    VERB(MED, "CG: rr: nu %i nv %i\n", nu, nv);
    //printf("rr: %i\n", n);
    alloc_vec(vv, n);
    alloc_matc(mm, n, n);
    alloc_matc(tt, n, n);
    orthonormalize(u, 0, n, subset);
    //orthonormalize(u, nu, nv, subset);
    for(j=0; j<n; j++) {
      QLA_Real nrm;
      r_eq_norm2_V(&nrm, u[j], subset);
      //printf("nrm = %g\n", nrm);
      V_eq_V(p, u[j], subset);
      linop(Mp, p, subset);
      for(i=j; i<n; i++) {
	c_eq_V_dot_V(&zz, u[i], Mp, subset);
	//printf("%i %i %g %g\n", i, j, QLA_real(zz), QLA_imag(zz));
	setc(mm, i, j, zz);
	QLA_c_eq_ca(zz, zz);
	setc(mm, j, i, zz);
      }
    }
    //printf(" mm:\n");
    //gsl_matrix_complex_fprintf(stdout, mm.mat, "  %10g");
    eigsc(mm, vv, tt, n);
    //printf(" tt:\n");
    //gsl_matrix_complex_fprintf(stdout, tt.mat, "  %10g");
    rotate_vecsc(u, tt, subset);
    for(i=0; i<n; i++) {
      l[i] = vv.vec->data[i];
#if 0
      QLA_Real ev, evi, n2, rn2;
      Vector *r;
      create_V(r);
      V_eq_V(p, u[i], subset);
      linop(Mp, p, subset);
      r_eq_re_V_dot_V(&ev, u[i], Mp, subset);
      r_eq_norm2_V(&n2, u[i], subset);
      ev = ev/n2;
      evi = -1/ev;
      V_eq_r_times_V_plus_V(r, &evi, Mp, u[i], subset);
      r_eq_norm2_V(&rn2, r, subset);
      printf(" %i %-10g %-10g %-10g %-10g\n", i, l[i], ev, rn2, n2);
      destroy_V(r);
#endif
    }
    free_vec(vv);
    free_matc(mm);
    free_matc(tt);
    if(n==1) l[1] = 0;
    VERB(MED, "CG: ev[0] %g ev[1] %g ev[%i] %g\n", l[0], l[1], n-1, l[n-1]);
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

void
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
  QLA_Real rsq, oldrsq, pkp;
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
  int nev = eigcg->nev;
  int m = eigcg->m;
  int numax = eigcg->numax;
  int nu = eigcg->nu;
  int nv = eigcg->nv;

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

  //fprintf(stderr, "numax %i m %i nev %i\n", numax, m, nev);
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
  vecr td, td1, td2;
  matr t, t1, ta, ta1, ta2, tb, tb0, tb1, tb2, tb3, y, y0, y1, h;
  alloc_vec(td, m);
  alloc_mat(t, m, m);
  alloc_mat(ta, m, m);
  alloc_mat(tb, m, m);
  alloc_mat(y, m, 2*nev);
  alloc_mat(h, 2*nev, 2*nev);
  subvec(td1, td, 0, m-1);
  subvec(td2, td, 0, 2*nev);
  submat(t1, t, 0, 0, m-1, m-1);
  submat(ta1, ta, 0, 0, m-1, m-1);
  submat(ta2, ta, 0, 0, 2*nev, 2*nev);
  submat(tb0, tb, 0, 0, m, nev);
  submat(tb1, tb, 0, 0, m-1, m-1);
  submat(tb2, tb, m-1, 0, 1, nev);
  submat(tb3, tb, 0, 0, m, 2*nev);
  submat(y0, y, 0, 0, m, nev);
  submat(y1, y, 0, nev, m, nev);
  // === end eigCG struct
  TRACE;

  create_V(r);
  create_V(Mp);

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "CG: rsqstop = %g\n", rsqstop);
  rsq = 0;
  oldrsq = rsq;

  if(nv>nev) {
    rayleighRitz(u+nu, l+nu, 0, nv, linop, p, Mp, subset);
    nv = nev;
  }
  if(nv==nev) {
    //nv = 1;
    rayleighRitz(u, l, nu, nv, linop, p, Mp, subset);
    nu += nv; nv = 0;
  }
  V_eq_V(p, out, subset);
  linop(Mp, p, subset);
  V_eq_V_minus_V(r, in, Mp, subset);
  initcg(p, u, l, nu, r, subset);
  V_peq_V(out, p, subset);
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

      if(nv>nev) {
	rayleighRitz(v, l+nu, 0, nv, linop, p, Mp, subset);
	nv = nev;
      }
      if(nv==nev) {
	//nv = 1;
	rayleighRitz(u, l, nu, nv, linop, p, Mp, subset);
	nu += nv; nv = 0;
      }
      if(nu+m>numax) nu = numax - m;
      v = u + nu;
      nv = 0;
      zeror(t);
      TRACE;

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

#define setv() { \
      QLA_Real s = 1/sqrt(rsq); \
      V_eq_r_times_V(v[nv], &s, r, subset); \
      nv++; \
    }
    if(nv==0) {
      setr(t, 0, 0, 1./a);
      setv();
    } else if(nv<m) {
      QLA_Real ts;
      setr(t, nv, nv, 1./a + b/a0);
      ts = -sqrt(b)/a0;
      setr(t, nv, nv-1, ts);
      setr(t, nv-1, nv, ts);
      setv();
    } else {
      Vector *sp, *sMp;
      create_V(sp);
      create_V(sMp);
      V_eq_V(sp, p, subset);
      V_eq_V(sMp, Mp, subset);
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
      printf("cmp1:\n");
      gett(ta, v, m, linop, p, Mp, subset);
      compmr(ta, t);
      eigsr(ta, td, tb, nev);
      printf("t:\n");
      printvr(td);
#endif
      copyr(ta, t);
      eigsr(ta, td, tb, nev);
      //printf("tm:\n");
      //printvr(td);
      copyr(y0, tb0);
      copyr(ta1, t1);
      eigsr(ta1, td1, tb1, nev);
      zeror(tb2);
      copyr(y1, tb0);
      //mmultanr(h, y, y);
      //printf("y^T y:\n");
      //printmr(h);
      orthocolsr(y1, y0);
      //mmultanr(h, y, y);
      //printf("y^T y:\n");
      //printmr(h);
      TRACE;
      mmultnnr(tb3, t, y);
      TRACE;
      mmultanr(h, y, tb3);
      TRACE;
      eigsr(h, td2, ta2, 2*nev);
      //printf("td2:\n");
      //printvr(td2);
      TRACE;
      mmultnnr(tb3, y, ta2);
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
      rotate_vecs(v, tb3, subset);
      TRACE;

      nv = 2*nev;
      setv();
      TRACE;
      zeror(t);
      V_eq_V(p, v[nv-1], subset);
      linop(Mp, p, subset);
      for(i=0; i<nv-1; i++) {
	QLA_Real rr;
	r_eq_re_V_dot_V(&rr, v[i], Mp, subset);
	setr(t, i, nv-1, rr);
	setr(t, nv-1, i, rr);
	setr(t, i, i, getvr(td2,i));
      }
      setr(t, nv-1, nv-1, 1./a + b/a0);
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
  eigcg->nu = nu;
  eigcg->nv = nv;

#if 0
  for(i=0; i<numax; i++) {
    destroy_V(u[i]);
  }
  free(u);
  free(l);
#endif

  destroy_V(r);
  destroy_V(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

  return QOP_SUCCESS;
}
