#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

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

static void
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

static void
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

static void
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

static gsl_eigen_symmv_workspace *
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

static gsl_eigen_hermv_workspace *
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
