#include <stddef.h>
#include <math.h>

//#define sscal sscal_
#define dscal dscal_
#if 1
#define sscal sscal_
#define saxpy saxpy_
#define dsdot dsdot_
#define ccopy ccopy_
#define scnrm2 scnrm2_
#define dscal dscal_
#define daxpy daxpy_
#define ddot ddot_
#define zcopy zcopy_
#define dznrm2 dznrm2_
#define dgemm dgemm_
#define dsyev dsyev_
#define zheev zheev_
#endif

extern void sscal(int *, float *, float *, int *);
extern void saxpy(int *, float *, float *, int *, float *, int *);
extern double dsdot(int *, float *, int *, float *, int *);
extern void ccopy(int *, QLA_F_Complex *, int *, QLA_F_Complex *, int *);
extern float scnrm2(int *, QLA_F_Complex *, int *);

extern void dscal(int *, double *, double *, int *);
extern void daxpy(int *, double *, double *, int *, double *, int *);
extern double ddot(int *, double *, int *, double *, int *);
extern void zcopy(int *, QLA_D_Complex *, int *, QLA_D_Complex *, int *);
extern double dznrm2(int *, QLA_D_Complex *, int *);


extern int dgemm(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
		 double *a, int *lda, double *b, int *ldb, double *beta, double *c,
		 int *ldc);
extern int dsyev(char *jobz, char *uplo, int *n, double *a,
		 int *lda, double *w, double *work,
		 int *lwork, int *info);
extern int zheev(char *jobz, char *uplo, int *n, _Complex double *a,
		 int *lda, double *w, _Complex double *work,
		 int *lwork, double *rwork, int *info);


#if QOP_Precision == 1

#define v_eq_zero(r, n) {int ii,n2=2*(n); for(ii=0; ii<n2; ii++) ((float *)(r))[ii] = 0.0; }
#define v_eq_v(r, a, n) {int one=1; ccopy(&n, a, &one, r, &one);}
#define v_peq_v(r, a, n) {float p_one=1.0; int one=1,n2=2*(n); saxpy(&n2, &p_one, (float *)a, &one, (float *)r, &one);}
#define v_meq_v(r, a, n) {float m_one=-1.0; int one=1,n2=2*(n); saxpy(&n2, &m_one, (float *)a, &one, (float *)r, &one);}
#define v_eq_v_minus_v(r, a, b, n) { v_eq_v(r, a, n); v_meq_v(r, b, n); }
#define v_teq_r(r, s, n) {float ss=s; int one=1,n2=2*(n); sscal(&n2, &ss, (float *)r, &one);}
#define v_peq_r_times_v(r, s, a, n) {float ss=s; int one=1,n2=2*(n); saxpy(&n2, &ss, (float *)a, &one, (float *)r, &one);}
#define v_meq_r_times_v(r, s, a, n) {float ms=-s; int one=1,n2=2*(n); saxpy(&n2, &ms, (float *)a, &one, (float *)r, &one);}
//#define norm2_v(a, n) ({double lnrm2; int one=1; lnrm2 = scnrm2(&n, a, &one); lnrm2 *= lnrm2; QMP_sum_double(&lnrm2); lnrm2;})
#define norm2_v(a, n) ({double lnrm2; int one=1,n2=2*(n); lnrm2 = dsdot(&n2, (float *)a, &one, (float *)a, &one); QMP_sum_double(&lnrm2); lnrm2;})
#define re_v_dot_v(a, b, n) ({double ldot; int one=1,n2=2*(n); ldot = dsdot(&n2, (float *)a, &one, (float *)b, &one); QMP_sum_double(&ldot); ldot;})

#else

#define v_eq_zero(r, n) {int ii,n2=2*(n); for(ii=0; ii<n2; ii++) ((double *)(r))[ii] = 0.0; }
#define v_eq_v(r, a, n) {int one=1; zcopy(&n, a, &one, r, &one);}
#define v_peq_v(r, a, n) {double p_one=1.0; int one=1,n2=2*(n); daxpy(&n2, &p_one, (double *)a, &one, (double *)r, &one);}
#define v_meq_v(r, a, n) {double m_one=-1.0; int one=1,n2=2*(n); daxpy(&n2, &m_one, (double *)a, &one, (double *)r, &one);}
#define v_eq_v_minus_v(r, a, b, n) { v_eq_v(r, a, n); v_meq_v(r, b, n); }
#define v_teq_r(r, s, n) {double ss=s; int one=1,n2=2*(n); dscal(&n2, &ss, (double *)r, &one);}
#define v_peq_r_times_v(r, s, a, n) {double ss=s; int one=1,n2=2*(n); daxpy(&n2, &ss, (double *)a, &one, (double *)r, &one);}
#define v_meq_r_times_v(r, s, a, n) {double ms=-s; int one=1,n2=2*(n); daxpy(&n2, &ms, (double *)a, &one, (double *)r, &one);}
//#define norm2_v(a, n) ({double lnrm2; int one=1; lnrm2 = dznrm2(&n, a, &one); lnrm2 *= lnrm2; QMP_sum_double(&lnrm2); lnrm2;})
#define norm2_v(a, n) ({double lnrm2; int one=1,n2=2*(n); lnrm2 = ddot(&n2, (double *)a, &one, (double *)a, &one); QMP_sum_double(&lnrm2); lnrm2;})
#define re_v_dot_v(a, b, n) ({double ldot; int one=1,n2=2*(n); ldot = ddot(&n2, (double *)a, &one, (double *)b, &one); QMP_sum_double(&ldot); ldot;})

#endif

#define v_eq_r_times_v(r, s, a, n) { v_eq_v(r,a,n); v_teq_r(r,s,n); }




typedef struct {
  int size;
  int stride;
  double *data;
} dvec;

typedef struct {
  int size1; // nrows
  int size2; // ncols
  int tda;   // row stride
  double *data;
} dmat;

typedef struct {
  int size;
  int stride;
  _Complex double *data;
} zvec;

typedef struct {
  int size1; // nrows
  int size2; // ncols
  int tda;   // row stride
  _Complex double *data;
} zmat;

#define zmat_alloc(m,r,c) { (m).size1 = r; (m).size2 = c; (m).tda = (m).size2; (m).data = (_Complex double *)malloc((m).size1*(m).size2*sizeof(_Complex double)); }
#define zmat_free(m) free((m).data)
#define zmat_get(m,k1,k2) (m).data[(m).tda*(k1)+(k2)]
#define zmat_set(m,k1,k2,z) (m).data[(m).tda*(k1)+(k2)] = (z)
#define zmat_copy(m1,m2) { for(int _i=0; _i<(m1).size1; _i++) for(int _j=0; _j<(m1).size2; _j++) dmat_set(m1,_i,_j,dmat_get(m2,_i,_j)); }
#define zmat_copya(m1,m2) { for(int _i=0; _i<(m1).size1; _i++) for(int _j=0; _j<(m1).size2; _j++) dmat_set(m1,_i,_j,conj(dmat_get(m2,_j,_i))); }
#define zmat_print(m) {int _i,_j; for(_i=0; _i<(m).size1; _i++) for(_j=0; _j<(m).size2; _j++) printf("  %i %i %-12g %-12g\n", _i, _j, creal(zmat_get(m,_i,_j)), cimag(zmat_get(m,_i,_j)));}
#define zeigs(m,d,v,n) {int info; zheev_tmp *t = zewrk((m).size1); zheev("V","L",&(m).size1,(m).data,&(m).tda,(d).data,t->work,&t->lwork,t->rwork,&info);zmat_copya(v,m);}


#define dvec_alloc(v,n) { (v).size = n; (v).stride = 1; (v).data = (double *)malloc((v).size*sizeof(double)); }
#define dmat_alloc(m,r,c) { (m).size1 = r; (m).size2 = c; (m).tda = (m).size2; (m).data = (double *)malloc((m).size1*(m).size2*sizeof(double)); }
#define dvec_free(v) free((v).data)
#define dmat_free(m) free((m).data)

#define dvec_get(v,k) (v).data[(v).stride*(k)]
#define dmat_get(m,k1,k2) (m).data[(m).tda*(k1)+(k2)]
#define dvec_set(v,k,r) (v).data[(v).stride*(k)] = (r)
#define dmat_set(m,k1,k2,r) (m).data[(m).tda*(k1)+(k2)] = (r)

#define dvec_zero(v) { for(int _i=0; _i<(v).size; _i++) dvec_set(v,_i,0.); }
#define dmat_zero(m) { for(int _i=0; _i<(m).size1; _i++) for(int _j=0; _j<(m).size2; _j++) dmat_set(m,_i,_j,0.); }

#define dvec_copy(v1,v2) { for(int _i=0; _i<(v1).size; _i++) dvec_set(v1,_i,dvec_get(v2,_i)); }
#define dmat_copy(m1,m2) { for(int _i=0; _i<(m1).size1; _i++) for(int _j=0; _j<(m1).size2; _j++) dmat_set(m1,_i,_j,dmat_get(m2,_i,_j)); }
#define dmat_copya(m1,m2) { for(int _i=0; _i<(m1).size1; _i++) for(int _j=0; _j<(m1).size2; _j++) dmat_set(m1,_i,_j,dmat_get(m2,_j,_i)); }

#define dsubvec(v1,v2,k,n) { (v1).size = n; (v1).stride = (v2).stride; (v1).data = (v2).data + (k); }
#define dcolvec(v,m,k) { (v).size = (m).size1; (v).stride = (m).tda; (v).data = (m).data + (k); }
#define dsubmat(m1,m2,k1,k2,n1,n2) { (m1).size1 = n1; (m1).size2 = n2; (m1).tda = (m2).tda; (m1).data = (m2).data + ((k1)*(m1).tda) + (k2); }

#define dmd(m) {printf("%i %i %i\n", m.size1, m.size2, m.tda);}
#define dmulmnn(m1,m2,m3) {double _zero=0,_one=1; dgemm("N","N",&(m1).size2,&(m1).size1,&(m3).size1,&_one,(m3).data,&(m3).tda,(m2).data,&(m2).tda,&_zero,(m1).data,&(m1).tda);}
#define dmulman(m1,m2,m3) {double _zero=0,_one=1; dgemm("N","T",&(m1).size2,&(m1).size1,&(m3).size1,&_one,(m3).data,&(m3).tda,(m2).data,&(m2).tda,&_zero,(m1).data,&(m1).tda);}

#define dorthocols(y1,y0) dorthocols_func(&y1,&y0)


#define dv_dot_v(v1,v2) ddot(&(v1).size,(v1).data,&(v1).stride,(v2).data,&(v2).stride)
#define dv_peq_r_times_v(v1,r,v2) daxpy(&(v1).size,&r,(v2).data,&(v2).stride,(v1).data,&(v1).stride)
#define dnrm2_v(v) ddot(&(v).size,(v).data,&(v).stride,(v).data,&(v).stride)
#define dv_teq_r(v,r) dscal(&(v).size,&r,(v).data,&(v).stride)

#define deigs(m,d,v,n) {int info; dsyev_tmp *t = ewrk((m).size1); dsyev("V","L",&(m).size1,(m).data,&(m).tda,(d).data,t->work,&t->lwork,&info);dmat_copya(v,m);}
#define dmat_comp(a,b) dmat_comp_func(&a,&b)
#define dvec_print(v) {int _i; for(_i=0; _i<(v).size; _i++) printf("  %-12g\n", dvec_gat(v,_i));}
#define dmat_print(m) {int _i,_j; for(_i=0; _i<(m).size1; _i++) for(_j=0; _j<(m).size2; _j++) printf("  %i %i %-12g\n", _i, _j, dmat_get(m,_i,_j));}


#define printmr(m) printmr_func(m.mat)


#define alloc_vecc(v,n) { (v).vec = gsl_vector_complex_alloc(n); }
#define alloc_matc(m,r,c) { (m).mat = gsl_matrix_complex_alloc(r,c); }
#define free_vecc(v) gsl_vector_complex_free((v).vec)
#define free_matc(m) gsl_matrix_complex_free((m).mat)
#define setc(m,k1,k2,r) {gsl_complex zt; GSL_SET_COMPLEX(&zt,QLA_real(r),QLA_imag(r));gsl_matrix_complex_set((m).mat,k1,k2,zt);}
#define eigsc(m,d,v,n) {gsl_eigen_hermv((m).mat,(d).vec,(v).mat,ewrkc((m).mat->size1)); gsl_eigen_hermv_sort((d).vec,(v).mat,GSL_EIGEN_SORT_VAL_ASC);}


#if 0
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
#endif

static void
dmat_comp_func(dmat *a, dmat *b)
{
  int i, j, n1, n2;
  n1 = a->size1;
  n2 = a->size2;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {
      double av = dmat_get(*a, i, j);
      double bv = dmat_get(*b, i, j);
      printf("%12g %12g %12g\n", av, bv, 0.5*fabs(av-bv)/(av+bv));
    }
  }
}

// orthonormalize v against u
static void
dorthocols_func(dmat *v, dmat *u)
{
  int i, j, ni, nj;
  ni = v->size2;
  nj = u->size2;
  //printf("dorth %i %i\n", ni, nj);
  //dmat_print(*v);
  //printf("u\n");
  //dmat_print(*u);
  for(i=0; i<ni; i++) {
    dvec vv;
    dcolvec(vv, *v, i);
    //dvec_print(vv);
    for(j=0; j<nj; j++) {
      double s;
      dvec uu;
      dcolvec(uu, *u, j);
      //dvec_print(uu);
      s = - dv_dot_v(uu, vv);
      dv_peq_r_times_v(vv, s, uu);
      //dvec_print(vv);
    }
    for(j=0; j<i; j++) {
      double s;
      dvec uu;
      dcolvec(uu, *v, j);
      s = - dv_dot_v(uu, vv);
      dv_peq_r_times_v(vv, s, uu);
      //dvec_print(vv);
    }
    double f;
    f = 1./sqrt(dnrm2_v(vv));
    dv_teq_r(vv, f);
    //printf(" %i %g\n", i, f);
    //dvec_print(vv);
  }
}

typedef struct {
  double *work;
  int lwork;
  int n;
} dsyev_tmp;

static dsyev_tmp *
ewrk(int n)
{
  static dsyev_tmp w;
  static int s=0;
  if(n>s) {
    if(s) free(w.work);
    s = n;
    w.n = n;
    w.lwork = 3*n;
    w.work = (double *) malloc(w.lwork*sizeof(double));
  }
  return &w;
}

typedef struct {
  _Complex double *work;
  double *rwork;
  int lwork;
  int n;
} zheev_tmp;

static zheev_tmp *
zewrk(int n)
{
  static zheev_tmp w;
  static int s=0;
  if(n>s) {
    if(s) { free(w.work); free(w.rwork); }
    s = n;
    w.n = n;
    w.lwork = 3*n;
    w.work = (_Complex double *) malloc(w.lwork*sizeof(_Complex double));
    w.rwork = (double *) malloc(3*n*sizeof(double));
  }
  return &w;
}
