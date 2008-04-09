#if 0
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

#if QOP_Precision == 1

#define v_eq_zero(r, n) {float zero=0.0; int one=1,n2=2*(n); sscal(&n2, &zero, (float *)r, &one);}
#define v_eq_v(r, a, n) {int one=1; ccopy(&n, a, &one, r, &one);}
#define v_peq_v(r, a, n) {float p_one=1.0; int one=1,n2=2*(n); saxpy(&n2, &p_one, (float *)a, &one, (float *)r, &one);}
#define v_meq_v(r, a, n) {float m_one=-1.0; int one=1,n2=2*(n); saxpy(&n2, &m_one, (float *)a, &one, (float *)r, &one);}
#define v_eq_v_minus_v(r, a, b, n) { v_eq_v(r, a, n); v_meq_v(r, b, n); }
#define v_teq_r(r, s, n) {float ss=s; int one=1,n2=2*(n); sscal(&n2, &ss, (float *)r, &one);}
#define v_peq_r_times_v(r, s, a, n) {float ss=s; int one=1,n2=2*(n); saxpy(&n2, &ss, (float *)a, &one, (float *)r, &one);}
#define v_meq_r_times_v(r, s, a, n) {float ms=-s; int one=1,n2=2*(n); saxpy(&n2, &ms, (float *)a, &one, (float *)r, &one);}
#define norm2_v(a, n) ({double lnrm2; int one=1; lnrm2 = scnrm2(&n, a, &one); lnrm2 *= lnrm2; QMP_sum_double(&lnrm2); lnrm2;})
#define re_v_dot_v(a, b, n) ({double ldot; int one=1,n2=2*(n); ldot = dsdot(&n2, (float *)a, &one, (float *)b, &one); QMP_sum_double(&ldot); ldot;})

#else

#define v_eq_zero(r, n) {double zero=0.0; int one=1,n2=2*(n); dscal(&n2, &zero, (double *)r, &one);}
#define v_eq_v(r, a, n) {int one=1; zcopy(&n, a, &one, r, &one);}
#define v_peq_v(r, a, n) {double p_one=1.0; int one=1,n2=2*(n); daxpy(&n2, &p_one, (double *)a, &one, (double *)r, &one);}
#define v_meq_v(r, a, n) {double m_one=-1.0; int one=1,n2=2*(n); daxpy(&n2, &m_one, (double *)a, &one, (double *)r, &one);}
#define v_eq_v_minus_v(r, a, b, n) { v_eq_v(r, a, n); v_meq_v(r, b, n); }
#define v_teq_r(r, s, n) {double ss=s; int one=1,n2=2*(n); dscal(&n2, &ss, (double *)r, &one);}
#define v_peq_r_times_v(r, s, a, n) {double ss=s; int one=1,n2=2*(n); daxpy(&n2, &ss, (double *)a, &one, (double *)r, &one);}
#define v_meq_r_times_v(r, s, a, n) {double ms=-s; int one=1,n2=2*(n); daxpy(&n2, &ms, (double *)a, &one, (double *)r, &one);}
#define norm2_v(a, n) ({double lnrm2; int one=1; lnrm2 = dznrm2(&n, a, &one); lnrm2 *= lnrm2; QMP_sum_double(&lnrm2); lnrm2;})
#define re_v_dot_v(a, b, n) ({double ldot; int one=1,n2=2*(n); ldot = ddot(&n2, (double *)a, &one, (double *)b, &one); QMP_sum_double(&ldot); ldot;})

#endif

typedef void (*linop_blas_t)(QLA_Complex *out, QLA_Complex *in, void *args);

/* regular CG */
QOP_status_t
QOPPC(invert_cg)(linop_blas_t linop,
		 void *args,
		 QOP_invert_arg_t *inv_arg,
		 QOP_resid_arg_t *res_arg,
		 QLA_Complex *out,
		 QLA_Complex *in,
		 int n)
{
  QLA_Real a, b;
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  QLA_Complex *r, *p, *Mp;
  int iteration=0, total_iterations=0, nrestart=-1;
  int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;

  r = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  p = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  Mp = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));

  insq = norm2_v(in, n);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "CG: rsqstop = %g\n", rsqstop);
  rsq = 0;
  oldrsq = rsq;

  while(1) {

    if( (total_iterations==0) ||
	(iteration>=restart_iterations) ||
	(total_iterations>=max_iterations) ||
	(rsq<rsqstop) ) {  /* only way out */

      if( (total_iterations>=max_iterations) ||
	  (nrestart>=max_restarts) ) break;
      nrestart++;

      linop(Mp, out, args);
      iteration = 1;
      total_iterations++;

      v_eq_v_minus_v(r, in, Mp, n);
      rsq = norm2_v(r, n);
      VERB(LOW, "CG: (re)start: iter %i rsq = %g\n", total_iterations, rsq);
      if( (rsq<rsqstop) ||
	  (total_iterations>=max_iterations) ) break;

      v_eq_v(p, r, n);
      //b = 1.0;

    } else {

      //r_eq_re_V_dot_V(&b, Mp, r, subset);
      //b = -a*b/oldrsq;
      b = rsq / oldrsq;
      //v_eq_r_times_v_plus_v(p, b, p, r, n);
      v_teq_r(p, b, n);
      v_peq_v(p, r, n);
      //QLA_Real s = 1.0 / b;
      //v_peq_r_times_v(p, s, r, n);

    }
    oldrsq = rsq;

    linop(Mp, p, args);
    iteration++;
    total_iterations++;

    pkp = re_v_dot_v(p, Mp, n);

    a = rsq / pkp;
    v_peq_r_times_v(out, a, p, n);
    v_meq_r_times_v(r, a, Mp, n);

    //QLA_Real s = a / b;
    //v_peq_r_times_v(out, s, p, n);
    //v_meq_r_times_v(r, s, Mp, n);

    rsq = norm2_v(r, n);
    VERB(MED, "CG: iter %i rsq = %g\n", total_iterations, rsq);
  }

  free(r);
  free(p);
  free(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

  return QOP_SUCCESS;
}

/* milti-shift CG */
QOP_status_t
QOPPC(invert_cgms)(linop_blas_t linop,
		   void *args,
		   QOP_invert_arg_t *inv_arg,
		   QOP_resid_arg_t **res_arg,
		   QLA_Real *shifts,
		   int nshifts,
		   QLA_Complex **out,
		   QLA_Complex *in,
		   int n)
{
  QLA_Real a[nshifts], b[nshifts];
  QLA_Real bo[nshifts], z[nshifts], zo[nshifts], zn[nshifts];
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  int iteration=0, i, imin;
  QLA_Complex *r, *Mp, *pm[nshifts];

  imin = 0;
  for(i=1; i<nshifts; i++) if(shifts[i]<shifts[imin]) imin = i;

  r = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  Mp = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  for(i=0; i<nshifts; i++) {
    pm[i] = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  }

  insq = norm2_v(in, n);
  v_eq_v(r, in, n);
  for(i=0; i<nshifts; i++) {
    v_eq_zero(out[i], n);
    v_eq_v(pm[i], r, n);
    zo[i] = z[i] = 1;
    bo[i] = -1;
    a[i] = 0;
  }

  rsqstop = res_arg[imin]->rsqmin * insq;
  VERB(LOW, "CGMS: rsqstop = %g\n", rsqstop);
  rsq = insq;
  //printf("start %g\n", rsq);

  while(1) {
    oldrsq = rsq;

    linop(Mp, pm[imin], args);
    //if(shifts[imin]!=0.0) V_peq_r_times_V(Mp, shifts+imin, p, subset);
    iteration++;

    pkp = re_v_dot_v(pm[imin], Mp, n);

    b[imin] = rsq / pkp;
    zn[imin] = 1;
    for(i=0; i<nshifts; i++) {
      if(i!=imin) {
	double c1;
	zn[i] = z[i]*zo[i]*bo[imin];
	c1 = b[imin]*a[imin]*(zo[i]-z[i]);
	c1 += zo[i]*bo[imin]*(1+shifts[i]*b[imin]);
	if(c1!=0.0) zn[i] /= c1;
	else zn[i] = 0;
	if(z[i]!=0.0) b[i] = b[imin]*zn[i]/z[i];
	else zn[i] = b[i] = 0;
      }
    }

    for(i=0; i<nshifts; i++) {
      v_peq_r_times_v(out[i], b[i], pm[i], n);
    }

    v_meq_r_times_v(r, b[imin], Mp, n);
    rsq = norm2_v(r, n);
    VERB(MED, "CGMS: iter %i rsq = %g\n", iteration, rsq);

    if( (iteration%inv_arg->restart==0) ||
	(iteration>=inv_arg->max_iter) ||
	(rsq<rsqstop) ) {  /* only way out */
      break;
    }

    a[imin] = rsq / oldrsq;
    for(i=0; i<nshifts; i++) {
      if(i!=imin) {
	double c2 = z[i]*b[imin];
	if(c2!=0.0) a[i] = a[imin]*zn[i]*b[i]/c2;
	else a[i] = 0;
      }
    }

    //V_eq_r_times_V_plus_V(p, a+imin, p, r, subset);
    v_teq_r(pm[imin], a[imin], n);
    v_peq_v(pm[imin], r, n);
    for(i=0; i<nshifts; i++) {
      if(i!=imin) {
	//v_eq_r_times_v(Mp, zn[i], r, n);
	//v_eq_r_times_v_plus_v(pm[i], a+i, pm[i], Mp, subset);
	v_teq_r(pm[i], a[i], n);
	v_peq_r_times_v(pm[i], zn[i], r, n);
      }
    }

    for(i=0; i<nshifts; i++) {
      bo[i] = b[i];
      zo[i] = z[i];
      z[i] = zn[i];
    }
  }

  free(r);
  free(Mp);
  for(i=0; i<nshifts; i++) {
    free(pm[i]);
  }

  for(i=0; i<nshifts; i++) {
    res_arg[i]->final_rsq = rsq/insq;
    res_arg[i]->final_iter = iteration;
    res_arg[i]->final_restart = 0;
  }
  VERB(MED, "CGMS: done: iter %i rsq = %g\n", iteration, rsq);

  return QOP_SUCCESS;
}

/***************************************************************************/

typedef struct {
  Vector *in, *out;
  QDP_Subset subset;
  QOPPCV(linop_t) *linop;
} linop_args_t;

static void
linop_blas(QLA_Complex *out, QLA_Complex *in, void *args)
{
  linop_args_t *a = (linop_args_t *)args;
  insert_packed_V(a->in, (void *)in, a->subset);
  //insert_packed_V(a->out, out, a->subset);
  a->linop(a->out, a->in, a->subset);
  extract_packed_V((void *)out, a->out, a->subset);
  //int n = QDP_subset_len(a->subset)*csize_V;
  //v_eq_v(out, in, n);
}

/* regular CG */
QOP_status_t
QOPPCV(invert_cg)(QOPPCV(linop_t) *linop,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  Vector *out,
		  Vector *in,
		  Vector *p,
		  QDP_Subset subset
		  vIndexDef)
{
  QLA_Complex *bout, *bin;
  linop_args_t args;
  int n = QDP_subset_len(subset)*csize_V;
  QOP_status_t st;

  args.linop = linop;
  //create_V(args.in);
  args.in = p;
  create_V(args.out);
  args.subset = subset;

  bin = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  bout = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));

  extract_packed_V((void *)bin, in, subset);
  extract_packed_V((void *)bout, out, subset);

  //V_eq_V(p, out, subset);
  //linop(args.out, p, subset);
  //V_meq_V(args.out, in, subset);
  //QLA_Real rsq;
  //r_eq_norm2_V(&rsq, args.out, subset);
  //VERB(LOW, "CG: (re)start: iter 0 rsq = %g\n", rsq);

  //linop_blas(bout, bout, &args);
  //v_meq_v(bout, bin, n);
  //VERB(LOW, "CG: (re)start: iter 0 rsq = %g\n", norm2_v(bout, n));

  st = QOPPC(invert_cg)(linop_blas, (void *)&args, inv_arg, res_arg, bout, bin, n);

  insert_packed_V(out, (void *)bout, subset);

  free(bin);
  free(bout);
  //destroy_V(args.in);
  destroy_V(args.out);
  return st;
}

/* milti-shift CG */
QOP_status_t
QOPPCV(invert_cgms)(QOPPCV(linop_t) *linop,
		    QOP_invert_arg_t *inv_arg,
		    QOP_resid_arg_t **res_arg,
		    QLA_Real *shifts,
		    int nshifts,
		    Vector **out,
		    Vector *in,
		    Vector *p,
		    QDP_Subset subset
		    vIndexDef)
{
  QLA_Complex **bout, *bin;
  linop_args_t args;
  int n = QDP_subset_len(subset)*csize_V;
  QOP_status_t st;
  int i;

  args.linop = linop;
  args.in = p;
  create_V(args.out);
  args.subset = subset;

  bin = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  extract_packed_V((void *)bin, in, subset);
  bout = (QLA_Complex **) malloc(nshifts*sizeof(QLA_Complex *));
  for(i=0; i<nshifts; i++) {
    bout[i] = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
    //extract_packed_V((void *)bout, out, subset);
  }

  st = QOPPC(invert_cgms)(linop_blas, (void *)&args, inv_arg, res_arg, shifts, nshifts, bout, bin, n);

  for(i=0; i<nshifts; i++) {
    insert_packed_V(out[i], (void *)bout[i], subset);
  }

  free(bin);
  for(i=0; i<nshifts; i++) {
    free(bout[i]);
  }
  free(bout);
  destroy_V(args.out);
  return st;
}
