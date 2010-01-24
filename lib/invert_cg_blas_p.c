#include "linalg.h"

typedef void (*linop_blas_t)(QLA_Complex *out, QLA_Complex *in, void *args);

/* regular CG */
QOP_status_t
QOPPCV(invert_cg_blas)(linop_blas_t linop,
		       void *args,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       QLA_Complex *out,
		       QLA_Complex *in,
		       int n)
{
  QLA_Real a;
  QLA_Real rsq, oldrsq, pkp, relnorm2;
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
  relnorm2 = 1.;
  oldrsq = rsq;

  /* Default output values unless reassigned */
  res_arg->final_rsq = 0;
  res_arg->final_rel = 0;
  res_arg->final_iter = 0;
  res_arg->final_restart = 0;
  
  /* Special case of exactly zero source */
  if(insq == 0.){
    VERB(LOW, "CG: exiting because of zero source\n");
    free(Mp);
    free(p);
    free(r);
    v_eq_zero(out, n);
    
    return QOP_SUCCESS;
  }

  while(1) {

    if( (total_iterations==0) ||
	(iteration>=restart_iterations) ||
	(total_iterations>=max_iterations) ||
	((rsqstop <= 0 || rsq<rsqstop) &&
	 (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ){
      /* only way out */

      /* stop when we exhaust iterations */
      if( (total_iterations>=max_iterations) ||
	  (nrestart>=max_restarts) ) break;

      /* otherwise, restart */
      nrestart++;

      /* compute true residual */
      linop(Mp, out, args);
      iteration = 1;
      total_iterations++;

      v_eq_v_minus_v(r, in, Mp, n);
      rsq = norm2_v(r, n);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_v(r, out, n);

      VERB(LOW, "CG: (re)start: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);

      /* stop here if converged */
      if( ((rsqstop <= 0 || rsq<rsqstop) &&
	   (res_arg->relmin <= 0 || relnorm2<res_arg->relmin)) ||
	  (total_iterations>=max_iterations) ) break;

      //v_eq_v(p, r, n);
      QLA_Real s = 1.0 / rsq;
      v_eq_r_times_v(p, s, r, n);

    } else {

      //r_eq_re_V_dot_V(&b, Mp, r, subset);
      //b = -a*b/oldrsq;
      //b *= rsq / oldrsq;
      //b = rsq;
      //v_eq_r_times_v_plus_v(p, b, p, r, n);
      //v_teq_r(p, b, n);
      //v_peq_v(p, r, n);
      QLA_Real s = 1.0 / rsq;
      v_peq_r_times_v(p, s, r, n);

    }

    linop(Mp, p, args);
    iteration++;
    total_iterations++;

    pkp = re_v_dot_v(p, Mp, n);

    //a = rsq / pkp;
    //v_peq_r_times_v(out, a, p, n);
    //v_meq_r_times_v(r, a, Mp, n);

    //QLA_Real s = a / b;
    a = 1.0 / pkp;
    v_peq_r_times_v(out, a, p, n);
    v_meq_r_times_v(r, a, Mp, n);

    oldrsq = rsq;
    rsq = norm2_v(r, n);

    /* compute FNAL norm if requested */
    if(res_arg->relmin > 0)
      relnorm2 = relnorm2_v(r, out, n);

    VERB(HI, "CG: iter %i rsq = %g rel = %g\n", 
	 total_iterations, rsq, relnorm2);
  }
  VERB(LOW, "CG: done: iter %i rsq = %g rel = %g\n", 
       total_iterations, rsq, relnorm2);

  free(r);
  free(p);
  free(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

  return QOP_SUCCESS;
}

/* milti-shift CG */
QOP_status_t
QOPPCV(invert_cgms_blas)(linop_blas_t linop,
			 void *args,
			 QOP_invert_arg_t *inv_arg,
			 QOP_resid_arg_t **res_arg,
			 QLA_Real *shifts,
			 int nshifts,
			 QLA_Complex **out,
			 QLA_Complex *in,
			 int n)
{
  double a[nshifts], b[nshifts], s[nshifts];
  double bo[nshifts], z[nshifts], zo[nshifts], zn[nshifts];
  double rsq, oldrsq, pkp, relnorm2;
  double insq;
  double rsqstop;
  QLA_Real t;
  int iteration=0, i, imin, imax;
  QLA_Complex *r, *Mp, *pm[nshifts];

  imin = 0;
  for(i=1; i<nshifts; i++) if(shifts[i]<shifts[imin]) imin = i;

  imax = 0;
  for(i=1; i<nshifts; i++) if(shifts[i]>shifts[imax]) imax = i;

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
    s[i] = 1;
  }

  rsqstop = res_arg[imin]->rsqmin * insq;
  VERB(LOW, "CGMS: rsqstop = %g\n", rsqstop);
  rsq = insq;
  relnorm2 = 0;
  //printf("start %g\n", rsq);

  while(1) {
    oldrsq = rsq;

    linop(Mp, pm[imin], args);
    //if(shifts[imin]!=0.0) V_peq_r_times_V(Mp, shifts+imin, p, subset);
    iteration++;

    pkp = s[imin]*s[imin] * re_v_dot_v(pm[imin], Mp, n);

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
      //v_peq_r_times_v(out[i], b[i], pm[i], n);
      t = b[i] * s[i];
      v_peq_r_times_v(out[i], t, pm[i], n);
    }

    //v_meq_r_times_v(r, b[imin], Mp, n);
    t = b[imin] * s[imin];
    v_meq_r_times_v(r, t, Mp, n);
    rsq = norm2_v(r, n);

    /* compute FNAL norm if requested */
    /* here we look at the largest shift, since the FNAL norm is
       most stringent for that case */
    if(res_arg[imax]->relmin > 0){
      t = b[imax] * s[imax];
      v_meq_r_times_v(r, t, Mp, n);
      relnorm2 = relnorm2_v(r, out[imax], n);
    }

    VERB(HI, "CGMS: iter %i rsq = %g rel = %g\n", iteration, rsq,
	 relnorm2);

    if( (iteration%inv_arg->restart==0) ||
	(iteration>=inv_arg->max_iter) ||
	((rsqstop <= 0 || rsq<rsqstop) &&
	 (res_arg[imax]->relmin <= 0 || relnorm2<res_arg[imax]->relmin)) ) {
      /* only way out */
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
    //v_teq_r(pm[imin], a[imin], n);
    //v_peq_v(pm[imin], r, n);
    s[imin] *= a[imin];
    t = 1.0/s[imin];
    v_peq_r_times_v(pm[imin], t, r, n);
    for(i=0; i<nshifts; i++) {
      if(i!=imin) {
	//v_eq_r_times_v(Mp, zn[i], r, n);
	//v_eq_r_times_v_plus_v(pm[i], a+i, pm[i], Mp, subset);
	//v_teq_r(pm[i], a[i], n);
	//v_peq_r_times_v(pm[i], zn[i], r, n);
	s[i] *= a[i];
	t = zn[i]/s[i];
	v_peq_r_times_v(pm[i], t, r, n);
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
    res_arg[i]->final_rel = relnorm2;
    res_arg[i]->final_iter = iteration;
    res_arg[i]->final_restart = 0;
  }
  VERB(LOW, "CGMS: done: iter %i rsq = %g\n", iteration, rsq);

  return QOP_SUCCESS;
}

/***************************************************************************/

typedef struct {
  Vector *in, *out;
  QDP_Subset subset;
  QOPPCV(linop_t) *linop;
  int _n;
} linop_args_t;

static void
linop_blas(QLA_Complex *out, QLA_Complex *in, void *args)
{
  linop_args_t *a = (linop_args_t *)args;
  int _n = a->_n;
  insert_packed_V(a->in, in, a->subset);
  //insert_packed_V(a->out, out, a->subset);
  a->linop(a->out, a->in, a->subset);
  extract_packed_V(out, a->out, a->subset);
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
  args._n = _N;

  bin = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));
  bout = (QLA_Complex *) malloc(n*sizeof(QLA_Complex));

  extract_packed_V(bin, in, subset);
  extract_packed_V(bout, out, subset);

  //V_eq_V(p, out, subset);
  //linop(args.out, p, subset);
  //V_meq_V(args.out, in, subset);
  //QLA_Real rsq;
  //r_eq_norm2_V(&rsq, args.out, subset);
  //VERB(LOW, "CG: (re)start: iter 0 rsq = %g\n", rsq);

  //linop_blas(bout, bout, &args);
  //v_meq_v(bout, bin, n);
  //VERB(LOW, "CG: (re)start: iter 0 rsq = %g\n", norm2_v(bout, n));

  st = QOPPCV(invert_cg_blas)(linop_blas, (void *)&args, inv_arg, res_arg, bout, bin, n);

  insert_packed_V(out, bout, subset);

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

  st = QOPPCV(invert_cgms_blas)(linop_blas, (void *)&args, inv_arg, res_arg, shifts, nshifts, bout, bin, n);

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
