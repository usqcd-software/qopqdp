QOP_status_t
QOPPCV(invert_cg)(QOPPCV(linop_t) *linop,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  Vector *out,
		  Vector *in,
		  Vector *p,
		  QDP_Subset subset)
{
  QLA_Real a, b;
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  int iteration=0;
  Vector *r, *Mp;

  r = create_V();
  Mp = create_V();

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  rsq = 0;
  oldrsq = rsq;

  while(1) {

    if( (iteration%inv_arg->restart==0) ||
	(iteration>=inv_arg->max_iter) ||
	(rsq<rsqstop) ) {  /* only way out */

      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration++;

      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);
      if( (rsq<rsqstop) || (iteration>=inv_arg->max_iter) ) break;
      V_eq_V(p, r, subset);

    } else {

      //r_eq_re_V_dot_V(&b, Mp, r, subset);
      //b = -a*b/oldrsq;
      b = rsq / oldrsq;
      V_eq_r_times_V_plus_V(p, &b, p, r, subset);

    }
    oldrsq = rsq;

    linop(Mp, p, subset);
    iteration++;

    r_eq_re_V_dot_V(&pkp, p, Mp, subset);

    a = rsq / pkp;

    V_peq_r_times_V(out, &a, p, subset);
    V_meq_r_times_V(r, &a, Mp, subset);
    r_eq_norm2_V(&rsq, r, subset);
  }

  destroy_V(r);
  destroy_V(Mp);

  res_arg->final_rsq = rsq;
  res_arg->final_iter = iteration;

  return QOP_SUCCESS;
}

#define sum(res, arr) res=0; for(i=0; i<n; i++) res += arr[i];

QOP_status_t
QOPPCV(invert_cgv)(QOPPCV(linopv_t) *linop,
		   QOP_invert_arg_t *inv_arg,
		   QOP_resid_arg_t *res_arg,
		   Vector *out[],
		   Vector *in[],
		   QDP_Subset subset,
		   int n)
{
  QLA_Real a[n], b[n];
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq, vreal[n];
  QLA_Real rsqstop;
  int i, iteration=0;
  Vector *r[n], *p[n], *Mp[n];

  for(i=0; i<n; i++) {
    r[i] = create_V();
    p[i] = create_V();
    Mp[i] = create_V();
  }

  r_veq_norm2_V(vreal, in, subset, n);
  sum(insq, vreal);
  rsqstop = res_arg->rsqmin * insq;
  rsq = 0;
  oldrsq = rsq;

  while(1) {

    if( (iteration%inv_arg->restart==0) ||
	(iteration>=inv_arg->max_iter) ||
	(rsq<rsqstop) ) {  /* only way out */

      V_veq_V(p, out, subset, n);
      linop(Mp, p, subset);
      iteration++;

      V_veq_V_minus_V(r, in, Mp, subset, n);
      r_veq_norm2_V(vreal, r, subset, n);
      sum(rsq, vreal);
      if( (rsq<rsqstop) || (iteration>=inv_arg->max_iter) ) break;
      V_veq_V(p, r, subset, n);

    } else {

      b[0] = rsq / oldrsq;
      for(i=1; i<n; i++) b[i] = b[0];
      V_veq_r_times_V_plus_V(p, b, p, r, subset, n);

    }
    oldrsq = rsq;

    linop(Mp, p, subset);
    iteration++;

    r_veq_re_V_dot_V(vreal, p, Mp, subset, n);
    sum(pkp, vreal);

    a[0] = rsq / pkp;
    for(i=1; i<n; i++) a[i] = a[0];

    V_vpeq_r_times_V(out, a, p, subset, n);
    V_vmeq_r_times_V(r, a, Mp, subset, n);
    r_veq_norm2_V(vreal, r, subset, n);
    sum(rsq, vreal);
  }

  for(i=0; i<n; i++) {
    destroy_V(r[i]);
    destroy_V(p[i]);
    destroy_V(Mp[i]);
  }

  res_arg->final_rsq = rsq;
  res_arg->final_iter = iteration;

  return QOP_SUCCESS;
}
