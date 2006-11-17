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
  QLA_Real a[nshifts], b[nshifts];
  QLA_Real bo[nshifts], z[nshifts], zo[nshifts], zn[nshifts];
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real insq;
  QLA_Real rsqstop;
  int iteration=0, i, imin;
  Vector *r, *Mp, *pm[nshifts];

  imin = 0;
  for(i=1; i<nshifts; i++) if(shifts[i]<shifts[imin]) imin = i;

  r = create_V();
  Mp = create_V();
  for(i=0; i<nshifts; i++) {
    if(i==imin) pm[i] = p;
    else pm[i] = create_V();
  }

  r_eq_norm2_V(&insq, in, subset);
  V_eq_V(r, in, subset);
  for(i=0; i<nshifts; i++) {
    V_eq_zero(out[i], subset);
    V_eq_V(pm[i], r, subset);
    zo[i] = z[i] = 1;
    bo[i] = -1;
    a[i] = 0;
  }

  rsqstop = res_arg[imin]->rsqmin * insq;
  rsq = insq;
  //printf("start %g\n", rsq);

  while(1) {
    oldrsq = rsq;

    linop(Mp, p, subset);
    if(shifts[imin]!=0.0) V_peq_r_times_V(Mp, shifts+imin, p, subset);
    iteration++;

    r_eq_re_V_dot_V(&pkp, p, Mp, subset);

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
      V_peq_r_times_V(out[i], b+i, pm[i], subset);
    }

    V_meq_r_times_V(r, b+imin, Mp, subset);
    r_eq_norm2_V(&rsq, r, subset);
    //printf("      %g\n", rsq);

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

    V_eq_r_times_V_plus_V(p, a+imin, p, r, subset);
    for(i=0; i<nshifts; i++) {
      if(i!=imin) {
	V_eq_r_times_V(Mp, zn+i, r, subset);
	V_eq_r_times_V_plus_V(pm[i], a+i, pm[i], Mp, subset);
      }
    }

    for(i=0; i<nshifts; i++) {
      bo[i] = b[i];
      zo[i] = z[i];
      z[i] = zn[i];
    }
  }

  destroy_V(r);
  destroy_V(Mp);
  for(i=0; i<nshifts; i++) {
    if(i!=imin) destroy_V(pm[i]);
  }

  for(i=0; i<nshifts; i++) {
    res_arg[i]->final_rsq = rsq;
    res_arg[i]->final_iter = iteration;
  }

  return QOP_SUCCESS;
}
