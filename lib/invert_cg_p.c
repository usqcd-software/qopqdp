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
#define NC QDP_get_nc(first_qdp_object(in))
  QLA_D_Real a, b;
  QLA_D_Real rsq, oldrsq, pkp, relnorm2;
  QLA_D_Real insq;
  QLA_D_Real rsqstop;
  Vector *r, *Mp;
  int iteration=0, total_iterations=0, nrestart=-1;
  int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;

  /* Default output values unless reassigned */
  res_arg->final_rsq = 0;
  res_arg->final_rel = 0;
  res_arg->final_iter = 0;
  res_arg->final_restart = 0;

  r_eq_norm2_V(&insq, in, subset);
  /* Special case of exactly zero source */
  if(insq == 0.) {
    VERB(LOW, "CG: exiting because of zero source\n");
    V_eq_zero(out, subset);
    return QOP_SUCCESS;
  }

  create_V(r);
  create_V(Mp);

  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "CG: rsqmin = %g relmin = %g\n", res_arg->rsqmin, res_arg->relmin);
  VERB(LOW, "CG: rsqstop = %g\n", rsqstop);
  rsq = 0;
  relnorm2 = 1.;
  oldrsq = rsq;
  while(1) {

    if( (total_iterations==0) ||
	(iteration>=restart_iterations) ||
	(total_iterations>=max_iterations) ||
	(rsq<rsqstop) ||
	(relnorm2<res_arg->relmin) ) {
      /* only way out */

      /* stop when we exhaust iterations */
      if( (total_iterations>=max_iterations) ||
	  (nrestart>=max_restarts) ) break;

      /* otherwise, restart */
      nrestart++;

      /* compute true residual */
      V_eq_V(p, out, subset);
      linop(Mp, p, subset);
      iteration = 1;
      total_iterations++;

      V_eq_V_minus_V(r, in, Mp, subset);
      r_eq_norm2_V(&rsq, r, subset);

      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);

      VERB(LOW, "CG: (re)start: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);

      /* stop here if converged */
      if( (rsq<rsqstop) ||
	  (relnorm2<res_arg->relmin) ||
	  (total_iterations>=max_iterations) ) break;

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
    total_iterations++;

    r_eq_re_V_dot_V(&pkp, p, Mp, subset);
    if(pkp<=0) break;  // loss of precision in calculating pkp
    a = rsq / pkp;

    V_peq_r_times_V(out, &a, p, subset);
    V_meq_r_times_V(r, &a, Mp, subset);
    r_eq_norm2_V(&rsq, r, subset);

    /* compute FNAL norm if requested */
    if(res_arg->relmin > 0)
      relnorm2 = relnorm2_V(r, out, subset);

    VERB(HI, "CG: iter %i rsq = %g rel = %g\n", 
	 total_iterations, rsq, relnorm2);
  }
  VERB(LOW, "CG: done: iter %i rsq = %g rel = %g\n", 
       total_iterations, rsq, relnorm2);

  destroy_V(r);
  destroy_V(Mp);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_iter = total_iterations;
  res_arg->final_restart = nrestart;

  return QOP_SUCCESS;
#undef NC
}

/* multi-shift CG */
QOP_status_t
QOPPCV(invert_cgms)(QOPPCV(linopn_t) *linop,
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
#define NC QDP_get_nc(first_qdp_object(in))
  QLA_D_Real a[nshifts], b[nshifts];
  QLA_D_Real bo[nshifts], z[nshifts], zo[nshifts], zn[nshifts];
  QLA_D_Real rsq, oldrsq, pkp, relnorm2;
  QLA_D_Real insq;
  QLA_D_Real rsqstop;
  int iteration=0, imin, imax;
  Vector *r, *Mp, *pm[nshifts];

  r_eq_norm2_V(&insq, in, subset);
  /* Special case - source is exactly zero */
  if(insq == 0.) {
    for(int i=0; i<nshifts; i++) {
      V_eq_zero(out[i], subset);
      res_arg[i]->final_rsq = 0;
      res_arg[i]->final_rel = 0;
      res_arg[i]->final_iter = 0;
      res_arg[i]->final_restart = 0;
    }
    return QOP_SUCCESS;
  }

  imin = 0;
  for(int i=1; i<nshifts; i++) if(shifts[i]<shifts[imin]) imin = i;

  imax = 0;
  for(int i=1; i<nshifts; i++) if(shifts[i]>shifts[imax]) imax = i;

  create_V(r);
  create_V(Mp);
  for(int i=0; i<nshifts; i++) {
    if(i==imin) pm[i] = p;
    else create_V(pm[i]);
  }

  V_eq_V(r, in, subset);
  for(int i=0; i<nshifts; i++) {
    V_eq_zero(out[i], subset);
    V_eq_V(pm[i], r, subset);
    zo[i] = z[i] = 1;
    bo[i] = -1;
    a[i] = 0;
  }

  rsqstop = res_arg[imin]->rsqmin * insq;
  VERB(MED, "CGMS: rsqmin = %g relmin = %g\n", res_arg[imin]->rsqmin,
       res_arg[imin]->relmin);
  VERB(MED, "CGMS: rsqstop = %g\n", rsqstop);
  rsq = insq;
  relnorm2 = 1;
  //printf("start %g\n", rsq);
  double enrm=0, oldenrm;

  while(1) {
    oldrsq = rsq;

    pkp = linop(Mp, p, subset);
    //if(shifts[imin]!=0.0) V_peq_r_times_V(Mp, shifts+imin, p, subset);
    iteration++;

    //r_eq_re_V_dot_V(&pkp, p, Mp, subset);
    if(pkp<=0) break;  // loss of precision in calculating pkp

    b[imin] = rsq / pkp;
    zn[imin] = 1;
    for(int i=0; i<nshifts; i++) {
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

    for(int i=0; i<nshifts; i++) {
      V_peq_r_times_V(out[i], b+i, pm[i], subset);
    }

    V_meq_r_times_V(r, b+imin, Mp, subset);
    r_eq_norm2_V(&rsq, r, subset);

    /* compute FNAL norm if requested */
    /* here we look at the largest shift, since the FNAL norm is
       most stringent for that case */
    if(res_arg[imax]->relmin > 0) {
      V_meq_r_times_V(r, b+imax, Mp, subset);
      relnorm2 = relnorm2_V(r, out[imax], subset);
    }

    VERB(HI, "CGMS: iter %i rsq = %g rel = %g\n", iteration, rsq, relnorm2);

    if( (iteration%inv_arg->restart==0) ||
	(iteration>=inv_arg->max_iter) ||
	(rsq<rsqstop) ||
	(relnorm2<res_arg[imax]->relmin) ) {
      /* only way out */
      break;
    }

    if(iteration%100==0) {
      double xAx = linop(Mp, out[imin], subset);
      double xb;
      r_eq_re_V_dot_V(&xb, out[imin], in, subset);
      oldenrm = enrm;
      enrm = 2*xb - xAx;
      VERB(MED, "enrm: %g\n", enrm);
      if(enrm<oldenrm) {
	VERB(LOW, "enrm: %g < oldenrm: %g\n", enrm, oldenrm);
	break;
      }
    }

    int doReliable = (iteration%20==-1);
    if(doReliable) { // reliable update
      linop(r, out[imin], subset);
      iteration++;
      V_eq_V_minus_V(r, in, r, subset);
      //QOP_printf0("rsq: %g\n", rsq);
      r_eq_norm2_V(&rsq, r, subset);
      //QOP_printf0("rsq: %g\n", rsq);
      //r_eq_re_V_dot_V(&rsq, p, r, subset);
      //QOP_printf0("rsq2: %g\n", rsq);
      // make p orthogonal to r, p <- p + c r, c = -r.p/r.r
      QLA_D_Complex rp, c;
      c_eq_V_dot_V(&rp, r, p, subset);
      QLA_c_eq_c_div_r(c, rp, -rsq);
      V_peq_c_times_V(p, &c, r, subset);
      //c_eq_V_dot_V(&rp, r, p, subset);
      //QOP_printf0("rp: %g %g\n", QLA_real(rp), QLA_imag(rp));
      //QLA_c_eq_c_div_r(c, rp, -rsq);
      //V_peq_c_times_V(p, &c, r, subset);
      //c_eq_V_dot_V(&rp, r, p, subset);
      //QOP_printf0("rp: %g %g\n", QLA_real(rp), QLA_imag(rp));
      //double pAr;
      //r_eq_re_V_dot_V(&pAr, Mp, r, subset);
      //a[imin] = -pAr / pkp;
      a[imin] = rsq / oldrsq;
    } else {
      a[imin] = rsq / oldrsq;
    }
    for(int i=0; i<nshifts; i++) {
      if(i!=imin) {
	double c2 = z[i]*b[imin];
	if(c2!=0.0) a[i] = a[imin]*zn[i]*b[i]/c2;
	else a[i] = 0;
      }
    }

    V_eq_r_times_V_plus_V(p, a+imin, p, r, subset);
    for(int i=0; i<nshifts; i++) {
      if(i!=imin) {
	V_eq_r_times_V(Mp, zn+i, r, subset);
	V_eq_r_times_V_plus_V(pm[i], a+i, pm[i], Mp, subset);
      }
    }

    for(int i=0; i<nshifts; i++) {
      bo[i] = b[i];
      zo[i] = z[i];
      z[i] = zn[i];
    }
  }

  destroy_V(r);
  destroy_V(Mp);
  for(int i=0; i<nshifts; i++) {
    if(i!=imin) destroy_V(pm[i]);
  }

  for(int i=0; i<nshifts; i++) {
    res_arg[i]->final_rsq = rsq/insq;
    res_arg[i]->final_rel = relnorm2;
    res_arg[i]->final_iter = iteration;
    res_arg[i]->final_restart = 0;
  }
  VERB(MED, "CGMS: done: iter %i rsq = %g\n", iteration, rsq);

  return QOP_SUCCESS;
#undef NC
}
