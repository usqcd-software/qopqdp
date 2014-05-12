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
  int iteration=0;
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

  int imin=0, irsq=0, irel=0;
  for(int i=1; i<nshifts; i++) {
    if(shifts[i]<shifts[imin]) imin = i;
    if(res_arg[i]->rsqmin>=res_arg[irsq]->rsqmin) irsq = i;
    if(res_arg[i]->relmin>res_arg[irel]->relmin) irel = i;
  }
  double rsqstop = res_arg[irsq]->rsqmin * insq;
  double relstop = res_arg[irel]->relmin;

  create_V(r);
  create_V(Mp);
  V_eq_V(r, in, subset);
  for(int i=0; i<nshifts; i++) {
    if(i==imin) pm[i] = p;
    else create_V(pm[i]);
    V_eq_zero(out[i], subset);
    V_eq_V(pm[i], r, subset);
    zo[i] = z[i] = 1;
    bo[i] = -1;
    a[i] = 0;
  }

  VERB(MED, "CGMS: rsqmin: %g  relmin: %g  in2: %g  rsqstop: %g\n",
       res_arg[irsq]->rsqmin, relstop, insq, rsqstop);
  rsq = insq;
  relnorm2 = 1;
  int rsqCheckInterval = 100;
  int nextRsqCheck = rsqCheckInterval;
  double trueRsqFac = 1;
  int relCheckInterval = 100;
  int nextRelCheck = relCheckInterval;
  double trueRelFac = 1;

  double yAy = 0;
  Vector *y, *Ay, *ry;
  create_V(y);
  create_V(Ay);
  create_V(ry);
  V_eq_zero(y, subset);
  V_eq_zero(Ay, subset);
  V_eq_V(ry, in, subset);

  while(1) {
    oldrsq = rsq;

    pkp = linop(Mp, p, subset);
    //if(shifts[imin]!=0.0) V_peq_r_times_V(Mp, shifts+imin, p, subset);
    iteration++;

    //r_eq_re_V_dot_V(&pkp, p, Mp, subset);
    if(pkp<=0) break;  // loss of precision in calculating pkp

    //#define CGMS_FLEX1
#ifndef CGMS_FLEX1
    b[imin] = rsq / pkp;
    V_peq_r_times_V(out[imin], b+imin, p, subset);
    V_meq_r_times_V(r, b+imin, Mp, subset);
#else
    {
      QLA_D_Complex pr, bb;
      c_eq_V_dot_V(&pr, p, r, subset);
      QLA_c_eq_c_div_r(bb, pr, pkp);
      V_peq_c_times_V(out[imin], &bb, p, subset);
      V_meq_c_times_V(r, &bb, Mp, subset);
      b[imin] = QLA_real(bb);
    }
#endif
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
      if(i!=imin) {
	V_peq_r_times_V(out[i], b+i, pm[i], subset);
      }
    }

    //#define CGMS_FLEX2
#ifndef CGMS_FLEX2
    r_eq_norm2_V(&rsq, r, subset);
#else
    QLA_D_Complex aa;
    {
      QLA_Complex ca[2];
      Vector *va[2] = {r, Mp}, *ra[2] = {r, r};
      c_veq_V_dot_V(ca, va, ra, subset, 2);
      rsq = QLA_real(ca[0]);
      QLA_c_eq_c_div_r(aa, ca[1], -pkp);
    }
#endif
    double rsqc = zn[irsq]*zn[irsq]*rsq;

    if(relstop > 0) { /* compute FNAL norm if requested */
      relnorm2 = zn[irel]*zn[irel]*relnorm2_V(r, out[irel], subset);
    }

    VERB(HI, "CGMS: iter %i rsq = %g rel = %g\n", iteration, rsqc, relnorm2);

    if( iteration>=inv_arg->restart || iteration>=inv_arg->max_iter ) break;
    if( rsqc*trueRsqFac<rsqstop ) nextRsqCheck = iteration;
    if( relnorm2*trueRelFac<relstop ) nextRelCheck = iteration;

    if(nextRsqCheck<=iteration || (relstop>0 && nextRelCheck<=iteration)) {
      // update improved solution
      if(yAy>0) {  // A-orthogonalize x to y: x = x - (y.A.x/y.A.y) y
	QLA_D_Complex yAx, c;
	c_eq_V_dot_V(&yAx, Ay, out[imin], subset);
	QLA_c_eq_c_div_r(c, yAx, -yAy);
	//QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
	V_peq_c_times_V(out[imin], &c, y, subset);
	//V_eq_c_times_V_plus_V(x, &c, y, out[imin], subset);
      } else {
	//V_eq_V(y, out[irsq], subset);
      }
      // y += c * x, r -= c A x  (c = x.r/x.A.x)
      double xAx = linop(Mp, out[imin], subset);
      iteration++;
      if(xAx<=0) break;
      QLA_D_Complex xr, c;
      c_eq_V_dot_V(&xr, out[imin], ry, subset);
      QLA_c_eq_c_div_r(c, xr, xAx);
      //QLA_c_eq_r(c, 1);
      //QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
      V_peq_c_times_V(y, &c, out[imin], subset);
      //V_meq_c_times_V(ry, &c, Mp, subset);
      yAy = linop(Ay, y, subset);
      iteration++;
      V_eq_V_minus_V(ry, in, Ay, subset);
      V_eq_zero(out[imin], subset);
    }

    if(nextRsqCheck<=iteration) {
      // calculate true residual
      double tr2, en;
      if(irsq==imin) {
	r_eq_norm2_V(&tr2, ry, subset);
	V_eq_V_plus_V(Mp, in, ry, subset);
	r_eq_re_V_dot_V(&en, y, Mp, subset);
      } else {
	double xAx = linop(Mp, out[irsq], subset);
	iteration++;
	if(xAx<=0) break;
	V_meq_V(Mp, in, subset);
	r_eq_norm2_V(&tr2, Mp, subset);
	// calculate energy norm of residual
	V_meq_V(Mp, in, subset);
	r_eq_re_V_dot_V(&en, out[irsq], Mp, subset);
	en = -en;
      }
      VERB(MED, "CGMS: %i: rsq: %g  tr2: %g  enrm: %.16g\n",
	   iteration, rsq, tr2, en);
      if(tr2<=rsqstop) { rsq=tr2; break; }
      if(tr2>2*rsq) break;
      trueRsqFac = tr2/rsq;
      nextRsqCheck += rsqCheckInterval;
    }
    if(relstop>0 && nextRelCheck<=iteration) {
      // calculate true residual
      double xAx = linop(Mp, out[irel], subset);
      iteration++;
      if(xAx<=0) break;
      V_meq_V(Mp, in, subset);
      double rel2 = relnorm2_V(Mp, out[irel], subset);
      VERB(MED, "%i: rel2: %g  trel2: %g\n", iteration, relnorm2, rel2);
      if(rel2<=relstop) { relnorm2 = rel2; break; }
      //if(rel2>1.4*relnorm2) break;
      trueRelFac = rel2/relnorm2;
      nextRelCheck += relCheckInterval;
      relnorm2 = rel2;
    }

#ifndef CGMS_FLEX2
    a[imin] = rsq / oldrsq;
    V_eq_r_times_V_plus_V(p, a+imin, p, r, subset);
#else
    a[imin] = QLA_real(aa);
    V_eq_c_times_V(Mp, &aa, p, subset);
    V_eq_V_plus_V(p, r, Mp, subset);
#endif
    for(int i=0; i<nshifts; i++) {
      if(i!=imin) {
	double c2 = z[i]*b[imin];
	if(c2!=0.0) a[i] = a[imin]*zn[i]*b[i]/c2;
	else a[i] = 0;
      }
    }

    //V_eq_r_times_V_plus_V(p, a+imin, p, r, subset);
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
  V_eq_V(out[imin], y, subset);

  destroy_V(y);
  destroy_V(Ay);
  destroy_V(ry);

  destroy_V(r);
  destroy_V(Mp);
  for(int i=0; i<nshifts; i++) {
    if(i!=imin) destroy_V(pm[i]);
  }

  for(int i=0; i<nshifts; i++) {
    double s = zn[i]*zn[i];
    res_arg[i]->final_rsq = s*trueRsqFac*rsq/insq;
    res_arg[i]->final_rel = s*trueRelFac*relnorm2;
    res_arg[i]->final_iter = iteration;
    res_arg[i]->final_restart = 0;
  }
  VERB(MED, "CGMS: done: iter: %i  rsq: %g  rsq/insq: %g\n",
       iteration, rsq, rsq/insq);

  return QOP_SUCCESS;
#undef NC
}

#if 0
/* multi-shift CG - experimental */
QOP_status_t
QOPPCV(invert_cgms2)(QOPPCV(linopn_t) *linop,
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
  double enrm=0, oldenrm, yAy = 0;
  Vector *y, *Ay, *ry, *x;
  create_V(y);
  create_V(Ay);
  create_V(ry);
  create_V(x);
  V_eq_zero(y, subset);
  V_eq_zero(Ay, subset);
  V_eq_V(ry, in, subset);

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

    // p <- r + aa p | aa = -p.A.r/p.A.p
    //QLA_D_Complex pAr, aa;
    //c_eq_V_dot_V(&pAr, Mp, r, subset);
    //QLA_c_eq_c_div_r(aa, pAr, -pkp);

    if(iteration%100==0) {
      if(yAy>0) {  // A-orthogonalize x to y: x = x - (y.A.x/y.A.y) y
	QLA_D_Complex yAx, c;
	c_eq_V_dot_V(&yAx, Ay, out[imin], subset);
	QLA_c_eq_c_div_r(c, yAx, -yAy);
	//QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
	//V_peq_c_times_V(out[imin], &c, y, subset);
	V_eq_c_times_V_plus_V(x, &c, y, out[imin], subset);
      } else {
	V_eq_V(x, out[imin], subset);
      }
      double enrmout;
      int done = 0;
      {  // y += c * x, r -= c A x  (c = x.r/x.A.x)
	double xAx = linop(Mp, x, subset);
	if(xAx<=0) break;
	QLA_D_Complex xr, c;
	c_eq_V_dot_V(&xr, x, ry, subset);
	QLA_c_eq_c_div_r(c, xr, xAx);
	//QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
	V_peq_c_times_V(y, &c, x, subset);
	QLA_c_eqm_c(c, c);
	V_peq_c_times_V(ry, &c, Mp, subset);
	//V_eq_zero(out[imin], subset);
	r_eq_norm2_V(&xAx, ry, subset);
	// enorm(out)
	enrmout = linop(Mp, out[imin], subset);
	r_eq_re_V_dot_V(&yAy, out[imin], in, subset);
	enrmout = 2*yAy - enrmout;
	V_meq_V(Mp, in, subset);
	double tr2;
	r_eq_norm2_V(&tr2, Mp, subset);
	VERB(MED, "%i: ry2: %g  r2: %g  tr2: %g  r2/tr2: %g\n", iteration, xAx, rsq, tr2, rsq/tr2);
	if(rsq<0.7*tr2) done = 1;
      }
      {
	yAy = linop(Ay, y, subset);
	double yb;
	r_eq_re_V_dot_V(&yb, y, in, subset);
	oldenrm = enrm;
	//enrm = 2*yb - yAy;
	double yr; r_eq_re_V_dot_V(&yr, y, ry, subset);
	enrm = yb + yr;
	VERB(MED, "%i: enrm: %.16g  diff: %g  enrmo: %.16g\n", iteration, enrm, enrm-oldenrm, enrmout);
	if(enrm<=oldenrm) {
	  VERB(LOW, "enrm: %.6f < oldenrm: %.6f\n", enrm, oldenrm);
	  break;
	}
      }
      if(done) break;
    }

    int doReliable = (iteration%10==-1);
    if(doReliable) { // reliable update
      linop(r, out[imin], subset);
      //iteration++;
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
      //pkp = linop(Mp, p, subset);
      //c_eq_V_dot_V(&pAr, Mp, r, subset);
      //QLA_c_eq_c_div_r(aa, pAr, -pkp);
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
    //V_eq_c_times_V(Mp, &aa, p, subset);
    //V_eq_V_plus_V(p, r, Mp, subset);
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
#if 0
  if(yAy>0) {  // A-orthogonalize x to y: x = x - (y.A.x/y.A.y) y
    QLA_D_Complex yAx, c;
    c_eq_V_dot_V(&yAx, Ay, out[imin], subset);
    QLA_c_eq_c_div_r(c, yAx, -yAy);
    //QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
    V_peq_c_times_V(out[imin], &c, y, subset);
  }
  {  // y += c * x, r -= c A x  (c = x.r/x.A.x)
    double xAx = linop(Mp, out[imin], subset);
    if(xAx>0) {
      QLA_D_Complex xr, c;
      c_eq_V_dot_V(&xr, out[imin], ry, subset);
      QLA_c_eq_c_div_r(c, xr, xAx);
      //QOP_printf0("c: %g %g\n", QLA_real(c), QLA_imag(c));
      V_peq_c_times_V(y, &c, out[imin], subset);
    } else {
      VERB(LOW, "xAx<=0: %g\n", xAx);
    }
  }
  V_eq_V(out[imin], y, subset);
#endif
  destroy_V(y);
  destroy_V(Ay);
  destroy_V(ry);
  destroy_V(x);

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
#endif
