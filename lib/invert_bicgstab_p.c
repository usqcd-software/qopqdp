#define printf0 QOP_printf0

QOP_status_t
QOPPCV(invert_bicgstab)(QOPPCV(linop_t) *linop,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg,
			Vector *out,
			Vector *in,
			Vector *p,
			Vector *r,
			QDP_Subset subset
			vIndexDef)
{
#define NC QDP_get_nc(first_qdp_object(in))
  QDP_Lattice *lat = get_lattice_V(r);
  QLA_D_Complex rho0, rho1;
  QLA_D_Complex alpha, beta, omega;
  QLA_D_Complex ctmp1, ctmp2;
  QLA_D_Real t2;
  QLA_D_Real rsq, relnorm2;
  QLA_D_Real insq;
  QLA_D_Real rsqstop;
  Vector *r0, *t, *v;
  int iteration=0, total_iterations=0, nrestart=-1;
  int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;

  // avoid compiler warning about unitialized variables
  QLA_c_eq_r(alpha, 0);
  QLA_c_eq_r(omega, 0);
  QLA_c_eq_r(rho1, 0);

  create_V(r0);
  create_V(t);
  create_V(v);

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "BICG: rsqstop = %g relmin = %g\n", rsqstop, res_arg->relmin);
  rsq = 0;
  relnorm2 = 1.;

  /* Default output values unless reassigned */
  res_arg->final_rsq = 0;
  res_arg->final_rel = 0;
  res_arg->final_iter = 0;
  res_arg->final_restart = 0;

  /* Special case of exactly zero source */
  if(insq == 0.){
    VERB(LOW, "BICG: exiting because of zero source\n");
    destroy_V(r0);
    destroy_V(t);
    destroy_V(v);
    V_eq_zero(out, subset);
    
    return QOP_SUCCESS;
  }

  while(1) {

    if( (total_iterations==0) ||
	(iteration>=restart_iterations) ||
	(total_iterations>=max_iterations) ||
	(rsqstop>0 && rsq<rsqstop) ||
	(res_arg->relmin>0 && relnorm2<res_arg->relmin) ) {
      /* only way out */

      /* stop when we exhaust iterations */
      if( (total_iterations>=max_iterations) ||
	  (nrestart>=max_restarts) ) break;

      /* otherwise, restart */
      nrestart++;

      /* compute true residual */
      V_eq_V(p, out, subset);
      linop(r, p, subset);
      iteration = 1;
      total_iterations++;

      V_eq_V_minus_V(r, in, r, subset);
      V_eq_V(r0,r,subset);
      r_eq_norm2_V(&rsq, r, subset);
      /* compute FNAL norm if requested */
      if(res_arg->relmin > 0)
	relnorm2 = relnorm2_V(r, out, subset);
      VERB(LOW, "BICG: (re)start: iter %i rsq = %g rel = %g\n", 
	   total_iterations, rsq, relnorm2);

      /* stop here if converged */
      if( rsq<rsqstop ||
	  relnorm2<res_arg->relmin ||
	  total_iterations>=max_iterations ) break;

      V_eq_zero(v, subset);
      V_eq_zero(p, subset);

      QLA_c_eq_r(rho0,1.0);
      QLA_c_eq_r(alpha,1.0);
      QLA_c_eq_r(omega,1.0);

      QLA_c_eq_r(rho1,rsq);
      V_eq_V(p, r, subset);

    } else {

      QLA_c_eq_c(rho0,rho1);
      c_eq_V_dot_V(&rho1, r0, r, subset);
      if(QLA_norm2_c(rho1)==0.) break;
      QLA_c_eq_c_times_c(ctmp1,rho1,alpha);
      QLA_c_eq_c_times_c(ctmp2,rho0,omega);
      QLA_D_c_eq_c_div_c(beta,ctmp1,ctmp2);
#define C2(x) QLA_real(x), QLA_imag(x)
      //printf0("beta = %g %g  ctmp1 = %g %g  ctmp2 = %g %g\n", C2(beta), C2(ctmp1), C2(ctmp2));

      QLA_c_eqm_c(omega, omega);
      //V_eq_c_times_V_plus_V(t, &omega, v, p, subset);
      //V_eq_c_times_V_plus_V(p, &beta, t, r, subset);
      {
 	QLA_Complex tc[2];
	Vector *tr[2], *ta[2], *tb[2];
	tr[0] = t; QLA_c_eq_c(tc[0], omega); ta[0] = v; tb[0] = p;
	tr[1] = p; QLA_c_eq_c(tc[1], beta);  ta[1] = t; tb[1] = r;
	V_veq_c_times_V_plus_V(tr, tc, ta, tb, subset, 2);
      }
    }

    linop(v, p, subset);
    iteration++;
    total_iterations++;

    c_eq_V_dot_V(&ctmp1,r0,v,subset);
    if(QLA_norm2_c(ctmp1)==0.) break;
    QLA_D_c_eq_c_div_c(alpha, rho1, ctmp1);
    //printf0("alpha = %g %g  rho1 = %g %g  ctmp1 = %g %g\n", C2(alpha), C2(rho1), C2(ctmp1));
    V_meq_c_times_V(r, &alpha, v, subset);

    linop(t, r, subset);

    c_eq_V_dot_V(&omega, t, r, subset);
    r_eq_norm2_V(&t2, t, subset);
    QLA_c_eq_c_div_r(omega, omega, t2);
    //printf0("omega = %g %g  t2 = %g\n", C2(beta), t2);

    /* update the solution and residuals */
#if 0
    V_peq_c_times_V(out, &omega, r, subset);
    V_peq_c_times_V(out, &alpha, p, subset);
    V_meq_c_times_V(r, &omega, t, subset);
#else
    {
      QLA_Complex tc[3];
      Vector *tr[3], *ta[3];
      tr[0] = out; QLA_c_eq_c(tc[0],omega);  ta[0] = r;
      tr[1] = out; QLA_c_eq_c(tc[1],alpha);  ta[1] = p;
      tr[2] = r;   QLA_c_eqm_c(tc[2],omega); ta[2] = t;
      V_vpeq_c_times_V(tr, tc, ta, subset, 3);
    }
#endif

    r_eq_norm2_V(&rsq, r, subset);

    /* compute FNAL norm if requested */
    if(res_arg->relmin > 0)
      relnorm2 = relnorm2_V(r, out, subset);

    if(QLA_norm2_c(omega)==0.) break;

    VERB(HI, "BICG: iter %i rsq = %g rel = %g\n", total_iterations, rsq, 
	 relnorm2);
  }
  VERB(LOW, "BICG: done: iter %i rsq = %g rel = %g\n", 
       total_iterations, rsq, relnorm2);

  destroy_V(r0);
  destroy_V(t);
  destroy_V(v);

  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_iter = iteration;

  return QOP_SUCCESS;
#undef NC
}
