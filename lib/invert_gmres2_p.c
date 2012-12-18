#define printf0 QOP_printf0

#define NMAX 32
static int nv=0, nvalloc=0;
static int level[NMAX];
static Vector *vec[NMAX];
static Vector *Avec[NMAX];
static QLA_D_Real Avn[NMAX];
static QLA_D_Complex alpha[NMAX];
static QLA_D_Complex *beta[NMAX];
#if QOP_Colors == 'N'
static int gnc;
#define NC gnc
#endif

static void
init(void)
{
  nv = 0;
  nvalloc = 0;
}

static void
fini(void)
{
  for(int i=0; i<nvalloc; i++) {
    destroy_V(vec[i]);
    destroy_V(Avec[i]);
    free(beta[i]);
  }
  nv = 0;
  nvalloc = 0;
}

static void
combine(int n, QDP_Subset subset)
{
  QLA_D_Complex c;
  QLA_c_eq_c_div_c(c, alpha[n-1], alpha[n]);
  V_peq_c_times_V(Avec[n], &c, Avec[n-1], subset);
  { Vector *tv = Avec[n]; Avec[n] = Avec[n-1]; Avec[n-1] = tv; }
  Avn[n-1] = Avn[n] + QLA_norm2_c(c) * Avn[n-1];
  QLA_c_peq_c(c, beta[n][n-1]);
  V_peq_c_times_V(vec[n], &c, vec[n-1], subset);
  { Vector *tv = vec[n]; vec[n] = vec[n-1]; vec[n-1] = tv; }
  QLA_c_eq_c(alpha[n-1], alpha[n]);
  for(int i=0; i<n-1; i++) {
    QLA_c_peq_c_times_c(beta[n][i], c, beta[n-1][i]);
    QLA_c_eq_c(beta[n-1][i], beta[n][i]);
  }
}

static void
addvec(QDP_Subset subset)
{
  if(nv>1) {
    int n = nv-1;
    while(n>0 && level[n-1]<=level[n]+2) {
      combine(n, subset);
      n--;
      level[n]++;
    }
    nv = n + 1;
  }

  level[nv] = 1;
  nv++;
  if(nv>nvalloc) {
    for(int i=nvalloc; i<nv; i++) {
      create_V(vec[i]);
      create_V(Avec[i]);
      beta[i] = malloc(i*sizeof(*beta[i]));
    }
    nvalloc = nv;
  }
}

static void
orth(QDP_Subset subset)
{
  int n = nv-1;
#if 0
  for(int i=0; i<n; i++) {
    QLA_D_Complex z;
    c_eq_V_dot_V(&z, Avec[i], Avec[n], subset);
    QLA_c_eq_c_div_r(beta[n][i], z, -Avn[i]);
    V_peq_c_times_V(Avec[n], &beta[n][i], Avec[i], subset);
    //V_peq_c_times_V(vec[n], &beta[n][i], vec[i], subset);
  }
#else
  QLA_Complex z[n];
  Vector /* *vp[n],*/ *vAp[n];
  for(int i=0; i<n; i++) { /*vp[i] = vec[n];*/ vAp[i] = Avec[n]; }
  c_veq_V_dot_V(z, Avec, vAp, subset, n);
  for(int i=0; i<n; i++) { QLA_c_eq_c_div_r(beta[n][i], z[i], -Avn[i]); QLA_c_eq_c(z[i], beta[n][i]); }
  V_vpeq_c_times_V(vAp, z, Avec, subset, n);
  //V_vpeq_c_times_V(vp, z, vec, subset, n);
#endif
}

static void
getx(Vector *x, QDP_Subset subset)
{
  QLA_D_Complex b[nv];
  for(int i=0; i<nv-1; i++) QLA_c_eq_c(b[i], beta[nv-1][i]);
  V_eq_c_times_V(x, &alpha[nv-1], vec[nv-1], subset);
  for(int i=nv-2; i>=0; i--) {
    QLA_D_Complex c;
    QLA_c_eq_c_plus_c(c, alpha[i], b[i]);
    V_peq_c_times_V(x, &c, vec[i], subset);
    for(int j=0; j<i; j++) {
      QLA_c_peq_c_times_c(b[j], c, beta[i][j]);
    }
  }
}

QOP_status_t
QOPPCV(invert_gmres2)(QOPPCV(linop_t) *linop,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      Vector *out,
		      Vector *in,
		      Vector *r,
		      QDP_Subset subset
		      vIndexDef)
{
  QLA_D_Real rsq, insq, rsqstop, relnorm2=0;
  int iteration=0, total_iterations=0;//, nrestart=-1;
  //int restart_iterations=inv_arg->restart;
  int max_iterations=inv_arg->max_iter;
  int max_restarts=inv_arg->max_restarts;
  if(max_restarts<0) max_restarts = 5;
#if QOP_Colors == 'N'
  gnc = QDP_get_nc(first_qdp_object(in));
#endif

  /* Default output values unless reassigned */
  res_arg->final_rsq = 0;
  res_arg->final_rel = 0;
  res_arg->final_iter = 0;
  res_arg->final_restart = 0;

  r_eq_norm2_V(&insq, in, subset);
  rsqstop = res_arg->rsqmin * insq;
  VERB(LOW, "GMRES2: rsqstop = %g\n", rsqstop);

  /* Special case of exactly zero source */
  if(insq == 0.) {
    VERB(LOW, "GMRES2: exiting because of zero source\n");
    V_eq_zero(out, subset);
    return QOP_SUCCESS;
  }

  init();
  V_eq_V(r, in, subset);
  r_eq_norm2_V(&rsq, r, subset);

  while(rsq>rsqstop && total_iterations<max_iterations) {
    iteration++;
    total_iterations++;
    addvec(subset);
    V_eq_V(vec[nv-1], r, subset);
    linop(Avec[nv-1], r, subset);
    orth(subset);

    //r_eq_norm2_V(&rsq, r, subset);
    //r_eq_norm2_V(&Avn[nv-1], Avec[nv-1], subset);
    {
      QLA_Real rr[2];
      Vector *vv[2];
      vv[0] = r; vv[1] = Avec[nv-1];
      r_veq_norm2_V(rr, vv, subset, 2);
      rsq = rr[0]; Avn[nv-1] = rr[1];
    }
    {
      QLA_D_Complex ctmp;
      c_eq_V_dot_V(&ctmp, Avec[nv-1], r, subset);
      QLA_c_eq_c_div_r(alpha[nv-1], ctmp, Avn[nv-1]);
    }
    V_meq_c_times_V(r, &alpha[nv-1], Avec[nv-1], subset);
    //r_eq_norm2_V(&rsq, r, subset);
    rsq -= QLA_norm2_c(alpha[nv-1])*Avn[nv-1]; // estimate

    VERB(HI, "GMRES2: iter %i rsq = %g rel = %g\n", total_iterations, rsq, 
	 relnorm2);
#if 0
    {
      getx(r, subset);
      linop(out, r, subset);
      V_eq_V_minus_V(r, in, out, subset);
      r_eq_norm2_V(&rsq, r, subset);
      VERB(HI, "GMRES2: iter %i rsq = %g rel = %g\n", total_iterations, rsq, 
	   relnorm2);
    }
#endif
  }
  VERB(LOW, "GMRES2: done: iter %i rsq = %g rel = %g\n", 
       total_iterations, rsq, relnorm2);

  getx(out, subset);
  fini();

  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_iter = iteration;

  return QOP_SUCCESS;
}
