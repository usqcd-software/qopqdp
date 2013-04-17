#include <qop_internal.h>

#ifdef USE_MG

#include <qdp_fn.h>
#include <qdp_dn.h>
#include <qdp_dfn.h>
#include "solvers.h"
#include <float.h>

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QCDP(x) QOP_F_##x
#define QCDPC(x) QOP_F3_##x
#define QDPN(x) QDP_FN_##x
#define QLAN(x) QLA_FN_##x
#else
#define QOPP(x) QOP_D_##x
#define QCDP(x) QOP_D_##x
#define QCDPC(x) QOP_D3_##x
#define QDPN(x) QDP_DN_##x
#define QLAN(x) QLA_DN_##x
#endif

#define printf0 QOP_printf0

void
QCDP(mgVcycle)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  double t0 = QDP_time();
  QCDP(MgVcycleArgs) *vc = (QCDP(MgVcycleArgs) *)args;
  QDP_Subset sub = vc->sub;
  int nv = vc->nv;
  QDPN(ColorVector) **r = vc->r;
  QDPN(ColorVector) **p = vc->p;
  QDPN(ColorVector) **Ap = vc->Ap;

  QLA_D_Real innrm=0, nrm;
  if(vc->verbose>0) {
    r_eq_norm2_V(&innrm, in);
    printf0("%*sVcycle in2 = %g\n", vc->indent, "", innrm);
  }

  { sub = QDP_all_L(QDPN(get_lattice_V)(r[0])); V_eq_zero(r); sub = vc->sub; }
  if(vc->npre>0) {
    if(vc->gcr) {
      QCDP(gcrSolve)(vc->gcr, out, in, vc->op, vc->opargs, NULL, NULL, FLT_EPSILON, vc->npre, vc->npre, vc->sub);
    } else {
      QCDP(cglsSolve)(vc->cgls, out, in, NULL, vc->sop, vc->sopargs, NULL, NULL, 1, 0, 0, vc->npre, vc->sub, vc->sub2);
    }
    V_eq_V(r, in);
    vc->op(Ap, out, sign, vc->opargs);
    V_meq_V(r, Ap);
    if(vc->verbose>0) {
      r_eq_norm2_V(&nrm, r);
      printf0("%*sVcycle r2 = %g  r2/in2 = %g\n", vc->indent, "", nrm, nrm/innrm);
    }
    { double t1 = QDP_time(); vc->tpre += t1 - t0; t0 = t1; }
    vc->cop(p, r, sign, vc->copargs);
    { double t1 = QDP_time(); vc->tcoarse += t1 - t0; t0 = t1; }
    V_peq_V(out, p);
  } else {
    V_eq_V(r, in);
    { double t1 = QDP_time(); vc->tpre += t1 - t0; t0 = t1; }
    vc->cop(out, r, sign, vc->copargs);
    { double t1 = QDP_time(); vc->tcoarse += t1 - t0; t0 = t1; }
  }

  //if(vc->npost) {
    vc->op(Ap, out, sign, vc->opargs);
    if(vc->s!=1) {
      QLA_D_Complex pAb, z;
      QLA_D_Real pAAp, s=vc->s;
      c_eq_V_dot_V(&pAb, Ap, in);
      r_eq_norm2_V(&pAAp, Ap);
      QLA_c_eq_c_times_r(z, pAb, (1-s)/pAAp);
      QLA_c_peq_r(z, s);
      V_eq_V(r, out);
      V_eq_c_times_V(out, &z, r);
      V_eq_V(r, Ap);
      V_eq_c_times_V(Ap, &z, r);
    }
    V_eq_V(r, in);
    V_meq_V(r, Ap);
  //}

  if(vc->verbose>0) {
    r_eq_norm2_V(&nrm, r);
    printf0("%*sVcycle r2 = %g  r2/in2 = %g\n", vc->indent, "", nrm, nrm/innrm);
    QLA_Real xr, xb;
    r_eq_re_V_dot_V(&xr, out, r);
    r_eq_re_V_dot_V(&xb, out, in);
    xr = -xr-xb;
    printf0(" VcAinv norm = %-13.8g\n", xr);
  }

  if(vc->npost) {
    if(vc->gcr) {
      QCDP(gcrSolve)(vc->gcr, p, r, vc->op, vc->opargs, NULL, NULL, FLT_EPSILON, vc->npost, vc->npost, vc->sub);
    } else {
      int level = (vc->indent/2)-1; // hack
      if(level==0) {
	QCDP(cglsSolve)(vc->cgls, p, r, NULL, vc->sop, vc->sopargs, NULL, NULL, vc->delta, 0, 0, vc->npost, vc->sub, vc->sub2);
      } else {
	QCDP(cgSolve)(vc->cgls, p, r, vc->sop, vc->sopargs, NULL, NULL, 0, vc->npost, vc->sub);
      }
    }
    V_peq_V(out, p);
  }

  if(vc->verbose>0) {
    V_eq_V(r, in);
    vc->op(Ap, out, sign, vc->opargs);
    V_meq_V(r, Ap);
    r_eq_norm2_V(&nrm, r);
    printf0("%*sVcycle r2 = %g  r2/in2 = %g\n", vc->indent, "", nrm, nrm/innrm);
    //QLA_Real xr, xb;
    //r_eq_re_V_dot_V(&xr, out, r);
    //r_eq_re_V_dot_V(&xb, out, in);
    //xr = -xr-xb;
    //printf0(" VcAinv norm = %-13.8g\n", xr);
  }
  { double t1 = QDP_time(); vc->tpost += t1 - t0; t0 = t1; }
  vc->count++;
}

#endif
