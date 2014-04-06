#include <qop_internal.h>

#include <string.h>
#include <math.h>
#include <qdp_fn.h>
#include <qdp_dn.h>
#include <qdp_dfn.h>
#include "solvers.h"

#if QOP_Precision == 'F'
#define QOPP(x) QOP_F_##x
#define QCDP(x) QOP_F_##x
#define QCDPC(x) QOP_F3_##x
#define QDPN(x) QDP_FN_##x
#else
#define QOPP(x) QOP_D_##x
#define QCDP(x) QOP_D_##x
#define QCDPC(x) QOP_D3_##x
#define QDPN(x) QDP_DN_##x
#endif

typedef QDPN(ColorVector) *CVP;

#define create(v) \
  CVP *v; \
  QOP_malloc(v, CVP, nv); \
  create_V(v); \
  bcg->v = v;

#define destroy(v) \
  CVP *v = bcg->v; \
  destroy_V(v); \
  QOP_free(v);

QOP_Bicgstab *
QOP_bicgstabInit(QDP_Lattice *lat, int nv, int nc)
{
  QOP_Bicgstab *bcg;
  QOP_malloc(bcg, QOP_Bicgstab, 1);

  bcg->nv = nv;
  bcg->nc = nc;

  create(r);
  create(p);
  create(t);
  create(v);
  create(Mx);

  bcg->verbose = 0;
  bcg->indent = 0;
  return bcg;
}

void
QOP_bicgstabFree(QOP_Bicgstab *bcg)
{
  int nv = bcg->nv;

  destroy(r);
  destroy(p);
  destroy(t);
  destroy(v);
  destroy(Mx);

  QOP_free(bcg);
}

#define seti(t) if(!strcmp(s,#t)) bcg->t = (int) v;
void
QOP_bicgstabSet(QOP_Bicgstab *bcg, char *s, double v)
{
  seti(verbose);
  seti(indent);
}

int
QOP_bicgstabSolve(QOP_Bicgstab *bcg, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
		   void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		   void *Aargs,
		   void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		   void *Margs,
		   double res, int itnlim, QDP_Subset sub)
{
  return QOP_bicgstabSolveS(bcg, x, b, Aop, Aargs, Mop, Margs, 1, res, itnlim, sub);
}

int
QOP_bicgstabSolveS(QOP_Bicgstab *bcg, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
		    void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		    void *Aargs,
		    void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		    void *Margs,
		    int sign, double res, int itnlim, QDP_Subset sub)
{
  QLA_D_Complex alpha, beta, omega, rho0, rho1, ctmp1, ctmp2;
  QLA_D_Real bsq, rsq, rsqstop, t2;
  int itn = 0;
  int nv = bcg->nv;

  CVP *r = bcg->r;
  CVP *p = bcg->p;
  CVP *t = bcg->t;
  CVP *v = bcg->v;
  CVP *Mx = bcg->Mx;
  CVP *r0 = b;

#define M(a,b) if(Mop) Mop(a,b,sign,Margs); else V_eq_V(a,b);
#define A(a,b) Aop(a,b,sign,Aargs)

  r_eq_norm2_V(&bsq, b);
  if(res<0) {
    rsqstop = res*res;
  } else {
    rsqstop = res*res*bsq;
  }

  V_eq_V(r, b);
  rsq = bsq;
  V_eq_V(p, r);
  V_eq_zero(x);

  QLA_c_eq_r(rho0, 1);
  QLA_c_eq_r(alpha, 1);
  QLA_c_eq_r(omega, 1);
  QLA_c_eq_r(rho1, rsq);

  if(bcg->verbose>0) QOP_printf0("%*s%i rsq = %g\n", bcg->indent, "", itn, rsq/bsq);
  while(rsq>rsqstop && itn<itnlim) {
    if(itn>0) {
      QLA_D_Complex ctmp1, ctmp2;
      QLA_c_eq_c(rho0, rho1);
      c_eq_V_dot_V(&rho1, r0, r);
#if 1
      QLA_c_eq_c_times_c(ctmp1,rho1,alpha);
      QLA_c_eq_c_times_c(ctmp2,rho0,omega);
      if(QLA_norm2_c(ctmp2)==0.) break;
      //QOP_printf0("ctmp2: %g\n", QLA_norm2_c(ctmp2));
      QLA_D_c_eq_c_div_c(beta,ctmp1,ctmp2);
#else
      if(QLA_norm2_c(rho1)==0.) break;
      QLA_c_eq_c_div_c(ctmp1,rho1,rho0);
      QLA_c_eq_c_div_c(ctmp2,alpha,omega);
      QLA_c_eq_c_times_c(beta,ctmp1,ctmp2);
#endif
      V_eq_V(t, p);
      V_meq_c_times_V(t, &omega, v);
      V_eq_V(p, r);
      V_peq_c_times_V(p, &beta, t);
    }
    itn++;

    M(Mx, p);
    A(v, Mx);

    c_eq_V_dot_V(&ctmp1, r0, v);
    if(QLA_norm2_c(ctmp1)==0.) break;
    //QOP_printf0("ctmp1: %g\n", QLA_norm2_c(ctmp1));
    QLA_D_c_eq_c_div_c(alpha, rho1, ctmp1);
    V_meq_c_times_V(r, &alpha, v);

    M(Mx, r);
    A(t, Mx);

    c_eq_V_dot_V(&ctmp2, t, r);
    r_eq_norm2_V(&t2, t);
    //QOP_printf0("t2: %g\n", t2);
    QLA_c_eq_c_div_r(omega, ctmp2, t2);

    V_peq_c_times_V(x, &omega, r);
    V_peq_c_times_V(x, &alpha, p);
    V_meq_c_times_V(r, &omega, t);

    r_eq_norm2_V(&rsq, r);
    if(QLA_norm2_c(omega)==0.) break;
    if(bcg->verbose>1) QOP_printf0("%*s%i rsq = %g\n", bcg->indent, "", itn, rsq/bsq);
  }
  if(bcg->verbose>0) QOP_printf0("%*s%i rsq = %g\n", bcg->indent, "", itn, rsq/bsq);

  return itn;
#undef A
#undef M
}

int
QCDP(bicgstabSolveA)(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
		     int sign, void *args)
{
  QCDP(BicgstabSolveArgs) *sa = (QCDP(BicgstabSolveArgs) *)args;
  return QCDP(bicgstabSolveS)(sa->bicgstab, out, in, sa->op, sa->opargs,
			      sa->pop, sa->popargs, sign, 
			      sa->res, sa->itmax, sa->sub);
}

int
QCDP(bicgstabSolveEo)(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
		      int sign, void *args)
{
  QLA_D_Real rsq, res;
  QCDP(BicgstabSolveArgs) *sa = (QCDP(BicgstabSolveArgs) *)args;
  int nv = sa->bicgstab->nv;
  QDP_Subset sub = QDP_all_L(QDPN(get_lattice_V)(in[0]));
  int nits;
  r_eq_norm2_V(&rsq, in);
  res = -sa->res*sqrt(rsq);
  sa->project(sa->ineo, in, sa->opargs);
  nits = QCDP(bicgstabSolve)(sa->bicgstab, sa->outeo, sa->ineo, sa->op,
			     sa->opargs, sa->pop, sa->popargs,
			     res, sa->itmax, sa->sub);
  sa->reconstruct(out, sa->outeo, in, sa->opargs);
  return nits;
}
