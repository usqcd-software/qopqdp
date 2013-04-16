#include <qop_internal.h>

#ifdef USE_MG

#include <string.h>
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
  cgls->v = v;

#define destroy(v) \
  CVP *v = cgls->v; \
  destroy_V(v); \
  QOP_free(v);

QOP_Cgls *
QOP_cglsInit(QDP_Lattice *lat, int nv, int nc)
{
  QOP_Cgls *cgls;
  QOP_malloc(cgls, QOP_Cgls, 1);

  cgls->nv = nv;
  cgls->nc = nc;

  create(r1);
  create(r2);
  create(z);
  create(p);
  create(Ap);
  create(Adr2);

  cgls->verbose = 0;
  cgls->indent = 0;
  return cgls;
}

void
QOP_cglsFree(QOP_Cgls *cgls)
{
  int nv = cgls->nv;

  destroy(r1);
  destroy(r2);
  destroy(z);
  destroy(p);
  destroy(Ap);
  destroy(Adr2);

  QOP_free(cgls);
}

#define seti(t) if(!strcmp(s,#t)) cgls->t = (int) v;
void
QOP_cglsSet(QOP_Cgls *cgls, char *s, double v)
{
  seti(verbose);
  seti(indent);
}

// solves   (A' A + delta) x = b1 + A' b2
// with r2 = b2 - A x
// and  r1 = b1 + A' r2 - delta x 
int
QOP_cglsSolve(QOP_Cgls *cgls, QDPN(ColorVector) *x[],
	      QDPN(ColorVector) *b1[], QDPN(ColorVector) *b2[],
	      void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
	      void *Aargs,
	      void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
	      void *Margs,
	      double delta, double res1, double res2, int itnlim, QDP_Subset sub1, QDP_Subset sub2)
{
  QLA_D_Real b1sq, b2sq, r1sq, r2sq, r1sqstop, r2sqstop, alpha, beta, p2=0., Ap2;
  int itn = 0;
  int nv = cgls->nv;
  QDP_Subset sub = sub1;

  CVP *r1 = cgls->r1;
  CVP *r2 = cgls->r2;
  CVP *z = cgls->z;
  CVP *p = cgls->p;
  CVP *Ap = cgls->Ap;
  CVP *Adr2 = cgls->Adr2;

#define M(a,b) if(Mop) Mop(a,b,1,Margs); else V_eq_V(a,b);
#define A(a,b) Aop(a,b,1,Aargs)
#define Ad(a,b) Aop(a,b,-1,Aargs)

  V_eq_zero(x);
  V_eq_zero(p);
  if(b2) {
    sub=sub2; V_eq_V(r2, b2); sub=sub1;
    if(b1) {
      Ad(Adr2, r2);
      V_eq_V_plus_V(r1, b1, Adr2);
    } else {
      Ad(r1, r2);
    }
  } else {
    sub=sub2; V_eq_zero(r2); sub=sub1;
    V_eq_V(r1, b1);
  }

  r_eq_norm2_V(&r1sq, r1);
  sub=sub2; r_eq_norm2_V(&r2sq, r2); sub=sub1;
  b1sq = r1sq;
  if(b1sq==0.) b1sq = 1.;
  b2sq = r2sq;
  if(b2sq==0.) b2sq = 1.;
  r1sqstop = res1*res1*b1sq;
  r2sqstop = res2*res2*b2sq;

  if(cgls->verbose>0) QOP_printf0("%*s%-3i r1sq = %-12g  r2sq = %g\n", cgls->indent, "", itn, r1sq/b1sq, r2sq/b2sq);
  do {
    itn++;
    M(z, r1);
    r_eq_re_V_dot_V(&beta, r1, z);
    beta = 1/beta;
    V_peq_r_times_V(p, &beta, z);
    //r_eq_re_V_dot_V(&beta, p, r1);
    //printf("%g\n", beta);
    A(Ap, p);
    if(delta!=0.) { r_eq_norm2_V(&p2, p); }
    sub=sub2; r_eq_norm2_V(&Ap2, Ap); sub=sub1;
    alpha = 1/(Ap2 + delta*p2);
    //QOP_printf0("alpha = %g\n", alpha);
    V_peq_r_times_V(x, &alpha, p);
    sub=sub2; V_meq_r_times_V(r2, &alpha, Ap); sub=sub1;
    if(b1) {
      Ad(Adr2, r2);
      V_eq_V_plus_V(r1, b1, Adr2);
    } else {
      Ad(r1, r2);
    }
    if(delta!=0.) V_meq_r_times_V(r1, &delta, x);
#if 0
    QLA_Complex zz;
    c_eq_V_dot_V(&zz, p, r1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
    c_eq_V_dot_V(&zz, p, b1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
#endif

    if(r1sqstop!=0.) { r_eq_norm2_V(&r1sq, r1); }
    if(r2sqstop!=0.) { sub=sub2; r_eq_norm2_V(&r2sq, r2); sub=sub1; }

    if(cgls->verbose>1) {
      if(b1) {
	static QLA_Real aino=0;
	QLA_Real xr, xb;
	r_eq_re_V_dot_V(&xr, x, r1);
	r_eq_re_V_dot_V(&xb, x, b1);
	xr = -xr-xb;
	QOP_printf0(" Ainv norm = %-13.8g    %g\n", xr, xr-aino);
	aino = xr;
      }
    }
    if(cgls->verbose>1) QOP_printf0("%*s%-3i r1sq = %-12g  r2sq = %g\n", cgls->indent, "", itn, r1sq/b1sq, r2sq/b2sq);
  } while(r1sq>=r1sqstop && r2sq>=r2sqstop && itn<itnlim);
  if(cgls->verbose>0) QOP_printf0("%*s%-3i r1sq = %-12g  r2sq = %g\n", cgls->indent, "", itn, r1sq/b1sq, r2sq/b2sq);

  return itn;
#undef A
#undef Ad
#undef M
}

// solves  A x = b  for HPD A
int
QOP_cgSolve(QOP_Cgls *cgls, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
	    void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
	    void *Aargs,
	    void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
	    void *Margs,
	    double res, int itnlim, QDP_Subset sub)
{
  QLA_D_Real bsq, rsq, rsqstop, alpha, beta, Ap2;
  int itn = 0;
  int nv = cgls->nv;

  CVP *r = cgls->r1;
  CVP *z = cgls->z;
  CVP *p = cgls->p;
  CVP *Ap = cgls->Ap;

#define A(a,b) Aop(a,b,1,Aargs)

  V_eq_zero(x);
  V_eq_zero(p);
  V_eq_V(r, b);

  r_eq_norm2_V(&rsq, r);
  bsq = rsq;
  rsqstop = res*res*bsq;

  if(cgls->verbose>0) QOP_printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);
  do {
    itn++;
    if(Mop) {
      Mop(z, r, 1, Margs);
      r_eq_re_V_dot_V(&beta, r, z);
      beta = 1/beta;
      V_peq_r_times_V(p, &beta, z);
    } else {
      beta = 1/rsq;
      V_peq_r_times_V(p, &beta, r);
    }
    //r_eq_re_V_dot_V(&beta, p, r1);
    //printf("%g\n", beta);
    A(Ap, p);
    r_eq_norm2_V(&Ap2, Ap);
    alpha = 1/Ap2;
    //QOP_printf0("alpha = %g\n", alpha);
    V_peq_r_times_V(x, &alpha, p);
    V_meq_r_times_V(r, &alpha, Ap);
#if 0
    QLA_Complex zz;
    c_eq_V_dot_V(&zz, p, r1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
    c_eq_V_dot_V(&zz, p, b1);
    QOP_printf0("%g %g\n", QLA_real(zz), QLA_imag(zz));
#endif

    r_eq_norm2_V(&rsq, r);

    if(cgls->verbose>1) QOP_printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);
  } while(rsq>=rsqstop && itn<itnlim);
  if(cgls->verbose>0) QOP_printf0("%*s%-3i rsq = %-12g\n", cgls->indent, "", itn, rsq/bsq);

  return itn;
#undef A
}

int
QOP_cglsSolveA(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
	       int sign, void *args)
{
  QOP_CglsSolveArgs *sa = (QOP_CglsSolveArgs *)args;
  return QOP_cglsSolve(sa->cgls, out, NULL, in, sa->op, sa->opargs, sa->pop,
		       sa->popargs,  0, 0, sa->res, sa->itmax,
		       sa->sub1, sa->sub2);
}

int
QOP_cglsSolveA1(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
		int sign, void *args)
{
  QOP_CglsSolveArgs *sa = (QOP_CglsSolveArgs *)args;
  return QOP_cglsSolve(sa->cgls, out, in, NULL, sa->op, sa->opargs, sa->pop,
		       sa->popargs, sa->delta, sa->res, 0, sa->itmax,
		       sa->sub1, sa->sub2);
}

int
QOP_cglsSolveEo(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
		int sign, void *args)
{
  QOP_CglsSolveArgs *sa = (QOP_CglsSolveArgs *)args;
  sa->project(sa->ineo, in, sa->opargs);
  int n= QOP_cglsSolve(sa->cgls, sa->outeo, NULL, sa->ineo, sa->op, sa->opargs,
			sa->pop, sa->popargs, 0, 0, sa->res, sa->itmax,
			sa->sub1, sa->sub2);
  sa->reconstruct(out, sa->outeo, in, sa->opargs);
  return n;
}

int
QOP_cglsSolveEo1(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
		 int sign, void *args)
{
  QOP_CglsSolveArgs *sa = (QOP_CglsSolveArgs *)args;
  sa->project(sa->ineo, in, sa->opargs);
  int n= QOP_cglsSolve(sa->cgls, sa->outeo, sa->ineo, NULL, sa->op, sa->opargs,
		       sa->pop, sa->popargs, sa->delta, sa->res, 0, sa->itmax,
		       sa->sub1, sa->sub2);
  sa->reconstruct(out, sa->outeo, in, sa->opargs);
  return n;
}

int
QOP_cgSolveA(QDPN(ColorVector) **out, QDPN(ColorVector) **in,
	     int sign, void *args)
{
  QOP_CglsSolveArgs *sa = (QOP_CglsSolveArgs *)args;
  return QOP_cgSolve(sa->cgls, out, in, sa->op, sa->opargs, sa->pop,
		     sa->popargs, sa->res, sa->itmax, sa->sub1);
}


#endif
