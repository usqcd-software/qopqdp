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

QOPP(Gcr) *
QCDP(gcrInit)(QDP_Lattice *lat, int nv, int nc, int ngcr)
{
  QOP_Gcr *gcr;
  QOP_malloc(gcr, QOPP(Gcr), 1);

  gcr->restart = 1;
  if(ngcr<0) { gcr->restart = 0; ngcr = -ngcr; }
  ngcr++;
  gcr->lat = lat;
  gcr->nv = nv;
  gcr->nc = nc;
  gcr->ngcr = ngcr;

  CVP *r, **Mr, **AMr;
  QOP_malloc(Mr, CVP *, ngcr);
  QOP_malloc(AMr, CVP *, ngcr);

  QOP_malloc(r, CVP, nv);
  create_V(r);
  for(int i=0; i<ngcr; i++) {
    QOP_malloc(Mr[i], CVP, nv);
    create_V(Mr[i]);
    { QDP_Subset sub=QDP_all_L(lat); V_eq_zero(Mr[i]); }
    QOP_malloc(AMr[i], CVP, nv);
    create_V(AMr[i]);
  }
  gcr->r = r;
  gcr->Mr = Mr;
  gcr->AMr = AMr;
  gcr->verbose = 0;
  gcr->reuse = 0;
  gcr->nsmooth = 0;

  return gcr;
}

void
QCDP(gcrFree)(QOP_Gcr *gcr)
{
  int nv = gcr->nv;
  int ngcr = gcr->ngcr;

  CVP *r = gcr->r;
  CVP **Mr = gcr->Mr;
  CVP **AMr = gcr->AMr;

  destroy_V(r);
  QOP_free(r);
  for(int i=0; i<ngcr; i++) {
    destroy_V(Mr[i]);
    QOP_free(Mr[i]);
    destroy_V(AMr[i]);
    QOP_free(AMr[i]);
  }
  QOP_free(Mr);
  QOP_free(AMr);
  QOP_free(gcr);
}

#define seti(t) if(!strcmp(s,#t)) gcr->t = (int) v;
void
QCDP(gcrSet)(QOP_Gcr *gcr, char *s, double v)
{
  seti(verbose);
  seti(indent);
  seti(reuse);
  seti(nsmooth);
}

static void
update(CVP *x, CVP *r, CVP **Mr, CVP **AMr, QLA_D_Real *AMrn,
       int k, int k0, int ngcr, int nv, QDP_Subset sub1,
       QDP_Subset sub2, double delta, double *Mrn, QLA_D_Complex *alp)
{
  QLA_D_Complex alpha, AMrr;
  QDP_Subset sub = sub1;

#if 1
#if 1
  for(int i=k0; i!=k; i=(i+1)%ngcr) {
    QLA_D_Complex beta;
    if(delta==0) {
      QLA_D_Complex z;
      c_eq_V_dot_V(&z, AMr[i], AMr[k]);
      QLA_c_eq_c_div_r(beta, z, -AMrn[i]);
    } else {
      QLA_D_Complex z1, z2;
      c_eq_V_dot_V(&z2, Mr[i], Mr[k]);
      sub=sub2; c_eq_V_dot_V(&z1, AMr[i], AMr[k]); sub=sub1;
      QLA_c_peq_r_times_c(z1, delta, z2);
      QLA_c_eq_c_div_r(beta, z1, -(AMrn[i]+delta*Mrn[i]));
    }
    V_peq_c_times_V(Mr[k], &beta, Mr[i]);
    sub=sub2; V_peq_c_times_V(AMr[k], &beta, AMr[i]); sub=sub1;
  }
#endif
#else
  QLA_D_Complex zv[nv*ngcr];
  CVP v1[nv*ngcr], v2[nv*ngcr];
  int l = 0;
  for(int i=k0; i!=k; i=(i+1)%ngcr) {
    for(int j=0; j<nv; j++) {
      v1[l] = AMr[i][j];
      v2[l] = AMr[k][j];
      l++;
    }
  }
  if(l>0) {
    QDPN(c_veq_V_dot_V)(zv, v1, v2, sub, l);
    l = 0;
    for(int i=k0; i!=k; i=(i+1)%ngcr) {
      QLA_D_Complex beta, z;
      QLA_c_eq_r(z, 0);
      for(int j=0; j<nv; j++) QLA_c_peq_c(z, zv[l++]);
      QLA_c_eq_c_div_r(beta, z, -AMrn[i]);
      V_peq_c_times_V(Mr[k], &beta, Mr[i]);
      V_peq_c_times_V(AMr[k], &beta, AMr[i]);
    }
  }
#endif

  if(delta==0) {
    r_eq_norm2_V(&AMrn[k], AMr[k]);
    c_eq_V_dot_V(&AMrr, AMr[k], r);
    QLA_c_eq_c_div_r(alpha, AMrr, AMrn[k]);
    V_peq_c_times_V(x, &alpha, Mr[k]);
    V_meq_c_times_V(r, &alpha, AMr[k]);
  } else {
    r_eq_norm2_V(&Mrn[k], Mr[k]);
    sub=sub2; r_eq_norm2_V(&AMrn[k], AMr[k]); sub=sub1;
    c_eq_V_dot_V(&AMrr, Mr[k], r);
    QLA_c_eq_c_div_r(alpha, AMrr, AMrn[k]+delta*Mrn[k]);
    V_peq_c_times_V(x, &alpha, Mr[k]);
    //V_meq_c_times_V(r, &alpha, AMr[k]);
    *alp = alpha;
  }
}

int
QCDP(gcrSolve)(QOP_Gcr *gcr, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
	       void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[],
			int sign, void *args),
	       void *Aargs,
	       void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[],
			int sign, void *args),
	       void *Margs,
	       double res, int itnlim, int ngcr, QDP_Subset sub)
{
  int nv = gcr->nv;
  //int restart = gcr->restart;
  //restart = 1;
  if(ngcr<0) { /*restart = 0;*/ ngcr = -ngcr; }
  ngcr++;
  if(ngcr>gcr->ngcr) {
    //printf0("error: ngcr (%i) > gcr->ncgr (%i)\n", ngcr, gcr->ngcr);
    //QDP_abort(1);
    ngcr = gcr->ngcr;
  }
  QLA_D_Real bsq, rsq, rsqstop, AMrn[ngcr];
  int itn = 0;
  int k = -1;
  int k0 = 0;

  CVP *r = gcr->r;
  CVP **Mr = gcr->Mr;
  CVP **AMr = gcr->AMr;

#define M(a,b) if(Mop) Mop(a,b,1,Margs); else V_eq_V(a,b);
#define A(a,b) Aop(a,b,1,Aargs)

  r_eq_norm2_V(&bsq, b);
  if(res>0) {
    rsqstop = res*res*bsq;
  } else {
    rsqstop = res*res;
  }
  //printf0("bsq = %g\n", bsq);
  //printf0("res = %g\n", res);
  //printf0("rsqstop = %g\n", rsqstop);

  V_eq_zero(x);
  V_eq_V(r, b);

  if(gcr->reuse) {
    for(k=0; k<ngcr; k++) {
      QLA_D_Real nrm;
      r_eq_norm2_V(&nrm, Mr[k]);
      if(nrm==0) break;
      A(AMr[k], Mr[k]);
      update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);
      if(gcr->verbose>1) {
	r_eq_norm2_V(&nrm, r);
	QOP_printf0("GCR reuse %i  rsq = %g\n", k, nrm);
      }
    }
    k--;
  }
  if(k<0) {
    k++;
    M(Mr[k], r);
    A(AMr[k], Mr[k]);
    update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);
    itn++;
  }

  r_eq_norm2_V(&rsq, r);
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  while(rsq>rsqstop && itn<itnlim) {
    k = (k+1)%ngcr;
    if(k==k0) k0 = (k0+1)%ngcr;

    if(itn%(gcr->nsmooth+1)==0) {
      M(Mr[k], r);
    } else {
      int km1 = (k-1+ngcr)%ngcr;
      //Aop(Mr[k],AMr[km1],-1,Aargs);
      //V_eq_V(Mr[k], r);
      V_eq_V(Mr[k], AMr[km1]);
    }
    A(AMr[k], Mr[k]);
    itn++;

    update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);

#if 0
    if(gcr->verbose>1) {
      static QLA_Real aino=0;
      QLA_Real xr, xb;
      r_eq_re_V_dot_V(&xr, x, r);
      r_eq_re_V_dot_V(&xb, x, b);
      xr = -xr-xb;
      QOP_printf0(" Ainv norm = %-13.8g    %g\n", xr, xr-aino);
      aino = xr;
    }
#endif

    r_eq_norm2_V(&rsq, r);
    if(gcr->verbose>1) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  }
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);

  //r_eq_norm2_V(&bsq, x);
  //QOP_printf0(" %i xsq = %g\n", itn, bsq);

  if(gcr->reuse) {
    k = (k+1)%ngcr;
    V_eq_V(Mr[k], x);
  }

  return itn;
#undef A
#undef M
}

// solves   (A' A + delta) x = b1 + A' b2
// with r2 = b2 - A x
// and  r1 = b1 + A' r2 - delta x 
int
QCDP(gcrlsSolve)(QOP_Gcr *gcr, QDPN(ColorVector) *x[],
		 QDPN(ColorVector) *b1[], QDPN(ColorVector) *b2[],
		 void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		 void *Aargs,
		 void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		 void *Margs,
		 double delta, double res1, double res2, int itnlim, int ngcr, QDP_Subset sub1, QDP_Subset sub2)
{
  int nv = gcr->nv;
  int nc = gcr->nc;
  QDP_Lattice *lat = gcr->lat;
  //int restart = gcr->restart;
  QDP_Subset sub = sub1;
  //restart = 1;
  if(ngcr<0) { /*restart = 0;*/ ngcr = -ngcr; }
  ngcr++;
  if(ngcr>gcr->ngcr) {
    //QOP_printf0("error: ngcr (%i) > gcr->ncgr (%i)\n", ngcr, gcr->ngcr);
    //QDP_abort(1);
    ngcr = gcr->ngcr;
  }
  QLA_D_Complex alpha;
  QLA_D_Real bsq, rsq, rsqstop, AMrn[ngcr], Mrn[ngcr];
  int itn = 0;
  int k = -1;
  int k0 = 0;

  CVP *r = gcr->r;
  CVP **Mr = gcr->Mr;
  CVP **AMr = gcr->AMr;
  CVP AAMr[nv];
  create_V(AAMr);

#define M(a,b) if(Mop) Mop(a,b,1,Margs); else V_eq_V(a,b);
#define A(a,b) Aop(a,b,1,Aargs)

  // assume b2==0 for now
  r_eq_norm2_V(&bsq, b1);
  if(res1>0) {
    rsqstop = res1*res1*bsq;
  } else {
    rsqstop = res1*res1;
  }
  //QOP_printf0("bsq = %g\n", bsq);
  //QOP_printf0("res = %g\n", res);
  //QOP_printf0("rsqstop = %g\n", rsqstop);

  V_eq_zero(x);
  V_eq_V(r, b1);

  if(gcr->reuse) {
    for(k=0; k<ngcr; k++) {
      QLA_D_Real nrm;
      r_eq_norm2_V(&nrm, Mr[k]);
      if(nrm==0) break;
      A(AMr[k], Mr[k]);
      update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub1, sub2, delta, Mrn, &alpha);
      Aop(AAMr, AMr[k], -1, Aargs);
      V_peq_r_times_V(AAMr, &delta, Mr[k]);
      V_meq_c_times_V(r, &alpha, AAMr);
      if(gcr->verbose>1) {
	r_eq_norm2_V(&nrm, r);
	QOP_printf0("GCR reuse %i  rsq = %g\n", k, nrm);
      }
    }
    k--;
  }
  if(k<0) {
    k++;
    M(Mr[k], r);
    A(AMr[k], Mr[k]);
    update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub1, sub2, delta, Mrn, &alpha);
    Aop(AAMr, AMr[k], -1, Aargs);
    V_peq_r_times_V(AAMr, &delta, Mr[k]);
    V_meq_c_times_V(r, &alpha, AAMr);
    itn++;
  }

  r_eq_norm2_V(&rsq, r);
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  while(rsq>rsqstop && itn<itnlim) {
    k = (k+1)%ngcr;
    if(k==k0) k0 = (k0+1)%ngcr;

    if(itn%(gcr->nsmooth+1)==0) {
      M(Mr[k], r);
    } else {
      //int km1 = (k-1+ngcr)%ngcr;
      //Aop(Mr[k],AMr[km1],-1,Aargs);
      V_eq_V(Mr[k], r);
      //V_eq_V(Mr[k], AMr[km1]);
      //V_eq_V(Mr[k], AAMr);
    }
    A(AMr[k], Mr[k]);
    itn++;
    update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub1, sub2, delta, Mrn, &alpha);
    //QLA_c_eq_r(alpha, 1);
    //V_peq_c_times_V(x, &alpha, Mr[k]);
    Aop(AAMr, AMr[k], -1, Aargs);
    V_peq_r_times_V(AAMr, &delta, Mr[k]);
    V_meq_c_times_V(r, &alpha, AAMr);
    if(gcr->verbose>1) {
      static QLA_Real aino=0;
      QOP_printf0(" alpha = %g\n", QLA_norm2_c(alpha));
      QLA_Real xr, xb;
      r_eq_re_V_dot_V(&xr, x, r);
      r_eq_re_V_dot_V(&xb, x, b1);
      xr = -xr-xb;
      QOP_printf0(" Ainv norm = %-13.8g    %g\n", xr, xr-aino);
      aino = xr;
    }

    r_eq_norm2_V(&rsq, r);
    if(gcr->verbose>1) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  }
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);

  //r_eq_norm2_V(&bsq, x);
  //QOP_printf0(" %i xsq = %g\n", itn, bsq);

  if(gcr->reuse) {
    k = (k+1)%ngcr;
    V_eq_V(Mr[k], x);
  }
  destroy_V(AAMr);

  return itn;
#undef A
#undef M
}

static void
update2(CVP *x, CVP *r, CVP **Mr, CVP **AMr, QLA_D_Real *AMrn,
	int k, int k0, int ngcr, int nv, QDP_Subset sub1,
	QDP_Subset sub2, double delta, double *Mrn, QLA_D_Complex *alp)
{
  QLA_D_Complex alpha, AMrr;
  QDP_Subset sub = sub1;

#if 1
#if 1
  for(int i=k0; i!=k; i=(i+1)%ngcr) {
    QLA_D_Complex beta;
    if(delta==0) {
      QLA_D_Complex z;
      c_eq_V_dot_V(&z, AMr[i], Mr[k]);
      QLA_c_eq_c_div_r(beta, z, -AMrn[i]);
    } else {
      QLA_D_Complex z1, z2;
      c_eq_V_dot_V(&z2, Mr[i], Mr[k]);
      sub=sub2; c_eq_V_dot_V(&z1, AMr[i], AMr[k]); sub=sub1;
      QLA_c_peq_r_times_c(z1, delta, z2);
      QLA_c_eq_c_div_r(beta, z1, -(AMrn[i]+delta*Mrn[i]));
    }
    V_peq_c_times_V(Mr[k], &beta, Mr[i]);
    sub=sub2; V_peq_c_times_V(AMr[k], &beta, AMr[i]); sub=sub1;
  }
#endif
#else
  QLA_D_Complex zv[nv*ngcr];
  CVP v1[nv*ngcr], v2[nv*ngcr];
  int l = 0;
  for(int i=k0; i!=k; i=(i+1)%ngcr) {
    for(int j=0; j<nv; j++) {
      v1[l] = AMr[i][j];
      v2[l] = AMr[k][j];
      l++;
    }
  }
  if(l>0) {
    QDPN(c_veq_V_dot_V)(zv, v1, v2, sub, l);
    l = 0;
    for(int i=k0; i!=k; i=(i+1)%ngcr) {
      QLA_D_Complex beta, z;
      QLA_c_eq_r(z, 0);
      for(int j=0; j<nv; j++) QLA_c_peq_c(z, zv[l++]);
      QLA_c_eq_c_div_r(beta, z, -AMrn[i]);
      V_peq_c_times_V(Mr[k], &beta, Mr[i]);
      V_peq_c_times_V(AMr[k], &beta, AMr[i]);
    }
  }
#endif

  if(delta==0) {
    r_eq_re_V_dot_V(&AMrn[k], Mr[k], AMr[k]);
    c_eq_V_dot_V(&AMrr, Mr[k], r);
    QLA_c_eq_c_div_r(alpha, AMrr, AMrn[k]);
    V_peq_c_times_V(x, &alpha, Mr[k]);
    V_meq_c_times_V(r, &alpha, AMr[k]);
  } else {
    r_eq_norm2_V(&Mrn[k], Mr[k]);
    sub=sub2; r_eq_norm2_V(&AMrn[k], AMr[k]); sub=sub1;
    c_eq_V_dot_V(&AMrr, Mr[k], r);
    QLA_c_eq_c_div_r(alpha, AMrr, AMrn[k]+delta*Mrn[k]);
    V_peq_c_times_V(x, &alpha, Mr[k]);
    //V_meq_c_times_V(r, &alpha, AMr[k]);
    *alp = alpha;
  }
}

int
QCDP(gcrcgSolve)(QOP_Gcr *gcr, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
		 void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		 void *Aargs,
		 void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[], int sign, void *args),
		 void *Margs,
		 double res, int itnlim, int ngcr, QDP_Subset sub)
{
  int nv = gcr->nv;
  //int restart = gcr->restart;
  //restart = 1;
  if(ngcr<0) { /*restart = 0;*/ ngcr = -ngcr; }
  ngcr++;
  if(ngcr>gcr->ngcr) {
    //QOP_printf0("error: ngcr (%i) > gcr->ncgr (%i)\n", ngcr, gcr->ngcr);
    //QDP_abort(1);
    ngcr = gcr->ngcr;
  }
  QLA_D_Real bsq, rsq, rsqstop, AMrn[ngcr];
  int itn = 0;
  int k = -1;
  int k0 = 0;

  CVP *r = gcr->r;
  CVP **Mr = gcr->Mr;
  CVP **AMr = gcr->AMr;

#define M(a,b) if(Mop) Mop(a,b,1,Margs); else V_eq_V(a,b);
#define A(a,b) Aop(a,b,1,Aargs)

  r_eq_norm2_V(&bsq, b);
  if(res>0) {
    rsqstop = res*res*bsq;
  } else {
    rsqstop = res*res;
  }
  //QOP_printf0("bsq = %g\n", bsq);
  //QOP_printf0("res = %g\n", res);
  //QOP_printf0("rsqstop = %g\n", rsqstop);

  V_eq_zero(x);
  V_eq_V(r, b);

  if(gcr->reuse) {
    for(k=0; k<ngcr; k++) {
      QLA_D_Real nrm;
      r_eq_norm2_V(&nrm, Mr[k]);
      if(nrm==0) break;
      A(AMr[k], Mr[k]);
      update2(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);
      if(gcr->verbose>1) {
	r_eq_norm2_V(&nrm, r);
	QOP_printf0("GCR reuse %i  rsq = %g\n", k, nrm);
      }
    }
    k--;
  }
  if(k<0) {
    k++;
    M(Mr[k], r);
    A(AMr[k], Mr[k]);
    update2(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);
    itn++;
  }

  r_eq_norm2_V(&rsq, r);
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  while(rsq>rsqstop && itn<itnlim) {
    itn++;
    k = (k+1)%ngcr;
    if(k==k0) k0 = (k0+1)%ngcr;

    M(Mr[k], r);
    A(AMr[k], Mr[k]);
    update2(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);

    r_eq_norm2_V(&rsq, r);
    if(gcr->verbose>1) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  }
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);

  //r_eq_norm2_V(&bsq, x);
  //QOP_printf0(" %i xsq = %g\n", itn, bsq);

  if(gcr->reuse) {
    k = (k+1)%ngcr;
    V_eq_V(Mr[k], x);
  }

  return itn;
#undef A
#undef M
}

#define V_eq_zeroe(x) { sub = sube; V_eq_zero(x); sub = suba; }
#define V_eq_zeroo(x) { sub = subo; V_eq_zero(x); sub = suba; }
#define V_eq_Ve(x,y) { sub = sube; V_eq_V(x,y); sub = suba; }
#define V_eq_Vo(x,y) { sub = subo; V_eq_V(x,y); sub = suba; }

int
QCDP(gcrSolve2)(QOP_Gcr *gcr, QDPN(ColorVector) *x[], QDPN(ColorVector) *b[],
		void Aop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[],
			 QOP_evenodd_t eo, void *args),
		void *Aargs,
		void Mop(QDPN(ColorVector) *Ax[], QDPN(ColorVector) *x[],
			 int sign, void *args),
		void *Margs,
		double res, int itnlim, int ngcr)
{
  QDP_Subset suba = QDP_all_L(gcr->lat);
  QDP_Subset sube = QDP_even_L(gcr->lat);
  QDP_Subset subo = QDP_odd_L(gcr->lat);
  QDP_Subset sub = suba;
  int nv = gcr->nv;
  //int restart = gcr->restart;
  //restart = 1;
  if(ngcr<0) { /*restart = 0;*/ ngcr = -ngcr; }
  ngcr++;
  ngcr *= 2;
  if(ngcr>gcr->ngcr) {
    //QOP_printf0("error: ngcr (%i) > gcr->ncgr (%i)\n", ngcr, gcr->ngcr);
    //QDP_abort(1);
    ngcr = gcr->ngcr;
    ngcr -= ngcr&1;
    if(ngcr==0) {
      QOP_printf0("error: ngcr (%i) == 0 (gcr->ngcr = %i)\n", ngcr, gcr->ngcr);
      QDP_abort(1);
    }
  }
  QLA_D_Real bsq, rsq, rsqstop, AMrn[ngcr];
  int itn = 0;
  int k = 0;
  int k0 = 0;

  CVP *r = gcr->r;
  CVP **Mr = gcr->Mr;
  CVP **AMr = gcr->AMr;

#define M(ae,ao,b) if(Mop) { Mop(ae,b,1,Margs); V_eq_Vo(ao,ae); V_eq_zeroo(ae); V_eq_zeroe(ao); } \
  else { V_eq_Ve(ae,b); V_eq_Vo(ao,b); V_eq_zeroo(ae); V_eq_zeroe(ao); }
#define Ae(a,b) Aop(a,b,QOP_EVEN,Aargs)
#define Ao(a,b) Aop(a,b,QOP_ODD,Aargs)

  r_eq_norm2_V(&bsq, b);
  if(res>0) {
    rsqstop = res*res*bsq;
  } else {
    rsqstop = res*res;
  }

  //QOP_printf0("bsq = %g\n", bsq);
  V_eq_zero(x);
  V_eq_V(r, b);

  M(Mr[k], Mr[k+1], r);
  Ae(AMr[k], Mr[k]);
  update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);
  Ao(AMr[k+1], Mr[k+1]);
  update(x, r, Mr, AMr, AMrn, k+1, k0, ngcr, nv, sub, sub, 0, NULL, NULL);

  r_eq_norm2_V(&rsq, r);
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  while(rsq>rsqstop && itn<itnlim) {
    itn++;
    k = (k+2)%ngcr;
    if(k==k0) k0 = (k0+2)%ngcr;

    M(Mr[k], Mr[k+1], r);

    Ae(AMr[k], Mr[k]);
    update(x, r, Mr, AMr, AMrn, k, k0, ngcr, nv, sub, sub, 0, NULL, NULL);

    Ao(AMr[k+1], Mr[k+1]);
    update(x, r, Mr, AMr, AMrn, k+1, k0, ngcr, nv, sub, sub, 0, NULL, NULL);

    r_eq_norm2_V(&rsq, r);
    if(gcr->verbose>1) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);
  }
  if(gcr->verbose>0) QOP_printf0("%*s%i rsq = %g\n", gcr->indent, "", itn, rsq/bsq);

  //r_eq_norm2_V(&bsq, x);
  //QOP_printf0(" %i xsq = %g\n", itn, bsq);

  return itn;
#undef A
#undef M
}

int
QCDP(gcrSolveA)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  QCDP(GcrSolveArgs) *sa = (QCDP(GcrSolveArgs) *)args;
  return QCDP(gcrSolve)(sa->gcr, out, in, sa->op, sa->opargs, sa->pop,
			sa->popargs, sa->res, sa->itmax, sa->ngcr, sa->sub);
}

int
QCDP(gcrcgSolveA)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  QCDP(GcrSolveArgs) *sa = (QCDP(GcrSolveArgs) *)args;
  return QCDP(gcrcgSolve)(sa->gcr, out, in, sa->op, sa->opargs, sa->pop,
			  sa->popargs, sa->res, sa->itmax, sa->ngcr, sa->sub);
}

int
QCDP(gcrSolveEo)(QDPN(ColorVector) **out, QDPN(ColorVector) **in, int sign, void *args)
{
  QCDP(GcrSolveArgs) *sa = (QCDP(GcrSolveArgs) *)args;
  sa->project(sa->ineo, in, sa->opargs);
  int n;
  n = QCDP(gcrSolve)(sa->gcr, sa->outeo, sa->ineo, sa->op, sa->opargs, sa->pop,
		     sa->popargs, sa->res, sa->itmax, sa->ngcr, sa->sub);
  sa->reconstruct(out, sa->outeo, in, sa->opargs);
  return n;
}

#endif
