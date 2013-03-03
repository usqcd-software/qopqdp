#include <qop_internal.h>

#define getU(tn, mu, nu) cacheshift(&ftmps[tn][nu], in[tn], nu, QDP_forward, 0)
#define getC(nu) cacheshift(&ctmps[nu], tc, nu, QDP_forward, !(ctn[nu]++))
#define shiftb(x, mu) cacheshift(&b ## x[mu], x, mu, QDP_backward, 1)
static QDP_ColorMatrix *
cacheshift(QDP_ColorMatrix **tmp, QDP_ColorMatrix *in, int mu, QDP_ShiftDir dir, int redo)
{
#define NC QDP_get_nc(in)
  QDP_ColorMatrix *r = *tmp;
  if(r==NULL) {
    r = *tmp = QDP_create_M();
    redo = 1;
  }
  if(redo) {
    QDP_M_eq_sM(r, in, QDP_neighbor[mu], dir, QDP_all);
  }
  return r;
#undef NC
}

// topdir = 1..nd
// sidedir = -nd..nd
// toplinknum,sidelinknum = 0..nin-1
void
QOP_staples(int nout, int nin, QDP_ColorMatrix *out[], QDP_ColorMatrix *in[],
	    int nstaples[], int *topdir[], int *sidedir[],
	    int *toplinknum[], int *sidelinknum[], QLA_Real *coef[])
{
#define NC QDP_get_nc(in[0])
  int nd = QDP_ndim();
  QDP_ColorMatrix *ftmps[nin][nd], *t1, *t2, *bt2[nd];
  for(int i=0; i<nin; i++)
    for(int j=0; j<nd; j++)
      ftmps[i][j] = NULL;
  for(int i=0; i<nd; i++) bt2[i] = NULL;
  t1 = QDP_create_M();
  t2 = QDP_create_M();

  for(int io=0; io<nout; io++) {
    //QOP_printf0("%i: ns: %i\n", io, nstaples[io]);
    for(int s=0; s<nstaples[io]; s++) {
      QLA_Real c = coef[io][s];
      int tn = toplinknum[io][s];
      int sdir = sidedir[io][s];
      //QOP_printf0(" %i:  sdir: %i  c: %g\n", s, sdir, c);
      if(sdir==0) {
	if(c==1) {
	  QDP_M_peq_M(out[io], in[tn], QDP_all);
	} else {
	  QDP_M_peq_r_times_M(out[io], &c, in[tn], QDP_all);
	}
      } else if(sdir>0) {
	int nu = sdir-1;
	int mu = topdir[io][s]-1;
	int sn = sidelinknum[io][s];
	QDP_ColorMatrix *Umunu = getU(tn, mu, nu);
	QDP_ColorMatrix *Unumu = getU(sn, nu, mu);
	QDP_M_eq_M_times_M(t1, in[sn], Umunu, QDP_all);
	if(c==1) {
	  QDP_M_peq_M_times_Ma(out[io], t1, Unumu, QDP_all);
	} else {
	  QDP_M_eq_M_times_Ma(t2, t1, Unumu, QDP_all);
	  QDP_M_peq_r_times_M(out[io], &c, t2, QDP_all);
	}
      } else {
	int nu = -sdir-1;
	int mu = topdir[io][s]-1;
	int sn = sidelinknum[io][s];
	QDP_ColorMatrix *Unumu = getU(sn, nu, mu);
	QDP_M_eq_M_times_M(t1, in[tn], Unumu, QDP_all);
	QDP_M_eq_Ma_times_M(t2, in[sn], t1, QDP_all);
	QDP_ColorMatrix *tb = shiftb(t2, nu);
	if(c==1) {
	  QDP_M_peq_M(out[io], tb, QDP_all);
	} else {
	  QDP_M_peq_r_times_M(out[io], &c, tb, QDP_all);
	}
	QDP_discard_M(tb);
      }
    }
  }

  for(int i=0; i<nin; i++)
    for(int j=0; j<nd; j++)
      if(ftmps[i][j]!=NULL) QDP_destroy_M(ftmps[i][j]);
  for(int i=0; i<nd; i++) if(bt2[i]!=NULL) QDP_destroy_M(bt2[i]);
  QDP_destroy_M(t1);
  QDP_destroy_M(t2);
#undef NC
}

// topdir = 1..nd
// sidedir = -nd..nd
// toplinknum,sidelinknum = 0..nin-1
void
QOP_staples_deriv(int nout, int nin, QDP_ColorMatrix *deriv[],
		  QDP_ColorMatrix *chain[], QDP_ColorMatrix *in[],
		  int nstaples[], int *topdir[], int *sidedir[],
		  int *toplinknum[], int *sidelinknum[], QLA_Real *coef[])
{
#define NC QDP_get_nc(in[0])
  int nd = QDP_ndim();
  QDP_ColorMatrix *ftmps[nin][nd], *t1, *t2, *t3, *t4, *tc, *bt2[nd], *bt3[nd], *ctmps[nd];
  int ctn[nd];
  for(int i=0; i<nin; i++)
    for(int j=0; j<nd; j++)
      ftmps[i][j] = NULL;
  for(int i=0; i<nd; i++) bt2[i] = bt3[i] = ctmps[i] = NULL;
  t1 = QDP_create_M();
  t2 = QDP_create_M();
  t3 = QDP_create_M();
  t4 = QDP_create_M();
  tc = QDP_create_M();

  // process in reverse in case calculated staples used as input for others
  for(int io=nout-1; io>=0; io--) {
    for(int i=0; i<nd; i++) {
      if(ctmps[i]) QDP_discard_M(ctmps[i]);
      ctn[i] = 0;
    }
    QDP_M_eq_M(tc, chain[io], QDP_all);
    for(int s=0; s<nstaples[io]; s++) {
      QLA_Real c = coef[io][s];
      int tn = toplinknum[io][s];
      int sdir = sidedir[io][s];
      //QOP_printf0("io: %i  s: %i  sdir: %i  tn: %i  c: %g\n", io, s, sdir, tn, c);
      if(sdir==0) {
	if(c==1) {
	  QDP_M_peq_M(deriv[tn], tc, QDP_all);
	} else {
	  QDP_M_peq_r_times_M(deriv[tn], &c, tc, QDP_all);
	}
      } else if(sdir>0) {
	int nu = sdir-1;
	int mu = topdir[io][s]-1;
	int sn = sidelinknum[io][s];
	//QOP_printf0("  mu: %i  nu: %i  sn: %i\n", mu, nu, sn);
	QDP_ColorMatrix *Umunu = getU(tn, mu, nu);
	QDP_ColorMatrix *Unumu = getU(sn, nu, mu);
	QDP_M_eq_M_times_M(t1, in[sn], Umunu, QDP_all);
	QDP_M_eq_Ma_times_M(t2, tc, t1, QDP_all);
	QDP_ColorMatrix *tb2 = shiftb(t2, mu);
	QDP_M_eq_M_times_M(t1, tc, Unumu, QDP_all);
	QDP_M_eq_Ma_times_M(t3, in[sn], t1, QDP_all);
	QDP_ColorMatrix *tb3 = shiftb(t3, nu);
	if(c==1) {
	  QDP_M_peq_M_times_Ma(deriv[sn], t1, Umunu, QDP_all);
	  QDP_M_peq_M(deriv[sn], tb2, QDP_all);
	  QDP_M_peq_M(deriv[tn], tb3, QDP_all);
	} else {
	  QDP_M_eq_M_times_Ma(t4, t1, Umunu, QDP_all);
	  QDP_M_peq_r_times_M(deriv[sn], &c, t4, QDP_all);
	  QDP_M_peq_r_times_M(deriv[sn], &c, tb2, QDP_all);
	  QDP_M_peq_r_times_M(deriv[tn], &c, tb3, QDP_all);
	}
	QDP_discard_M(tb2);
	QDP_discard_M(tb3);
      } else {
	int nu = -sdir-1;
	int mu = topdir[io][s]-1;
	int sn = sidelinknum[io][s];
	QDP_ColorMatrix *Cmunu = getC(nu);
	QDP_ColorMatrix *Unumu = getU(sn, nu, mu);
	QDP_M_eq_M_times_M(t1, in[sn], Cmunu, QDP_all);
	QDP_M_eq_Ma_times_M(t2, in[tn], t1, QDP_all);
	QDP_ColorMatrix *tb2 = shiftb(t2, mu);
	QDP_M_eq_M_times_M(t3, in[tn], Unumu, QDP_all);
	if(c==1) {
	  QDP_M_peq_M_times_Ma(deriv[tn], t1, Unumu, QDP_all);
	  QDP_M_peq_M_times_Ma(deriv[sn], t3, Cmunu, QDP_all);
	  QDP_M_peq_M(deriv[sn], tb2, QDP_all);
	} else {
	  QDP_M_eq_M_times_Ma(t4, t1, Unumu, QDP_all);
	  QDP_M_peq_r_times_M(deriv[tn], &c, t4, QDP_all);
	  QDP_M_eq_M_times_Ma(t4, t3, Cmunu, QDP_all);
	  QDP_M_peq_r_times_M(deriv[sn], &c, t4, QDP_all);
	  QDP_M_peq_r_times_M(deriv[sn], &c, tb2, QDP_all);
	}
	QDP_discard_M(tb2);
      }
    }
  }

  for(int i=0; i<nin; i++)
    for(int j=0; j<nd; j++)
      if(ftmps[i][j]!=NULL) QDP_destroy_M(ftmps[i][j]);
  for(int i=0; i<nd; i++) {
    if(bt2[i]!=NULL) QDP_destroy_M(bt2[i]);
    if(bt3[i]!=NULL) QDP_destroy_M(bt3[i]);
    if(ctmps[i]!=NULL) QDP_destroy_M(ctmps[i]);
  }
  QDP_destroy_M(t1);
  QDP_destroy_M(t2);
  QDP_destroy_M(t3);
  QDP_destroy_M(t4);
  QDP_destroy_M(tc);
#undef NC
}
