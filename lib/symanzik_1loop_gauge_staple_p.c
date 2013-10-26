#include <qop_internal.h>

static QDP_ColorMatrix *sm[2][8], *tm[2]={NULL,NULL};

static int
neighsubeo(int subl, int dir)
{
  return 1 - subl;
}

static int
neighsub32(int subl, int dir)
{
#if 0
  int i, s=0, r=0;
  for(i=QDP_ndim()-1; i>=0; i--) {
    s = 2*s + (coords[i]%2);
    r += coords[i]/2;
  }
  return 16*(r%2) + s;
#endif
  int bitnum = abs(dir)-1;
  int bitmask = 1<<bitnum;
  if( ((subl&bitmask)&&(dir>0)) || (!(subl&bitmask)&&(dir<0)) ) {
    subl ^= 16;
  }
  subl ^= bitmask;
  return subl;
}

//#define eosub(s) subset[s]
//#define eosub(s) QDP_even_and_odd[(s+(s>>1)+(s>>2)+(s>>3))&1]
#define eosub(s) QDP_all_L(lat)

static void
path_prod(QDP_ColorMatrix *u[], QDP_ColorMatrix *m, int path[], int len,
	  int subl, QDP_Subset subset[], int (*neighsubl)(int subl, int dir))
{
  int dir, sn;
  QDP_ShiftDir fb;
  QDP_ColorMatrix *p=NULL, *s=NULL;
  QDP_Lattice *lat = QDP_get_lattice_M(m);

  sn = 0;
  for(int i=0; i<len; i++) {
    dir = abs(path[i])-1;
    // if the path moves in the + dir then we shift from the backward dir
    fb = path[i]<0 ? QDP_forward : QDP_backward;
    if(fb==QDP_backward) { // path is moving in + dir
      if(i==0) {
	QDP_M_eq_Ma(tm[sn], u[dir], subset[subl]);
      } else {
	QDP_M_eq_Ma_times_M(tm[sn], u[dir], p, subset[subl]);
	QDP_discard_M(p);
      }
      subl = neighsubl(subl, path[i]);
      s = sm[sn][4+dir];
      QDP_discard_M(s);
      QDP_M_eq_sM(s, tm[sn], QDP_neighbor_L(lat)[dir], fb, eosub(subl));
      //p = t1; t1 = t2; t2 = p;
      sn = 1-sn;
      p = s;
    } else {
      if(i==0) {
	subl = neighsubl(subl, path[i]);
 	QDP_M_eq_M(tm[1-sn], u[dir], subset[subl]);
      } else {
	QDP_M_eq_M(tm[sn], p, subset[subl]);
	QDP_discard_M(p);
	subl = neighsubl(subl, path[i]);
	s = sm[sn][dir];
	QDP_discard_M(s);
	QDP_M_eq_sM(s, tm[sn], QDP_neighbor_L(lat)[dir], fb, eosub(subl));
	QDP_M_eq_M_times_M(tm[1-sn], u[dir], s, subset[subl]);
	QDP_discard_M(s);
      }
      p = tm[1-sn];
    }
  }
  QDP_M_eq_M(m, p, subset[subl]);
  QDP_discard_M(p);
  QDP_discard_M(s);
}

static void
get_staple_plaq(QDP_ColorMatrix *staple, int mu, QDP_ColorMatrix *u[],
		QOP_gauge_coeffs_t *coeffs,
		QDP_Subset subset, QDP_Subset osubset)
{
#define NC QDP_get_nc_M(staple)
  QDP_Lattice *lat = QDP_get_lattice_M(staple);
  int nd = QDP_ndim_L(lat);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);
  QLA_Real plaq = coeffs->plaquette;
  QLA_Real adpl = coeffs->adjoint_plaquette;
#if 1
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4, *temp5, *temp6;

  //temp1 = QDP_create_M();
  temp2 = QDP_create_M_L(lat);
  //temp3 = QDP_create_M();
  temp4 = QDP_create_M_L(lat);
  //temp5 = QDP_create_M();
  temp6 = QDP_create_M_L(lat);

  /* staple += u[nu](x) u[mu](x+nu) u*[nu](x+mu)
   *         + u*[nu](x-nu) u[mu](x-nu) u[nu](x-nu+mu) */
  for(int nu=0; nu<nd; nu++) {
    if (nu == mu) continue;

    temp1 = QDP_create_M_L(lat);
    temp3 = QDP_create_M_L(lat);
    temp5 = QDP_create_M_L(lat);

    QDP_M_eq_sM(temp1, u[nu], neighbor[mu], QDP_forward, QDP_all_L(lat));
    QDP_M_eq_Ma_times_M(temp2, u[nu], u[mu], osubset);
    QDP_M_eq_sM(temp3, u[mu], neighbor[nu], QDP_forward, subset);
    QDP_M_eq_M_times_M(temp4, temp2, temp1, osubset);
    QDP_M_eq_sM(temp5, temp4, neighbor[nu], QDP_backward, subset);
    QDP_M_eq_M_times_M(temp6, u[nu], temp3, subset);
    //QDP_M_peq_M_times_Ma(staple, temp6, temp1, subset);
    //QDP_M_peq_M(staple, temp5, subset);

    if(adpl==0) {
      QDP_M_peq_M_times_Ma(temp5, temp6, temp1, subset);
      QDP_M_peq_r_times_M(staple, &plaq, temp5, subset);
    } else {
      // FIXME
    }

    //QDP_discard_M(temp1);
    //QDP_discard_M(temp3);
    //QDP_discard_M(temp5);
    QDP_destroy_M(temp1);
    QDP_destroy_M(temp3);
    QDP_destroy_M(temp5);

  }  /* closes nu loop */

  //QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  //QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);
  //QDP_destroy_M(temp5);
  QDP_destroy_M(temp6);
#else
  QDP_ColorMatrix *t = QDP_create_M_L(lat);
  int nu, path[3];
  QDP_Subset subs[2];
  subs[0] = subset;
  subs[1] = osubset;
  for(nu=0; nu<nd; nu++) {
    if (nu == mu) continue;
    path[0] = 1+nu;
    path[1] = -(1+mu);
    path[2] = -(1+nu);
    path_prod(u, t, path, 3, 1, subs, neighsubeo);
    QDP_M_peq_M(staple, t, subset);
    path[0] = -(1+nu);
    path[1] = -(1+mu);
    path[2] = 1+nu;
    path_prod(u, t, path, 3, 1, subs, neighsubeo);
    QDP_M_peq_M(staple, t, subset);
  }
  QDP_destroy_M(t);
#endif
}

static void
get_staple_imp(QDP_ColorMatrix *staple, int mu, QDP_ColorMatrix **u,
	       int subl, QDP_Subset subs[], int (*neighsub)(int subl, int dir))
{
#define NC QDP_get_nc_M(staple)
  QDP_ColorMatrix *t;
  int path[5], bsubl;
  QDP_Lattice *lat = QDP_get_lattice_M(staple);
  int nd = QDP_ndim_L(lat);
  QLA_Real plaq = coeffs->plaquette;
  QLA_Real rect = coeffs->rectangle;
  QLA_Real pgm  = coeffs->parallelogram;
  QLA_Real adpl = coeffs->adjoint_plaquette;

  t = QDP_create_M_L(lat);
  for(mup=0; mup<2; mup++) {
    tm[mup] = QDP_create_M_L(lat);
    for(int nu=0; nu<8; nu++) {
      sm[mup][nu] = QDP_create_M_L(lat);
    }
  }

  //printf("test0\n");
  mup = 1 + mu;
  bsubl = neighsub(subl, mup);
  for(int nu=-nd; nu<=nd; nu++) {
    if ( nu==-mup || nu==0 || nu==mup ) continue;
    //s = QDP_create_M();
    path[0] = nu;
    path[1] = -mup;
    path[2] = -nu;
    path_prod(u, t, path, 3, bsubl, subs, neighsub);
    QDP_M_peq_M(staple, t, subs[subl]);
    //QDP_destroy_M(s);
  }
  //QDP_destroy_M(s);
  //printf("test1\n");
  //s = QDP_create_M();
  if(rect) {
    for(int nu=-nd; nu<=nd; nu++) {
      if ( nu==-mup || nu==0 || nu==mup ) continue;
      //s = QDP_create_M();
      path[0] = nu;
      path[1] = nu;
      path[2] = -mup;
      path[3] = -nu;
      path[4] = -nu;
      path_prod(u, t, path, 5, bsubl, subs, neighsub);
      QDP_M_peq_r_times_M(staple, &rect, t, subs[subl]);
      //QDP_destroy_M(s);
      //s = QDP_create_M();
      path[0] = nu;
      path[1] = -mup;
      path[2] = -mup;
      path[3] = -nu;
      path[4] = mup;
      path_prod(u, t, path, 5, bsubl, subs, neighsub);
      QDP_M_peq_r_times_M(staple, &rect, t, subs[subl]);
      //QDP_destroy_M(s);
      //s = QDP_create_M();
      path[0] = mup;
      path[1] = nu;
      path[2] = -mup;
      path[3] = -mup;
      path[4] = -nu;
      path_prod(u, t, path, 5, bsubl, subs, neighsub);
      QDP_M_peq_r_times_M(staple, &rect, t, subs[subl]);
      //QDP_destroy_M(s);
    }
  }
  if(pgm) {
    for(nu=-nd; nu<=nd; nu++) {
      if ( nu==-mup || nu==0 || nu==mup ) continue;
      for(int rho=-nd; rho<=nd; rho++) {
	if ( rho==-mup || rho==0 || rho==mup || rho==-nu || rho==nu ) continue;
	path[0] = nu;
	path[1] = rho;
	path[2] = -mup;
	path[3] = -nu;
	path[4] = -rho;
	path_prod(u, t, path, 5, bsubl, subs, neighsub);
	QDP_M_peq_r_times_M(staple, &para, t, subs[subl]);
      }
    }
  }
  QDP_destroy_M(t);

  for(int mup=0; mup<2; mup++) {
    for(int nu=0; nu<8; nu++) {
      QDP_destroy_M(sm[mup][nu]);
    }
    QDP_destroy_M(tm[mup]);
  }
}

void 
QOP_symanzik_1loop_gauge_staple_qdp(QOP_info_t *info,
				    QDP_ColorMatrix *staple,
				    int mu,
				    QDP_ColorMatrix *links[],
				    QOP_gauge_coeffs_t *coeffs,
				    QDP_Subset subs[], int subi)
{
  int imp = (coeffs->rectangle!=0)||(coeffs->parallelogram!=0);
  if(imp) {  // only works in 4d
#if 0
    QDP_Subset subs[2];
    subs[0] = QDP_all_L(lat);
    subs[1] = QDP_all_L(lat);
    get_staple_imp(staple, mu, u, 1, subs, neighsubeo);
#else
    get_staple_imp(staple, mu, links, subi, subs, neighsub32);
#endif
  } else { // should work in any dimension
    get_staple_plaq(staple, mu, links, coeffs, subs[subi], subs[1-subi]);
  }
}
