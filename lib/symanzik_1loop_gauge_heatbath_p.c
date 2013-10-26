#include <qop_internal.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//cbs = QOP_get_sub32(lat);
///get_staple(staple, mu, links, cbs, cb);
//su2_extract(r, m, i, j);
//su2_fill(&t, r, i, j);

static QLA_RandomState *rs;
static QLA_Real fac;

static void
su2_extract(NCPROT QLA_Real r[4], QLA_ColorMatrix(*m), int i, int j)
{
  QLA_Complex *a00, *a01, *a10, *a11;
  a00 = &QLA_elem_M(*m, i, i);
  a01 = &QLA_elem_M(*m, i, j);
  a10 = &QLA_elem_M(*m, j, i);
  a11 = &QLA_elem_M(*m, j, j);
  r[0] = QLA_real(*a00) + QLA_real(*a11);
  r[1] = QLA_imag(*a01) + QLA_imag(*a10);
  r[2] = QLA_real(*a01) - QLA_real(*a10);
  r[3] = QLA_imag(*a00) - QLA_imag(*a11);
}

static void
su2_fill(NCPROT QLA_ColorMatrix(*m), QLA_Real r[4], int i, int j)
{
  QLA_Complex z;
  QLA_c_eq_r(z, 1);
  QLA_M_eq_c(m, &z);

  QLA_c_eq_r_plus_i_r(z, r[0], r[3]);
  QLA_M_eq_elem_C(m, &z, i, i);

  QLA_c_eq_r_plus_i_r(z, r[2], r[1]);
  QLA_M_eq_elem_C(m, &z, i, j);

  r[2] = -r[2];
  QLA_c_eq_r_plus_i_r(z, r[2], r[1]);
  QLA_M_eq_elem_C(m, &z, j, i);

  r[3] = -r[3];
  QLA_c_eq_r_plus_i_r(z, r[0], r[3]);
  QLA_M_eq_elem_C(m, &z, j, j);
}

static void
get_hb2(QLA_Real *b, QLA_Real al, QLA_RandomState *srs)
{
  QLA_Real d, r, r2, rho, xr;

  if(al <= 2.0) {  /* creutz algorithm */
    QLA_Real xl, xd, a0;
    xl = exp(-2.0*al);
    xd = 1.0 - xl;
    for(int k=0; k<20; k++) {
      QLA_Real xr1, xr2, s;
      QLA_R_eq_random_S(&xr1, srs);
      QLA_R_eq_random_S(&xr2, srs);
      s = xl + xd*xr1;
      a0 = 1.0 + log(s)/al;
      if((1.0-a0*a0) > xr2*xr2) break;
    }
    d = 1.0 - a0;
  } else {  /* k-p algorithm */
    for(int k=0; k<20; k++) {
      QLA_Real xr1, xr2, xr3, xr4;
      QLA_R_eq_random_S(&xr1, srs);
      xr1 = log(xr1 + 1.e-10);
      QLA_R_eq_random_S(&xr2, srs);
      xr2 = log(xr2 + 1.e-10);
      QLA_R_eq_random_S(&xr3, srs);
      QLA_R_eq_random_S(&xr4, srs);
      xr3 = cos(2.0*M_PI*xr3);
      d = -(xr2 + xr1*xr3*xr3)/al;
      if((1.0 - 0.5*d) > xr4*xr4) break;
    }
  }

  b[0] = 1.0 - d;
  r2 = fabs(1.0 - b[0]*b[0]);
  r = sqrt(r2);

  QLA_R_eq_random_S(&xr, srs);
  b[3] = (2.0*xr - 1.0)*r;

  rho = sqrt(fabs(r2 - b[3]*b[3]));

  QLA_R_eq_random_S(&xr, srs);
  xr *= 2.0*M_PI;
  b[1] = rho*cos(xr);
  b[2] = rho*sin(xr);
}

static void
hb_func(NCPROT QLA_ColorMatrix(*m), int site)
{
  QLA_RandomState *srs = rs + site;
  if(QDP_Nc==1) {
    // FIXME
    //QLA_R_eq_random_S(&xr, srs);
    QLA_Complex z, zs, t;
    QLA_C_eq_elem_M(&z, m, 0, 0);
    QLA_c_eq_ca(zs, z);
    QLA_C_eq_C_divide_C(&t, &zs, &z);
    QLA_M_eq_elem_C(m, &t, 0, 0);
  } else {
    QLA_ColorMatrix(s);
    QLA_ColorMatrix(t);
    QLA_ColorMatrix(tt);
    QLA_Complex one;
    QLA_c_eq_r(one, 1);
    QLA_M_eq_c(&s, &one);

    /* Loop over SU(2) subgroup index */
    for(int i=0; i<QLA_Nc; i++) {
      for(int j=i+1; j<QLA_Nc; j++) {
	QLA_Real a[4], b[4], r[4], rn, rl;

	su2_extract(NCARG r, m, i, j);
	rn = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3] );
	rl = fac*rn;
	if(rn<1e-10) {
	  a[0] = 1; a[1] = a[2] = a[3] = 0;
	} else {
	  rn = 1/rn;
	  a[0] =  rn*r[0];
	  a[1] = -rn*r[1];
	  a[2] = -rn*r[2];
	  a[3] = -rn*r[3];
	}

	get_hb2(b, rl, srs);
	//b[0] = 1; b[1] = b[2] = b[3] = 0;

	r[0] = b[0]*a[0] - b[1]*a[1] - b[2]*a[2] - b[3]*a[3];
	r[1] = b[0]*a[1] + b[1]*a[0] - b[2]*a[3] + b[3]*a[2];
	r[2] = b[0]*a[2] + b[2]*a[0] - b[3]*a[1] + b[1]*a[3];
	r[3] = b[0]*a[3] + b[3]*a[0] - b[1]*a[2] + b[2]*a[1];

	su2_fill(NCARG &t, r, i, j);
	QLA_M_eq_M_times_M(&tt, &t, &s);
	QLA_M_eq_M(&s, &tt);
	QLA_M_eq_M_times_M(&tt, &t, m);
	QLA_M_eq_M(m, &tt);
      }
    }
    QLA_M_eq_M(m, &s);
  }
}

void 
QOP_symanzik_1loop_gauge_heatbath_qdp(QOP_info_t *info,
				      QDP_ColorMatrix *links[],
				      QLA_Real beta,
				      QOP_gauge_coeffs_t *coeffs,
				      QDP_RandomState *rs0)
{
#define NC QDP_get_nc(links[0])
  double dtime = QOP_time();
  double nflops = 0;
  fac = beta/QLA_Nc;
  QLA_Real plaq = fac*coeffs->plaquette;
  QLA_Real rect = fac*coeffs->rectangle;
  QLA_Real pgm  = fac*coeffs->parallelogram;
  QLA_Real adpl = fac*fac*coeffs->adjoint_plaquette;
  coeffs->adjoint_plaquette *= fac;
  int imp = (coeffs->rectangle!=0)||(coeffs->parallelogram!=0);
  QDP_Lattice *lat = QDP_get_lattice_M(links[0]);
  int nd = QDP_ndim_L(lat);
  QDP_Subset *cbs=QDP_even_and_odd_L(lat);
  int ncb = 2;
  if(imp) {
    ncb = 32;
    cbs = QOP_get_sub32(lat);
  }

  QDP_ColorMatrix *staple = QDP_create_M_L(lat);
  QDP_ColorMatrix *v = QDP_create_M_L(lat);
  QDP_ColorMatrix *tmp = QDP_create_M_L(lat);
  rs = QDP_expose_S(rs0);

  for(int cb=0; cb<ncb; cb++) {
    QDP_Subset subset = cbs[cb];
    for(int mu=0; mu<nd; mu++) {
      QDP_M_eq_zero(staple, subset);
      QOP_symanzik_1loop_gauge_staple_qdp(info, staple, mu, links, coeffs, cbs, cb);
      //heatbath(links[mu], staple, cbs[cb]);
      QDP_M_eq_M_times_Ma(v, links[mu], staple, subset);
      QDP_M_eq_funci(v, hb_func, subset);
      QDP_M_eq_M_times_M(tmp, v, links[mu], subset);
      QDP_M_eq_M(links[mu], tmp, subset);
    }
  }

  QDP_reset_S(rs0);
  QDP_destroy_M(tmp);
  QDP_destroy_M(v);
  QDP_destroy_M(staple);
  coeffs->adjoint_plaquette /= fac;

  info->final_sec = QOP_time() - dtime;
  info->final_flop = nflops*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
#undef NC
}
