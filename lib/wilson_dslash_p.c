/* adapted from MILC version 6 */

/**********************
** original comments **
**********************/
/******* d_congrad2.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 6 */
/* Wilson fermions */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
   the fermion spinors live on even sites only.  In other words, if
   Dslash_oe is the dslash operator with its source on even sites and
   its result on odd sites, etc.:

   without LU:
   M = 1 - kappa*( Dslash_eo + DSLASH_oe )
   with LU:
   M = 1 - kappa^2 * Dslash_eo * Dslash_oe
*/
/**************************
** end original comments **
**************************/

#include <string.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU
//#define CLOV_FUNC

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_wilson_inited;
extern int QOP_wilson_style;
extern int QOP_wilson_nsvec;
extern int QOP_wilson_nvec;
extern int QOP_wilson_cgtype;
extern int QOP_wilson_optnum;

static int old_style=-1;
//static int old_nsvec=-1;
//static int old_nvec=-1;
static int old_optnum=-1;

#define NTMPSUB 2
#define NTMP (3*NTMPSUB)
#define NHTMP 20
#define NDTMP 12
static int dslash_setup = 0;
static QDP_HalfFermion *htemp[NTMP][NHTMP];
static QDP_DiracFermion *dtemp[NTMP][NDTMP];
static QDP_DiracFermion *tin[NTMP];
#define tmpnum(eo,n) ((eo)+3*((n)-1))
#define tmpsub(eo,n) tin[tmpnum(eo,n)]

#define check_setup(flw) \
{ \
  if( (!dslash_setup) || (QOP_wilson_optnum != old_optnum) ) { \
    reset_temps(); \
  } \
  if( flw->dblstored != dblstore_style(QOP_wilson_style) ) { \
    double_store(flw); \
  } \
}

static void
free_temps(void)
{
  if(dslash_setup) {
    int i, j;

    for(i=0; i<NTMP; i++) {
      QDP_destroy_D(tin[i]);
    }

    if(shiftd_style(old_style)) {
      for(i=0; i<NTMP; i++) {
	for(j=0; j<NDTMP; j++) {
	  QDP_destroy_D(dtemp[i][j]);
	}
      }
    } else {
      for(i=0; i<NTMP; i++) {
	for(j=0; j<NHTMP; j++) {
	  QDP_destroy_H(htemp[i][j]);
	}
      }
    }
  }
  dslash_setup = 0;
}

static void
reset_temps(void)
{
  int i, j;

  free_temps();

  for(i=0; i<NTMP; i++) {
    tin[i] = QDP_create_D();
  }

  if(shiftd_style(QOP_wilson_style)) {
    for(i=0; i<NTMP; i++) {
      for(j=0; j<NDTMP; j++) {
	dtemp[i][j] = QDP_create_D();
      }
    }
  } else {
    for(i=0; i<NTMP; i++) {
      for(j=0; j<NHTMP; j++) {
	htemp[i][j] = QDP_create_H();
      }
    }
  }
  dslash_setup = 1;
  old_style = QOP_wilson_style;
  old_optnum = QOP_wilson_optnum;
}

static void
double_store(QOP_FermionLinksWilson *flw)
{
  int i;

  if(flw->dblstored) {
    for(i=0; i<4; i++) {
      QDP_destroy_M(flw->bcklinks[i]);
    }
    flw->dblstored = 0;
  }

  if(dblstore_style(QOP_wilson_style)) {
    for(i=0; i<4; i++) {
      flw->bcklinks[i] = QDP_create_M();
    }
    for(i=0; i<4; i++) {
      flw->dbllinks[2*i] = flw->links[i];
      flw->dbllinks[2*i+1] = flw->bcklinks[i];
    }
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<4; i++) {
      QDP_M_eq_sM(m, flw->links[i], QDP_neighbor[i], QDP_backward, QDP_all);
      QDP_M_eq_Ma(flw->bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    flw->dblstored = dblstore_style(QOP_wilson_style);
  }
}

QDP_DiracFermion *
QOPPC(wilson_dslash_get_tmp)(QOP_FermionLinksWilson *flw,
			     QOP_evenodd_t eo, int n)
{
  check_setup(flw);
  if(n>=1 && n<=NTMPSUB) return tmpsub(eo,n);
  else return NULL;
}

/* link routines */

static void
get_clov(QDP_DiracPropagator *clov, QDP_ColorMatrix *links[], double csw)
{
  QLA_DiracPropagator p;
  int i, j;
  QLA_P_eq_zero(&p);
  for(i=0; i<QLA_Nc; i++) for(j=0; j<QLA_Ns; j++)
    QLA_c_eq_r(QLA_elem_P(p, i, j, i, j), 1);
  QDP_P_eq_p(clov, &p, QDP_all);
}

static void
get_clovinv(REAL *clovinv, REAL *clov)
{
  memcpy(clovinv, clov, QDP_sites_on_node*2*6*6*sizeof(REAL));
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
{
  QOP_FermionLinksWilson *flw;
  QOP_GaugeField *gf;

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_wilson_convert_L_from_qdp(gf->links, NULL);

  if(clov!=NULL) {
    QOP_malloc(flw->clov, REAL, QDP_sites_on_node*2*6*6);
    QOP_malloc(flw->clovinv, REAL, QDP_sites_on_node*2*6*6);
    memcpy(flw->clov, clov, QDP_sites_on_node*2*6*6*sizeof(REAL));
    get_clovinv(flw->clovinv, flw->clov);
  } else {
    flw->clov = NULL;
    flw->clovinv = NULL;
  }
  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qdpclov = NULL;
  flw->qopgf = gf;

  return flw;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			   QOP_GaugeField *gauge)
{
  QOP_FermionLinksWilson *flw;
  QDP_DiracPropagator *clov;
  if(coeffs->clov_c==0.) {
    clov = NULL;
  } else {
    clov = QDP_create_P();
    get_clov(clov, gauge->links, coeffs->clov_c);
  }
  if(coeffs->aniso!=0.) {
    int i;
    for(i=0; i<3; i++) {
      QLA_Real f = coeffs->aniso;
      QDP_M_eq_r_times_M(gauge->links[i], &f, gauge->links[i], QDP_all);
    }
  }
  flw = QOP_wilson_create_L_from_qdp(gauge->links, clov);
  if(clov) QDP_destroy_P(clov);
  return flw;
}

void
QOP_wilson_extract_L_to_raw(REAL *links[], REAL *clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_extract_L_to_raw unimplemented.");
}

void
QOP_wilson_destroy_L(QOP_FermionLinksWilson *flw)
{
  int i;

  if(flw->qopgf) {
    QOP_destroy_G(flw->qopgf);
  } else {
    for(i=0; i<4; i++) QDP_destroy_M(flw->links[i]);
    free(flw->links);
  }
  if(flw->dblstored) {
    for(i=0; i<4; i++) QDP_destroy_M(flw->bcklinks[i]);
  }
  if(flw->qdpclov) {
    QDP_destroy_P(flw->qdpclov);
  }
  if(flw->clov) {
    free(flw->clov);
    free(flw->clovinv);
  }
  free(flw->bcklinks);
  free(flw->dbllinks);
  free(flw);
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], REAL *clov,
			      QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_convert_L_from_raw unimplemented");
  return NULL;
}

void
QOP_wilson_convert_L_to_raw(REAL ***links, REAL **clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  QOP_error("QOP_wilson_convert_L_to_raw unimplemented");
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			    QOP_GaugeField *gauge)
{
  QOP_error("QOP_wilson_convert_L_from_G unimplemented");
  return NULL;
}

QOP_GaugeField *
QOP_wilson_convert_L_to_G(QOP_FermionLinksWilson *links)
{
  QOP_error("QOP_wilson_convert_L_to_G unimplemented");
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_qdp(QDP_ColorMatrix *links[],
			     QDP_DiracPropagator *clov)
{
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_wilson_convert_L_from_qdp(newlinks, clov);
  flw->qdpclov = NULL;

  return flw;
}

void
QOP_wilson_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QDP_DiracPropagator *clov,
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_extract_L_to_qdp unimplemented");
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_qdp(QDP_ColorMatrix *links[],
			      QDP_DiracPropagator *clov)
{
  QOP_FermionLinksWilson *flw;
  int i;

  QOP_malloc(flw, QOPPC(FermionLinksWilson), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);
  if(clov!=NULL) {
    int size = QDP_sites_on_node*2*6*6;
    QOP_malloc(flw->clov, REAL, size);
    QOP_malloc(flw->clovinv, REAL, size);
    {
      QLA_DiracPropagator *dp;
      int x, b, i, ic, j, jc, is, js, k=0;
      dp = QDP_expose_P(clov);
      for(x=0; x<QDP_sites_on_node; x++) {
	for(b=0; b<2; b++) { // two chiral blocks
	  // first the diagonal
	  for(i=0; i<6; i++) {
	    ic = i/2;
	    is = 2*b + i%2;
	    flw->clov[k++] = QLA_real(QLA_elem_P(dp[x], ic, is, ic, is));
	  }
	  // now the offdiagonal
	  for(i=0; i<6; i++) {
	    ic = i/2;
	    is = 2*b + i%2;
	    for(j=i+1; j<6; j++) {
	      QLA_Complex z1, z2;
	      jc = j/2;
	      js = 2*b + j%2;
	      //QLA_c_eq_c_plus_ca(z1, QLA_elem_P(dp[x], ic, is, jc, js),
	      //                   QLA_elem_P(dp[x], jc, js, ic, is));
	      QLA_c_eq_c(z1, QLA_elem_P(dp[x], ic, is, jc, js));
	      QLA_c_peq_ca(z1, QLA_elem_P(dp[x], jc, js, ic, is));
	      QLA_c_eq_r_times_c(z2, 0.5, z1);
	      flw->clov[k++] = QLA_real(z2);
	      flw->clov[k++] = QLA_imag(z2);
	    }
	  }
	}
      }
      QDP_reset_P(clov);
    }
    get_clovinv(flw->clovinv, flw->clov);
  } else {
    flw->clov = NULL;
    flw->clovinv = NULL;
  }

  flw->dblstored = 0;
  for(i=0; i<4; i++) {
    flw->links[i] = links[i];
  }
  // scale links
  for(i=0; i<4; i++) {
    QLA_Real f = -0.5;
    QDP_M_eq_r_times_M(flw->links[i], &f, flw->links[i], QDP_all);
  }

  check_setup(flw);

  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qopgf = NULL;
  flw->qdpclov = clov;
  return flw;
}

void
QOP_wilson_convert_L_to_qdp(QDP_ColorMatrix ***links,
			    QDP_DiracPropagator **clov,
			    QOP_FermionLinksWilson *src)
{
  QOP_error("QOP_wilson_convert_L_to_qdp unimplemented");
}

/* dslash */

static void
clov(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out, QDP_DiracFermion *in,
     REAL kappa, QDP_DiracFermion *dsl, QDP_Subset subset);

static void
wilson_dslash0(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n);

static void
wilson_dslash1(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n);

#define wilson_hop(flw, dest, src, sign, eo) \
{ \
  QDP_DiracFermion *tsrc = src; \
  int _n = 1; \
  while(1) { \
    if(src==tmpsub(eo,_n)) break; \
    if(_n==NTMPSUB) { \
      _n = 1; \
      tsrc = tmpsub(eo,_n); \
      QDP_D_eq_D(tsrc, src, qdpsub(oppsub(eo))); \
      break; \
    } \
    _n++; \
  } \
  /*printf("%i %i\n", eo, _n);*/ \
  if(dblstore_style(QOP_wilson_style)) { \
    wilson_dslash1(flw, dest, tsrc, sign, eo, _n); \
  } else { \
    wilson_dslash0(flw, dest, tsrc, sign, eo, _n); \
  } \
}

void
QOP_wilson_dslash(QOP_info_t *info,
		  QOP_FermionLinksWilson *flw,
		  REAL kappa,
		  int sign,
		  QOP_DiracFermion *out,
		  QOP_DiracFermion *in,
		  QOP_evenodd_t eo_out,
		  QOP_evenodd_t eo_in)
{
  QOP_wilson_dslash_qdp(info,flw,kappa,sign,out->df,in->df,eo_out,eo_in);
}

void
QOP_wilson_dslash_qdp(QOP_info_t *info,
		      QOP_FermionLinksWilson *flw,
                      REAL kappa,
		      int sign,
		      QDP_DiracFermion *out,
		      QDP_DiracFermion *in,
		      QOP_evenodd_t eo_out,
		      QOP_evenodd_t eo_in)
{
  //printf("testd1\n");
  check_setup(flw);
  //printf("testd2\n");

  if(eo_in==eo_out) {
    if(eo_out==QOP_EVENODD) {
      //printf("testd21\n");
      wilson_hop(flw, out, in, sign, QOP_EVENODD);
      //printf("testd22\n");
      clov(flw, out, in, kappa, out, QDP_all);
      //printf("testd23\n");
    } else if(eo_out==QOP_EVEN) {
      clov(flw, out, in, kappa, NULL, QDP_even);
    } else {
      clov(flw, out, in, kappa, NULL, QDP_odd);
    }
  } else {
    if(eo_out==QOP_EVEN || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_ODD) {
	wilson_hop(flw, out, in, sign, QOP_EVEN);
      } else if(eo_in==QOP_EVEN) {
	clov(flw, out, in, kappa, NULL, QDP_even);
      } else {
	wilson_hop(flw, out, in, sign, QOP_EVEN);
	clov(flw, out, in, kappa, out, QDP_even);
      }
    }
    if(eo_out==QOP_ODD || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_EVEN) {
	wilson_hop(flw, out, in, sign, QOP_ODD);
      } else if(eo_in==QOP_ODD) {
	clov(flw, out, in, kappa, NULL, QDP_odd);
      } else {
	wilson_hop(flw, out, in, sign, QOP_ODD);
	clov(flw, out, in, kappa, out, QDP_odd);
      }
    }
  }
  //printf("testd3\n");
}

void
QOP_wilson_diaginv(QOP_info_t *info,
		   QOP_FermionLinksWilson *flw,
		   REAL kappa,
		   QOP_DiracFermion *out,
		   QOP_DiracFermion *in,
		   QOP_evenodd_t eo)
{
  QOP_wilson_diaginv_qdp(info,flw,kappa,out->df,in->df,eo);
}

void
QOP_wilson_diaginv_qdp(QOP_info_t *info,
		       QOP_FermionLinksWilson *flw,
		       REAL kappa,
		       QDP_DiracFermion *out,
		       QDP_DiracFermion *in,
		       QOP_evenodd_t eo)
{
  QLA_Real f = 2*kappa;
  QDP_D_eq_r_times_D(out, &f, in, qdpsub(eo));
}



QLA_Real clov_mkappa, *clov_clov;
QLA_DiracFermion *clov_dsl, *clov_in;

#define cmplx(x) (*((QLA_Complex *)(&(x))))

#ifdef CLOV_FUNC
static void
clov_func(QLA_DiracFermion *out, int coords[])
{
  int b, i, j, ic, jc, is, js, k, x, xb;
  QLA_DiracFermion *dsl;

  x = QDP_index(coords);  // site offset
  if(clov_dsl!=NULL) dsl = &clov_dsl[x]; else dsl = out;

  for(b=0; b<2; b++) {
    xb = 36*(2*x+b);  // chiral block offset (in REALs)
    for(i=0; i<6; i++) {
      QLA_Complex z;
      ic = i/2;
      is = 2*b + i%2;

      QLA_c_eq_r(z, 0.);

      //#if 0
      // lower triangular part comes from adjoint of upper
      k = xb + 6 + 2*(i-1); // block + skip diag + cmplx elem
      for(j=0; j<i; j++) {
	jc = j/2;
	js = 2*b + j%2;
	QLA_c_peq_ca_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	k += 2*(4 - j);
      }
      // diagonal part
      QLA_c_peq_r_times_c(z, clov_clov[xb], QLA_elem_D(clov_in[x],ic,is));
      // upper triangular part
      for(j=i+1; j<6; j++) {
	jc = j/2;
	js = 2*b + j%2;
	k += 2;
	QLA_c_peq_c_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
      }
      //#endif

      QLA_c_peq_r_times_c(z, clov_mkappa, QLA_elem_D(*dsl, ic, is));
      QLA_c_eq_c_plus_c(QLA_elem_D(*out, ic, is), z, QLA_elem_D(clov_in[x],ic,is));
    }
  }
}
#endif

static void
apply_clov(REAL *clov, QDP_DiracFermion *out, QLA_Real mkappa,
	   QDP_DiracFermion *dsl, QDP_DiracFermion *in, QDP_Subset subset)
{
#ifdef CLOV_FUNC
  clov_in = QDP_expose_D(in);
  clov_mkappa = mkappa;
  clov_clov = clov;
  if(dsl==out) clov_dsl = NULL; else clov_dsl = QDP_expose_D(dsl);
  QDP_D_eq_func(out, clov_func, subset);
  if(dsl!=out) QDP_reset_D(dsl);
  QDP_reset_D(in);
#else
  QLA_DiracFermion *clov_out;
  clov_out = QDP_expose_D(out);
  clov_in = QDP_expose_D(in);
  clov_mkappa = mkappa;
  clov_clov = clov;
  if(dsl==out) clov_dsl = clov_out;
  else clov_dsl = QDP_expose_D(dsl);
  {
    int x, start, end;
    if(subset==QDP_odd) start = QDP_subset_len(QDP_even);
    else start = 0;
    end = start + QDP_subset_len(subset);
    for(x=start; x<end; x++) {
      int b;
      QLA_DiracFermion *dsl;
      dsl = &clov_dsl[x];
      for(b=0; b<2; b++) {
	int xb;
	xb = 36*(2*x+b);  // chiral block offset (in REALs)
#if 0 // loop version
	int i, j, ic, jc, is, js, k;
	for(i=0; i<6; i++) {
	  QLA_Complex z;
	  ic = i/2;
	  is = 2*b + i%2;

	  QLA_c_eq_r(z, 0.);
	  //#if 0
	  // lower triangular part comes from adjoint of upper
	  k = xb + 6 + 2*(i-1); // block + skip diag + cmplx elem
	  for(j=0; j<i; j++) {
	    jc = j/2;
	    js = 2*b + j%2;
	    QLA_c_peq_ca_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	    k += 2*(4 - j);
	  }
	  // diagonal part
	  QLA_c_peq_r_times_c(z, clov_clov[xb], QLA_elem_D(clov_in[x],ic,is));
	  // upper triangular part
	  for(j=i+1; j<6; j++) {
	    jc = j/2;
	    js = 2*b + j%2;
	    k += 2;
	    QLA_c_peq_c_times_c(z, cmplx(clov_clov[k]), QLA_elem_D(clov_in[x],jc,js));
	  }
	  //#endif

	  QLA_c_peq_r_times_c(z, clov_mkappa, QLA_elem_D(*dsl, ic, is));
	  QLA_c_eq_c_plus_c(QLA_elem_D(clov_out[x], ic, is), z, QLA_elem_D(clov_in[x],ic,is));
	}
#else  // unrolled version
#define clov_diag(i) clov_clov[xb+i]
#define clov_offd(i) cmplx(clov_clov[xb+6+i])
#define src(i) QLA_elem_D(clov_in[x],i/2,2*b+i%2)
#define dsrc(i) QLA_elem_D(*dsl,i/2,2*b+i%2)
#define dest(i) QLA_elem_D(clov_out[x],i/2,2*b+i%2)
	{
	  QLA_Complex z;

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_r_times_c(z, clov_diag(0), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(0), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(1), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(2), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(3), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(4), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(0));
	  QLA_c_eq_c_plus_c(dest(0), z, src(0));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(0), src(0));
	  QLA_c_peq_r_times_c(z, clov_diag(1), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(5), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(6), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(7), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(8), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(1));
	  QLA_c_eq_c_plus_c(dest(1), z, src(1));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(1), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(5), src(1));
	  QLA_c_peq_r_times_c(z, clov_diag(2), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(9), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(10), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(11), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(2));
	  QLA_c_eq_c_plus_c(dest(2), z, src(2));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(2), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(6), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(9), src(2));
	  QLA_c_peq_r_times_c(z, clov_diag(3), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(12), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(13), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(3));
	  QLA_c_eq_c_plus_c(dest(3), z, src(3));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(3), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(7), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(10), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(12), src(3));
	  QLA_c_peq_r_times_c(z, clov_diag(4), src(4));
	  QLA_c_peq_c_times_c(z, clov_offd(14), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(4));
	  QLA_c_eq_c_plus_c(dest(4), z, src(4));

	  QLA_c_eq_r(z, 0.);
	  QLA_c_peq_c_times_c(z, clov_offd(4), src(0));
	  QLA_c_peq_c_times_c(z, clov_offd(8), src(1));
	  QLA_c_peq_c_times_c(z, clov_offd(11), src(2));
	  QLA_c_peq_c_times_c(z, clov_offd(13), src(3));
	  QLA_c_peq_c_times_c(z, clov_offd(14), src(4));
	  QLA_c_peq_r_times_c(z, clov_diag(5), src(5));
	  QLA_c_peq_r_times_c(z, clov_mkappa, dsrc(5));
	  QLA_c_eq_c_plus_c(dest(5), z, src(5));
	}
#endif
      }
    }
  }
  QDP_reset_D(in);
  QDP_reset_D(out);
  if(dsl!=out) QDP_reset_D(dsl);
#endif
}

static void
clov(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out, QDP_DiracFermion *in,
     REAL kappa, QDP_DiracFermion *dsl, QDP_Subset subset)
{
  QLA_Real m4 = 0.5/kappa;
  if(flw->clov==NULL) {
    if(dsl==NULL) {
      QDP_D_eq_r_times_D(out, &m4, in, subset);
    } else {
      QDP_D_eq_r_times_D_plus_D(out, &m4, in, dsl, subset);
    }
  } else {
    apply_clov(flw->clov, out, m4, dsl, in, subset);
  }
}

static void
clovinv(QOP_FermionLinksWilson *flw, QDP_DiracFermion *out,
	QDP_DiracFermion *in, REAL kappa, QDP_DiracFermion *dsl,
	QDP_Subset subset)
{
  QLA_Real m4 = 0.5/kappa;
  if(flw->clov==NULL) {
    QDP_D_eq_r_times_D_plus_D(out, &m4, dsl, in, subset);
  } else {
    apply_clov(flw->clovinv, out, m4, dsl, in, subset);
  }
}


/************ dslash *************/

/* Special dslash for use by congrad.  Uses restart_gather() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather, otherwise use restart_gather.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
static void
wilson_dslash0(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n)
{
  int mu, ntmp;
  QDP_DiracFermion *vsrc[4];
  QDP_DiracFermion *vdest[4];
  QDP_ShiftDir fwd[4], bck[4];
  int dir[4], sgn[4], msgn[4];
  QDP_Subset subset, othersubset;
  subset = qdpsub(eo);
  othersubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);

  sign = -sign;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vdest[mu] = dest;
    fwd[mu] = QDP_forward;
    bck[mu] = QDP_backward;
    dir[mu] = mu;
    sgn[mu] = sign;
    msgn[mu] = -sign;
  }
  //sgn[1] = -sign;
  //msgn[1] = sign;

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  printf0("dslash0\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_D_veq_sD\n");
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, QDP_neighbor+mu, fwd+mu, subset,
		   QOP_wilson_nsvec);
      printf0("end QDP_D_veq_sD\n");
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_D\n");
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, QDP_neighbor+mu, fwd+mu,
		   subset, QOP_wilson_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_wilson_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
#if 0
      QDP_H_veq_spproj_Ma_times_D(hf+mu, fwdlinks+mu, vsrc+mu,
                               dir+mu, msgn+mu, othersubset, QOP_wilson_nsvec);
      QDP_D_veq_sprecon_H(dtemp[ntmp]+8+mu, hf+mu,
                               dir+mu, msgn+mu, othersubset, QOP_wilson_nsvec);
#endif
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_Ma_times_D\n");
      QDP_H_veq_spproj_Ma_times_D(htemp[ntmp]+12+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+12+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  printf0("dslash0 - fwd\n");
  QDP_D_eq_zero(dest, subset);

  if(shiftd_style(QOP_wilson_style)) {
    //QDP_HalfFermion *hf[4];
    //for(mu=0; mu<4; mu++) hf[mu] = QDP_create_H();
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->links+mu, dtemp[ntmp]+mu,
				  dir+mu, sgn+mu, subset, QOP_wilson_nvec);
#if 0
      QDP_H_veq_spproj_M_times_D(hf+mu, flw->links+mu, dtemp[ntmp]+mu,
				 dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      QDP_D_vpeq_sprecon_H(vdest+mu, hf+mu,
			   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
#endif
    }
    //for(mu=0; mu<4; mu++) QDP_destroy_H(hf[mu]);
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->links+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_D(vdest+mu, dtemp[ntmp]+4+mu, subset, QOP_wilson_nvec);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_H(vdest+mu, htemp[ntmp]+4+mu, dir+mu, msgn+mu, subset,
			   QOP_wilson_nvec);
    }
  }

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
} /* end of dslash_special_qdp() */

/* Special dslash for use by congrad.  Uses restart_gather() when
   possible. Last argument is an integer, which will tell if
   gathers have been started.  If is_started=0,use
   start_gather, otherwise use restart_gather.
   Argument "tag" is a vector of a msg_tag *'s to use for
   the gathers.
   The calling program must clean up the gathers! */
static void
wilson_dslash1(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n)
{
  int mu, ntmp;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];
  QDP_Subset subset, othersubset;
  subset = qdpsub(eo);
  othersubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);

  sign = -sign;

  for(mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vsrc[mu+4] = src;
    vdest[mu] = dest;
    vdest[mu+4] = dest;
    dir[2*mu] = mu;
    dir[2*mu+1] = mu;
    sgn[2*mu] = sign;
    sgn[2*mu+1] = -sign;
    sh[2*mu] = QDP_neighbor[mu];
    sh[2*mu+1] = QDP_neighbor[mu];
    sd[2*mu] = QDP_forward;
    sd[2*mu+1] = QDP_backward;
  }
  //sgn[2] = -sign;
  //sgn[3] = sign;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  //printf0("ds1 1\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, sh+mu, sd+mu, subset,
		   QOP_wilson_nsvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, sh+mu, sd+mu, subset,
		   QOP_wilson_nsvec);
    }
  }
  //printf0("ds1 2\n");

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add to dest */

  QDP_D_eq_zero(dest, subset);
  //printf0("ds1 3\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->dbllinks+mu, dtemp[ntmp]+mu,
			          dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->dbllinks+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
  }
  //printf0("ds1 4\n");

  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
} /* end of dslash_special_qdp() */




#if 0


static void
wilson_mdslash1(QOP_FermionLinksWilson *flw,
		QDP_DiracFermion *out, QDP_DiracFermion *in,
		int sign, QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
#ifdef LU

  //printf0("here3\n");
  dslash_special_qdp(flw, tt2, in, sign, QDP_odd, 2);
  dslash_special_qdp(flw, out, tt2, sign, QDP_even, 3);
  QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);

#else

  //printf0("here6\n");
  dslash_special_qdp(flw, out, in, sign, QDP_all, 1);
  //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
  clov(flw, out, &mkappa, out, in, QDP_all);

#endif
}

static void
wilson_mdslash2(QOP_FermionLinksWilson *flw,
		QDP_DiracFermion *out, QDP_DiracFermion *in,
		QDP_Subset subset, QDP_Subset othersubset,
		QLA_Real mkappa)
{
#ifdef LU

  if(QOP_wilson_cgtype==1) {
    if(in==flw->cgp) {
      dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 0);
      dslash_special_qdp(flw, out, tt1, 1, QDP_even, 1);
      QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);
    } else {
      dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 2);
      dslash_special_qdp(flw, out, tt1, 1, QDP_even, 1);
      QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_even);
    }
  } else {
    dslash_special_qdp(flw, tt1, in, 1, QDP_odd, 0);
    dslash_special_qdp(flw, ttt, tt1, 1, QDP_even, 1);
    QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_even);

    dslash_special_qdp(flw, tt2, ttt, -1, QDP_odd, 2);
    dslash_special_qdp(flw, out, tt2, -1, QDP_even, 3);
    QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_even);
  }

#else

  if(QOP_wilson_cgtype==1) {
    if(in==flw->cgp) {
      dslash_special_qdp(flw, out, in, 1, QDP_all, 0);
      //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
      clov(flw, out, &mkappa, out, in, QDP_all);
    } else {
      dslash_special_qdp(flw, out, in, 1, QDP_all, 1);
      //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, in, QDP_all);
      clov(flw, out, &mkappa, out, in, QDP_all);
    }
  } else {
    //printf0("here6\n");
    dslash_special_qdp(flw, ttt, in, 1, QDP_all, 0);
    //QDP_D_eq_r_times_D_plus_D(ttt, &mkappa, ttt, in, QDP_all);
    clov(flw, ttt, &mkappa, ttt, in, QDP_all);
    //printf0("here7\n");
    dslash_special_qdp(flw, out, ttt, -1, QDP_all, 1);
    //QDP_D_eq_r_times_D_plus_D(out, &mkappa, out, ttt, QDP_all);
    clov(flw, out, &mkappa, out, ttt, QDP_all);
  }

#endif
}


#endif
