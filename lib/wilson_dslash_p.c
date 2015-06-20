/* conventions for raw packed clover term

// -0.5*[cs*F01*S01+ct*F23*S23+cs*F12*S12+ct*F03*S03+cs*F02*S02+ct*F13*S13]
// Fmn = [ U_mu(x) U_nu(x+mu) U^d_mu(x+nu) U^d_nu(x) + refl. - herm ]/8

the clover term is a REAL array of length 72*QDP_sites_on_node
on each site are 72 REALs made of 2 packed Hermitian matrices of 36 REALs each
each packed Hermitian matrix has the following structure:
6 REALs for the diagonal
30 REALs = 15 Complex for the lower triangle stored columnwise

the physical spin and color component order for the clover term is
given by the following order:

color  spin
0      0
0      1
1      0
1      1
2      0
2      1
0      2
0      3
1      2
1      3
2      2
2      3

In other words each chiral block (6x6 matrix) is ordered with color slowest
and the 2 components of spin for that block going fastest.

Explicitly, if c(ic,is,jc,js) is the clover term with ic,is the color and spin
of the row and jc,js the color and spin of the column then the order of
elements is:

real(c(0,0, 0,0))
real(c(0,1, 0,1))
real(c(1,0, 1,0))
real(c(1,1, 1,1))
real(c(2,0, 2,0))
real(c(2,1, 2,1))
c(0,1, 0,0)   (1,0)   6
c(1,0, 0,0)   (2,0)   8
c(1,1, 0,0)   (3,0)  10
c(2,0, 0,0)   (4,0)  12
c(2,1, 0,0)   (5,0)  14
c(1,0, 0,1)   (2,1)  16
c(1,1, 0,1)   (3,1)  18
c(2,0, 0,1)   (4,1)  20
c(2,1, 0,1)   (5,1)  22
c(1,1, 1,0)   (3,2)  24
c(2,0, 1,0)   (4,2)  26
c(2,1, 1,0)   (5,2)  28
c(2,0, 1,1)   (4,3)  30
c(2,1, 1,1)   (5,3)  32
c(2,1, 2,0)   (5,4)  34

and likewise for the second block by adding 2 to all is,js values above.

for the first block [is,js=0,1]
cl[2*ic+is] = Re c(ic,is,ic,is)
(ic,is,jc,js) -> (2*ic+is,2*jc+js) = (i,j) -> 2nc + 2*((2nc-(j+1)/2)*j+i-j-1)
cl[2nc-2+(4nc-3-j)*j+2i] = c(ic,is,jc,js)

(is,js) = (0,0) -> i=2ic,j=2jc
 jc = 1 -> 2nc-2+14+4ic = 18+4ic
 jc = 2 -> 2nc-2+20+4ic = 24+4ic

00   4 +  4*ic + 14*jc
11  14 +  4*ic + 10*jc
10  12 + 10*ic +  4*jc
01   6 +  4*ic + (18-4*jc)*jc

*/

	//#define DO_TRACE
#include <string.h>
//#include <math.h>
#include <qop_internal.h>

#define CLOV_REALS (8*QLA_Nc*QLA_Nc) // 2 packed (2Nc)x(2Nc) Hermitian matrices
#define CLOV_SIZE (CLOV_REALS*sizeof(REAL)) 

//#define printf0 QOP_printf0
#define printf0(...)

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_wilson_inited;
extern int QOP_wilson_style;
extern int QOP_wilson_nsvec;
extern int QOP_wilson_nvec;
extern int QOP_wilson_cgtype;
extern int QOP_wilson_optnum;

static int old_style=-1;
static int old_optnum=-1;

#define NTMPSUB 2
#define NTMP (3*NTMPSUB)
#define NHTMP 16
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
    reset_temps(flw);				       \
  } \
  if( flw->dblstored != dblstore_style(QOP_wilson_style) ) { \
    double_store(flw); \
  } \
}

static void
free_temps(void)
{
  if(dslash_setup) {
    for(int i=0; i<NTMP; i++) {
      QDP_destroy_D(tin[i]);
    }
    if(shiftd_style(old_style)) {
      for(int i=0; i<NTMP; i++) {
	for(int j=0; j<NDTMP; j++) {
	  QDP_destroy_D(dtemp[i][j]);
	}
      }
    } else {
      for(int i=0; i<NTMP; i++) {
	for(int j=0; j<NHTMP; j++) {
	  QDP_destroy_H(htemp[i][j]);
	}
      }
    }
  }
  dslash_setup = 0;
}

static void
reset_temps(QOP_FermionLinksWilson *flw)
{
#define NC QDP_get_nc(flw->links[0])
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);

  TRACE;
  //printf("nc: %i\n", QLA_Nc);
  free_temps();
  for(int i=0; i<NTMP; i++) {
    tin[i] = QDP_create_D_L(lat);
  }
  if(shiftd_style(QOP_wilson_style)) {
    for(int i=0; i<NTMP; i++) {
      for(int j=0; j<NDTMP; j++) {
	dtemp[i][j] = QDP_create_D_L(lat);
      }
    }
  } else {
    for(int i=0; i<NTMP; i++) {
      for(int j=0; j<NHTMP; j++) {
	htemp[i][j] = QDP_create_H_L(lat);
      }
    }
  }
  dslash_setup = 1;
  old_style = QOP_wilson_style;
  old_optnum = QOP_wilson_optnum;
#undef NC
}

static void
double_store(QOP_FermionLinksWilson *flw)
{
#define NC QDP_get_nc(flw->links[0])
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);

  if(flw->dblstored) {
    for(int i=0; i<4; i++) {
      QDP_destroy_M(flw->bcklinks[i]);
    }
    flw->dblstored = 0;
  }
  if(dblstore_style(QOP_wilson_style)) {
    for(int i=0; i<4; i++) {
      flw->bcklinks[i] = QDP_create_M_L(lat);
    }
    for(int i=0; i<4; i++) {
      flw->dbllinks[2*i] = flw->links[i];
      flw->dbllinks[2*i+1] = flw->bcklinks[i];
    }
    QDP_ColorMatrix *m = QDP_create_M_L(lat);
    for(int i=0; i<4; i++) {
      QDP_M_eq_sM(m, flw->links[i], neighbor[i], QDP_backward, all);
      QDP_M_eq_Ma(flw->bcklinks[i], m, all);
    }
    QDP_destroy_M(m);
    flw->dblstored = dblstore_style(QOP_wilson_style);
  }
#undef NC
}

QDP_DiracFermion *
QOP_wilson_dslash_get_tmp(QOP_FermionLinksWilson *flw,
			  QOP_evenodd_t eo, int n)
{
#define NC QDP_get_nc(flw->links[0])
  check_setup(flw);
  if(n>=1 && n<=NTMPSUB) return tmpsub(eo,n);
  else return NULL;
#undef NC
}

static void
f_mu_nu(QDP_ColorMatrix *fmn, QLA_Real scale, QDP_ColorMatrix *link[],
	int mu, int nu)
{
#define NC QDP_get_nc(fmn)
  QDP_ColorMatrix *temp1,*temp2,*temp3,*temp4,*tmat4;
  QDP_ColorMatrix *pqt0,*pqt1,*pqt2,*pqt3;
  QDP_ColorMatrix *pqt4;
  QDP_Lattice *lat = QDP_get_lattice_M(fmn);
  QDP_Subset all = QDP_all_L(lat);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);

  temp1 = QDP_create_M_L(lat);
  temp2 = QDP_create_M_L(lat);
  temp3 = QDP_create_M_L(lat);
  temp4 = QDP_create_M_L(lat);
  tmat4 = QDP_create_M_L(lat);
  pqt0  = QDP_create_M_L(lat);
  pqt1  = QDP_create_M_L(lat);
  pqt2  = QDP_create_M_L(lat);
  pqt3  = QDP_create_M_L(lat);
  pqt4  = QDP_create_M_L(lat);

  int order_flag = 0;
  if(mu>nu) {
    int i = mu;
    mu = nu;
    nu = i;
    order_flag=1;
  }

  /* Get pqt0 = U_nu(x+mu) : U_nu(x) from mu direction    */
  /* Get pqt1 = U_mu(x+nu) : U_mu(x) from nu direction    */
  QDP_M_eq_sM(pqt0,link[nu],neighbor[mu],QDP_forward,all);
  QDP_M_eq_sM(pqt1,link[mu],neighbor[nu],QDP_forward,all);

  /* creating a corner : temp1 = U^d_nu(x) U_mu(x) */
  QDP_M_eq_Ma_times_M(temp1,link[nu],link[mu],all);

  /* ------------------------------------------------- */
  /* creating a corner : fmn = [pqt0]*[pqt1^d]         */
  /*                     fmn = U_nu(x+mu) U^d_mu(x+nu) */
  QDP_M_eq_M_times_Ma(fmn,pqt0,pqt1,all);

  /* Create the following loops */
  /* U^d_nu(x) U_mu(x) U_nu(x+mu) U^d_mu(x+nu)  */
  QDP_M_eq_M_times_M(temp2,temp1,fmn,all);

  /* U_nu(x+mu) U^d_mu(x+nu) U^d_nu(x) U_mu(x) */
  QDP_M_eq_M_times_M(temp3,fmn,temp1,all);

  /* Creating the following +mu -nu plaquette */
  QDP_M_eq_sM(pqt2,temp2,neighbor[nu],QDP_backward,all);

  /* Creating the following -mu +nu plaquette */
  QDP_M_eq_sM(pqt3,temp3,neighbor[mu],QDP_backward,all);
  //QDP_discard_M(temp3); /* data in temp3 is no longer needed */

  /* creating +mu +nu plaquette and pit it in fmn      */
  /* tmat4 = U_mu(x) U_nu(x+mu) U^d_mu(x+nu)           */
  /* fmn   = U_mu(x) U_nu(x+mu) U^d_mu(x+nu) U^d_nu(x) */ 
  QDP_M_eq_M_times_M(tmat4,link[mu],fmn,all);
  QDP_M_eq_M_times_Ma(fmn,tmat4,link[nu],all);

  /* What is left is +mu -nu plaquette and adding them up */
  /* Right hand side of the clover field */
  QDP_M_eq_M_plus_M(fmn,pqt2,fmn,all);

  /* tmat4 = [(pqt[1]^d)   ] * [(tempmat1 )            ]  */
  /*       = [U^d_mu(x+nu) ] * [U^dagger_nu(x) *U_mu(x)]  */
  /* temp2 = [tmat4]*[pqt0 ]                            */
  /*       = [U^d_mu(x+nu) U^d_nu(x) U_mu(x)*[U_nu(x+mu)] */
  QDP_M_eq_Ma_times_M(tmat4,pqt1,temp1,all);
  QDP_M_eq_M_times_M(temp2,tmat4,pqt0,all);
  /* temp1 is a result of a shift and won't be needed   */
  QDP_discard_M(temp1);

  /* temp2 is now a plaquette -mu -nu and must be gathered */
  /* with displacement -mu-nu */
  /* pqt4 = U^d_mu(x-nu)U^d(x-mu-nu)U_mu(x-mu-nu)U_nu(x-nu) */
  QDP_M_eq_sM(temp4,temp2,neighbor[mu],QDP_backward,all);
  QDP_M_eq_sM(pqt4,temp4,neighbor[nu],QDP_backward,all);
  QDP_discard_M(temp2);

  /* Now gather -mu +nu  plaquette and add to fmn    */
  /* f_mn was the right hand side of the clover term */
  /* I add the third plaquette : f_mn = f_mn+pqt3    */
  /* U_nu(x) U^d(x-mu+nu) U^d_nu(x-mu) U_mu(x-mu)    */
  QDP_M_eq_M_plus_M(fmn,fmn,pqt3,all);

  /* finally add the last plaquette        */
  /* fmn = fmn+ pqt4                       */
  /* pqt4 is the last plaquette missing    */
  /* This completes the 4-plaquette        */
  QDP_M_eq_M_plus_M(fmn,fmn,pqt4,all);

  /* F_munu  is now 1/8 of f_mn-f_mn^dagger */
  /* QDP_T_eqop_Ta(Type *r, Type *a, subset) : r=conjugate(a) */
  /* tmat4  =Hermitian(fmn) */
  QDP_M_eq_Ma(tmat4,fmn,all);

  if(order_flag ==0){
    QDP_M_eq_M_minus_M(tmat4,fmn,tmat4,all);
  }
  else {
    QDP_M_eq_M_minus_M(tmat4,tmat4,fmn,all);
  }

  /* F_mn = 1/8 *tmat4 */
  /* QDP_T_eq_r_times_T(Type *r, QLA_real *a, Type *b, subset); */
  scale *= 0.125;
  QDP_M_eq_r_times_M(fmn, &scale, tmat4, all);

  QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);
  QDP_destroy_M(tmat4);
  QDP_destroy_M(pqt0);
  QDP_destroy_M(pqt1);
  QDP_destroy_M(pqt2);
  QDP_destroy_M(pqt3);
  QDP_destroy_M(pqt4);
#undef NC
}

static void
get_clov(QLA_Real *clov, QDP_ColorMatrix *link[], QLA_Real cs, QLA_Real ct)
{
#define NC QDP_get_nc(link[0])
  QDP_ColorMatrix *a0, *a1, *b0, *b1;
  QDP_ColorMatrix *f_mn;
  QDP_ColorMatrix *xm,*ym;
  QDP_Lattice *lat = QDP_get_lattice_M(link[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  QDP_Subset all = QDP_all_L(lat);
  int nc = NC;
  int bsize = 4*nc*nc;

  f_mn  = QDP_create_M_L(lat);
  a0  = QDP_create_M_L(lat);
  a1  = QDP_create_M_L(lat);
  xm = QDP_create_M_L(lat);
  ym = QDP_create_M_L(lat);

  f_mu_nu(f_mn, cs, link, 0, 1);         /* f_mn = cs*F_{01}        */
  QDP_M_eq_M(xm, f_mn, all);         /* xm = cs*F_{01}        */
  QDP_M_eq_M(ym, f_mn, all);         /* ym = cs*F_{01}        */

  f_mu_nu(f_mn, ct, link, 2, 3);         /* f_mn = ct*F_{23}        */ 
  QDP_M_meq_M(xm, f_mn, all);        /* xm = cs*F_{01}-ct*F_{23} */
  QDP_M_peq_M(ym, f_mn, all);        /* ym = cs*F_{01}+ct*F_{23} */

  QDP_M_eq_i_M(a0, xm, all);         /* a[0] = i(cs*F_{01}-ct*F_{23}) */ 
  QDP_M_eq_i_M(a1, ym, all);         /* a[1] = i(cs*F_{01}+ct*F_{23}) */

  QLA_ColorMatrix(*A[2]);
  A[0] = QDP_expose_M(a0);
  A[1] = QDP_expose_M(a1);
  // is,js = 0,0 and 1,1
#pragma omp parallel for
  for(int s=0; s<sites_on_node; s++) {
    for(int b=0; b<2; b++) { // block
      int k0 = bsize*(2*s+b);
      for(int ic=0; ic<QLA_Nc; ic++) {
	for(int is=0; is<2; is++) {
	  int ii = 2*ic + is;
	  QLA_Real sgn = 1-2*is;
	  clov[k0+ii] = sgn*QLA_real(QLA_elem_M(A[b][s],ic,ic));
	  for(int jc=0; jc<ic; jc++) { 
	    int jj = 2*jc + is;
	    int kk = 2*(nc-1) + (4*nc-3-jj)*jj + 2*ii;
	    clov[k0+kk  ] = sgn*QLA_real(QLA_elem_M(A[b][s],ic,jc));
	    clov[k0+kk+1] = sgn*QLA_imag(QLA_elem_M(A[b][s],ic,jc));
	  }
	}
      }
    }
  }
  QDP_reset_M(a0);
  QDP_reset_M(a1);
  QDP_destroy_M(a0);
  QDP_destroy_M(a1);

  b0 = QDP_create_M_L(lat);
  b1 = QDP_create_M_L(lat);

  f_mu_nu(xm, cs, link, 1, 2);               /* xm   = cs*F_{12}         */  
  f_mu_nu(f_mn, ct, link, 0, 3);             /* f_mn = ct*F_{03}         */

  QDP_M_eq_M_plus_M(ym, xm, f_mn, all);  /* ym = cs*F_{12}+ct*F_{03}    */
  QDP_M_meq_M(xm, f_mn, all);            /* xm = cs*F_{12}-ct*F_{03}    */

  QDP_M_eq_i_M(b0, xm, all);             /* b0 = i(cs*F_{12}-ct*F_{03}) */
  QDP_M_eq_i_M(b1, ym, all);             /* b1 = i(cs*F_{12}+ct*F_{03}) */

  f_mu_nu(f_mn, cs, link, 0, 2);             /* f_mn = cs*F_{02}           */
  QDP_M_meq_M(b0, f_mn, all);          /* b[0] = b[0]-cs*F_{02}      */
  QDP_M_meq_M(b1, f_mn, all);          /* b[1] = b[1]-cs*F_{02}      */

  f_mu_nu(f_mn, ct, link, 1, 3);             /* f_mn = ct*F_{13}           */
  QDP_M_meq_M(b0, f_mn, all);          /* b[0] = b[0]-ct*F_{13}      */
  QDP_M_peq_M(b1, f_mn, all);          /* b[1] = b[1]+ct*F_{13}      */

  /*  b[0] = i(cs*F_{12}-ct*F_{03})-(cs*F_{02}+ct*F_{13})    */
  /*  b[1] = i(cs*F_{12}+ct*F_{03})-(cs*F_{02}-ct*F_{13})    */

  QLA_ColorMatrix(*B[2]);
  B[0] = QDP_expose_M(b0);
  B[1] = QDP_expose_M(b1);
  // is,js = 1,0 (or 0,1)
#pragma omp parallel for
  for(int s=0; s<sites_on_node; s++) {
    for(int b=0; b<2; b++) { // block
      int k0 = bsize*(2*s+b);
      for(int ic=0; ic<QLA_Nc; ic++) {
	int ii = 2*ic + 1;
	for(int jc=0; jc<QLA_Nc; jc++) {
	  int jj = 2*jc + 0;
	  if(ii<jj) { // conjugate
	    int kk = 2*(nc-1) + (4*nc-3-ii)*ii + 2*jj;
	    clov[k0+kk  ] =  QLA_real(QLA_elem_M(B[b][s],ic,jc));
	    clov[k0+kk+1] = -QLA_imag(QLA_elem_M(B[b][s],ic,jc));
	  } else {
	    int kk = 2*(nc-1) + (4*nc-3-jj)*jj + 2*ii;
	    clov[k0+kk  ] = QLA_real(QLA_elem_M(B[b][s],ic,jc));
	    clov[k0+kk+1] = QLA_imag(QLA_elem_M(B[b][s],ic,jc));
	  }
	}
      }
    }
  }
  QDP_reset_M(b0);
  QDP_reset_M(b1);
  QDP_destroy_M(b0);
  QDP_destroy_M(b1);

  QDP_destroy_M(f_mn);
  QDP_destroy_M(xm);
  QDP_destroy_M(ym);
#undef NC
}

static void
clov_unpack(NCPROT QLA_Complex u[2*QLA_Nc][2*QLA_Nc], REAL *p)
{
  int i, j, k;
  for(i=0; i<2*QLA_Nc; i++) {
    QLA_c_eq_r(u[i][i], p[i]);
  }
  k = 2*QLA_Nc;
  for(i=0; i<2*QLA_Nc; i++) {
    for(j=i+1; j<2*QLA_Nc; j++) {
      QLA_c_eq_ca(u[i][j], *((QLA_Complex *)(p+k)));
      QLA_c_eq_c(u[j][i], *((QLA_Complex *)(p+k)));
      k += 2;
    }
  }
}

static void
clov_pack(NCPROT REAL *p, QLA_Complex u[2*QLA_Nc][2*QLA_Nc])
{
  int i, j, k;
  for(i=0; i<2*QLA_Nc; i++) {
    p[i] = QLA_real(u[i][i]);
  }
  k = 2*QLA_Nc;
  for(i=0; i<2*QLA_Nc; i++) {
    for(j=i+1; j<2*QLA_Nc; j++) {
      QLA_Complex z;
      QLA_c_eq_ca(z, u[i][j]);
      QLA_c_peq_c(z, u[j][i]);
      QLA_c_eq_r_times_c(*((QLA_Complex *)(p+k)), 0.5, z);
      k += 2;
    }
  }
}

static void
clov_invert(NCPROT QLA_Complex ci[2*QLA_Nc][2*QLA_Nc], QLA_Complex c[2*QLA_Nc][2*QLA_Nc])
{
  int i, j, k;
  for(i=0; i<2*QLA_Nc; i++) {
    for(j=0; j<2*QLA_Nc; j++) {
      QLA_c_eq_r(ci[i][j], 0);
    }
    QLA_c_eq_r(ci[i][i], 1);
  }
  for(k=0; k<2*QLA_Nc; k++) {
    QLA_Complex s;
    //QLA_c_eq_r_div_c(s, 1, c[k][k]);
    {
      QLA_Real r;
      r = 1/QLA_norm2_c(c[k][k]);
      QLA_c_eq_r_times_c(s, r, c[k][k]);
      QLA_c_eq_ca(s, s);
    }
    for(j=0; j<2*QLA_Nc; j++) {
      QLA_Complex t;
      QLA_c_eq_c_times_c(t, s, c[k][j]);
      QLA_c_eq_c(c[k][j], t);
    }
    for(j=0; j<2*QLA_Nc; j++) {
      QLA_Complex t;
      QLA_c_eq_c_times_c(t, s, ci[k][j]);
      QLA_c_eq_c(ci[k][j], t);
    }
    for(i=0; i<2*QLA_Nc; i++) {
      if(i==k) continue;
      QLA_c_eq_c(s, c[i][k]);
      for(j=0; j<2*QLA_Nc; j++) {
	QLA_c_meq_c_times_c(c[i][j], s, c[k][j]);
      }
      for(j=0; j<2*QLA_Nc; j++) {
	QLA_c_meq_c_times_c(ci[i][j], s, ci[k][j]);
      }
    }
  }
}
/*
static void
printcm(QLA_Complex *m, int nr, int nc)
{
  int i, j;
  for(i=0; i<nr; i++) {
    for(j=0; j<nc; j++) {
      printf0("%f %f ", QLA_real(m[i*nc+j]), QLA_imag(m[i*nc+j]));
    }
    printf0("\n");
  }
}
*/

static void
get_clovinv(QOP_FermionLinksWilson *flw, REAL kappa)
{
#define NC QDP_get_nc(flw->links[0])
  QLA_Real m4 = 0.5/kappa;
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  int i, j;
  if(flw->clovinv==NULL)
    QOP_malloc(flw->clovinv, REAL, sites_on_node*CLOV_REALS);
  for(i=0; i<2*sites_on_node; i++) {
    QLA_Complex cu[2*QLA_Nc][2*QLA_Nc], ciu[2*QLA_Nc][2*QLA_Nc];
    clov_unpack(NCARG cu, flw->clov+(CLOV_REALS/2)*i);
    for(j=0; j<2*QLA_Nc; j++) QLA_c_peq_r(cu[j][j], m4);
    //if(i==0) printcm(flw->clov, QLA_Nc, 2*QLA_Nc);
    //if(i==0) printcm(cu, 2*QLA_Nc, 2*QLA_Nc);
    clov_invert(NCARG ciu, cu);
    clov_pack(NCARG flw->clovinv+(CLOV_REALS/2)*i, ciu);
    //if(i==0) printcm(cu, 2*QLA_Nc, 2*QLA_Nc);
    //if(i==0) printcm(ciu, 2*QLA_Nc, 2*QLA_Nc);
    //if(i==0) printcm(flw->clovinv, QLA_Nc, 2*QLA_Nc);
  }
  flw->clovinvkappa = kappa;
#undef NC
}

#define NC int nc
QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(QDP_Lattice *lat, REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
#undef NC
{
#define NC nc
  QOP_FermionLinksWilson *flw;
  QOP_GaugeField *gf;
  int sites_on_node = QDP_sites_on_node_L(lat);
  
  WILSON_INVERT_BEGIN;

  gf = QOP_create_G_from_raw(lat, links, evenodd);
  flw = QOP_wilson_convert_L_from_qdp(gf->links, NULL);

  if(clov!=NULL) {
    QOP_malloc(flw->clov, REAL, sites_on_node*CLOV_REALS);
    memcpy(flw->clov, clov, sites_on_node*CLOV_SIZE);
  } else {
    flw->clov = NULL;
  }

  /* Keep pointer to allocated space so it can be freed */
  flw->qopgf = gf;

  WILSON_INVERT_END;
  return flw;
#undef NC
}

QOP_FermionLinksWilson *
QOP_wilson_initialize_gauge_L(void)
{
  QOP_FermionLinksWilson *flw;

  WILSON_INVERT_BEGIN;

  QOP_malloc(flw          , QOP_FermionLinksWilson, 1);
  QOP_malloc(flw->links   , QDP_ColorMatrix *     , 4);
  QOP_malloc(flw->bcklinks, QDP_ColorMatrix *     , 4);
  QOP_malloc(flw->dbllinks, QDP_ColorMatrix *     , 8);

  flw->dblstored = 0;
  flw->clov      = NULL;
  flw->rawclov   = NULL;
  flw->clovinv   = NULL;
  flw->rawlinks  = NULL;
  flw->qopgf     = NULL;
  flw->gauge     = NULL;
  flw->qdpclov   = NULL;
  flw->eigcg.u   = NULL;

  WILSON_INVERT_END;

  return flw;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_info_t *info, 
			   QOP_wilson_coeffs_t *coeffs,
			   QOP_GaugeField *gauge)
{ 
#define NC QDP_get_nc(gauge->links[0])
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix        *newlinks[4];
  QDP_Lattice *lat = QDP_get_lattice_M(gauge->links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  QDP_Subset all = QDP_all_L(lat);

  WILSON_INVERT_BEGIN;

  /* initialize FermionLinksWilson and allocate memory---- */
  flw = QOP_wilson_initialize_gauge_L();

  /* First create QDP Color Matrices */
  for(int i=0; i<4; i++) {
    newlinks[i] = QDP_create_M_L(lat);
    QDP_M_eq_M(newlinks[i], gauge->links[i], all);
  }

  /* get the clover coefficients and put them in flw->clow */
  if(coeffs->clov_s != 0 || coeffs->clov_t != 0) {
    int nreals = sites_on_node * CLOV_REALS;
    QOP_malloc(flw->clov, REAL, nreals);
    get_clov(flw->clov, newlinks, 0.5*coeffs->clov_s, 0.5*coeffs->clov_t);
  }

  /* Check the anisotropy -------------------------------  */
  if(coeffs->aniso != 0. && coeffs->aniso != 1.) {
    for(int i=0; i<3; i++) {
      QLA_Real f = coeffs->aniso;
      QDP_M_eq_r_times_M(newlinks[i], &f, newlinks[i], all);
    }
  }

  /* Scale the links ------------------------------------- */
  for(int i=0; i<4; i++) {
    QLA_Real f    = -0.5;
    QDP_M_eq_r_times_M(newlinks[i], &f, newlinks[i], all);
  }

  /* newlinks go to flw->links --------------------------- */
  for(int i=0; i<4; i++) {
    flw->links[i] = newlinks[i];
  }

  flw->gauge = gauge;

  WILSON_INVERT_END;
  return flw;
#undef NC
}

#if QOP_Precision == 'F'

/* Create a single-precision copy of double-precision Wilson fermion links */
QOP_F_FermionLinksWilson *
QOP_FD_wilson_create_L_from_L(QOP_D_FermionLinksWilson *flw_double)
{
#define NC QDP_get_nc(flw_double->links[0])
  QOP_F_FermionLinksWilson *flw_single;
  QDP_Lattice *lat = QDP_D_get_lattice_M(flw_double->links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  QDP_Subset all = QDP_all_L(lat);

  /* Create the parent struct */
  flw_single = QOP_wilson_initialize_gauge_L();

  /* Copy scalar values */
  flw_single->dblstored = flw_double->dblstored;
  flw_single->clovinvkappa = flw_double->clovinvkappa;

  /* Create and copy the modified gauge links */
  if(flw_double->links != NULL) {
    //QOP_malloc(flw_single->links, QDP_F3_ColorMatrix *, 4);
    for(int i=0; i<4; i++) {
      flw_single->links[i] = QDP_create_M_L(lat);
      QDP_FD_M_eq_M(flw_single->links[i], flw_double->links[i], all);
    }
  } else {
    QOP_printf0("QOP_FD_wilson_create_L_from_L: Error: missing gauge links\n");
  }

  /* Create and copy the modified backward gauge links */
  //QOP_malloc(flw_single->bcklinks, QDP_F3_ColorMatrix *, 4);
  for(int i=0; i<4; i++) {
    if(flw_double->dblstored != 0 && flw_double->bcklinks[i] != NULL) {
      flw_single->bcklinks[i] = QDP_create_M_L(lat);
      QDP_FD_M_eq_M(flw_single->bcklinks[i], flw_double->bcklinks[i], all);
    } else {
      flw_double->bcklinks[i] = NULL;
    }
  }

  /* Create and copy the double links */
  //QOP_malloc(flw_single->dbllinks, QDP_F3_ColorMatrix *, 8);
  for(int i=0; i<4; i++) {
    if(flw_double->dblstored != 0) {
      flw_single->dbllinks[2*i] = flw_single->links[i];
      flw_single->dbllinks[2*i+1] = flw_single->bcklinks[i];
    } else {
      flw_single->dbllinks[2*i] = NULL;
      flw_single->dbllinks[2*i+1] = NULL;
    }
  }

  /* Create and copy the clover term */
  if(flw_double->clov != NULL) {
    int size = sites_on_node*CLOV_REALS;
    QOP_malloc(flw_single->clov, QLA_Real, size);
    for(int k=0; k<size; k++) {
      flw_single->clov[k] = flw_double->clov[k];
    }
  } else {
    flw_single->clov = NULL;
  }

  /* Create and copy clovinv term */
  if(flw_double->clovinv != NULL){
    int size = sites_on_node*CLOV_REALS;
    QOP_malloc(flw_single->clovinv, QLA_Real, size);
    for(int k=0; k<size; k++) {
      flw_single->clovinv[k] = flw_double->clovinv[k];
    }
  } else {
    flw_single->clovinv = NULL;
  }

  return flw_single;
#undef NC
}

#endif // QOP_Precision == 'F'

void
QOP_wilson_extract_L_to_raw(REAL *links[], REAL *clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_extract_L_to_raw unimplemented.");
  WILSON_INVERT_END;
}

void
QOP_wilson_destroy_L(QOP_FermionLinksWilson *flw)
{
  WILSON_INVERT_BEGIN;

  if(flw->qopgf) {
    QOP_destroy_G(flw->qopgf);
  } else {
    for(int i=0; i<4; i++) QDP_destroy_M(flw->links[i]);
  }
  free(flw->links);
  if(flw->dblstored) {
    for(int i=0; i<4; i++) QDP_destroy_M(flw->bcklinks[i]);
  }
  if(flw->qdpclov) {
    QDP_destroy_P(flw->qdpclov);
  }
  if(flw->clov) free(flw->clov);
  if(flw->clovinv) free(flw->clovinv);
  free(flw->bcklinks);
  free(flw->dbllinks);
  if(flw->eigcg.u) {
    for(int i=0; i<flw->eigcg.numax; i++) {
      QDP_destroy_D(flw->eigcg.u[i]);
    }
    free(flw->eigcg.u);
    free(flw->eigcg.l);
  }
  free(flw);
  WILSON_INVERT_END;
}

#define NC int nc
QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
#undef NC
{
#define NC nc
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_from_raw unimplemented");
  WILSON_INVERT_END;
  return NULL;
#undef NC
}

void
QOP_wilson_convert_L_to_raw(REAL ***links, REAL **clov,
			    QOP_FermionLinksWilson *src, QOP_evenodd_t evenodd)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_to_raw unimplemented");
  WILSON_INVERT_END;
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_G(QOP_info_t *info, QOP_wilson_coeffs_t *coeffs,
			    QOP_GaugeField *gauge)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_from_G unimplemented");
  WILSON_INVERT_END;
  return NULL;
}

QOP_GaugeField *
QOP_wilson_convert_L_to_G(QOP_FermionLinksWilson *links)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_to_G unimplemented");
  WILSON_INVERT_END;
  return NULL;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_qdp(QDP_ColorMatrix *links[],
			     QDP_DiracPropagator *clov)
{
#define NC QDP_get_nc(links[0])
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix *newlinks[4];
  QDP_Lattice *lat = QDP_get_lattice_M(links[0]);
  QDP_Subset all = QDP_all_L(lat);

  WILSON_INVERT_BEGIN;

  for(int i=0; i<4; i++) {
    newlinks[i] = QDP_create_M_L(lat);
    QDP_M_eq_M(newlinks[i], links[i], all);
  }

  flw = QOP_wilson_convert_L_from_qdp(newlinks, clov);

  WILSON_INVERT_END;
  return flw;
#undef NC
}

void
QOP_wilson_extract_L_to_qdp(QDP_ColorMatrix *links[],
			    QDP_DiracPropagator *clov,
			    QOP_FermionLinksWilson *src)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_extract_L_to_qdp unimplemented");
  WILSON_INVERT_END;
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_qdp(QDP_ColorMatrix *links[],
			      QDP_DiracPropagator *clov)
{
#define NC QDP_get_nc(clov)
  WILSON_INVERT_BEGIN;
  QOP_FermionLinksWilson *flw = QOP_wilson_initialize_gauge_L();
  QDP_Lattice *lat = QDP_get_lattice_M(links[0]);
  int sites_on_node = QDP_sites_on_node_L(lat);
  QDP_Subset all = QDP_all_L(lat);

  if(clov!=NULL) {
    int size = sites_on_node*CLOV_REALS;
    QOP_malloc(flw->clov, REAL, size);
    QLA_DiracPropagator(*dp);
    dp = QDP_expose_P(clov);
#pragma omp parallel for
    for(int x=0; x<sites_on_node; x++) {
      for(int b=0; b<2; b++) { // two chiral blocks
	int k = (2*x+b)*(QLA_Nc*(2*QLA_Nc+1));
	// first the diagonal
	for(int i=0; i<2*QLA_Nc; i++) {
	  int ic = i/2;
	  int is = 2*b + i%2;
	  flw->clov[k++] = QLA_real(QLA_elem_P(dp[x], ic, is, ic, is));
	}
	// now the offdiagonal
	for(int i=0; i<2*QLA_Nc; i++) {
	  int ic = i/2;
	  int is = 2*b + i%2;
	  for(int j=i+1; j<2*QLA_Nc; j++) {
	    QLA_Complex z1, z2;
	    int jc = j/2;
	    int js = 2*b + j%2;
	    //QLA_c_eq_c_plus_ca(z1, QLA_elem_P(dp[x], ic, is, jc, js),
	    //                   QLA_elem_P(dp[x], jc, js, ic, is));
	    QLA_c_eq_c(z1, QLA_elem_P(dp[x], ic, is, jc, js));
	    QLA_c_peq_ca(z1, QLA_elem_P(dp[x], jc, js, ic, is));
	    QLA_c_eq_r_times_c(z2, 0.5, z1);
	    flw->clov[k++] = QLA_real(z2);
	    flw->clov[k++] = -QLA_imag(z2); // - since we now store lower tri
	  }
	}
      }
    }
    QDP_reset_P(clov);
  } else {
    flw->clov = NULL;
  }
  flw->clovinv = NULL;

  flw->dblstored = 0;
  for(int i=0; i<4; i++) {
    flw->links[i] = links[i];
  }
  // scale links
  for(int i=0; i<4; i++) {
    QLA_Real f = -0.5;
    QDP_M_eq_r_times_M(flw->links[i], &f, flw->links[i], all);
  }

  //check_setup(flw);

  WILSON_INVERT_END;
  return flw;
#undef NC
}

void
QOP_wilson_convert_L_to_qdp(QDP_ColorMatrix ***links,
			    QDP_DiracPropagator **clov,
			    QOP_FermionLinksWilson *src)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_to_qdp unimplemented");
  WILSON_INVERT_END;
}


/********************/
/* Dslash functions */
/********************/

static void
clov(QOP_FermionLinksWilson *flw, REAL kappa, QDP_DiracFermion *out,
     QDP_DiracFermion *in, QDP_Subset subset, int peq);

static void
apply_clov(REAL *clov, QLA_Real m4, QDP_DiracFermion *out,
	   QDP_DiracFermion *in, QDP_Subset subset);

static void
wilson_dslash0(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n);

static void
wilson_dslash1(QOP_FermionLinksWilson *flw,
	       QDP_DiracFermion *dest, QDP_DiracFermion *src,
	       int sign, QOP_evenodd_t eo, int n);

#define wilson_hop(flw, dest, src, sign, eo,lat)	\
{ \
  QDP_DiracFermion *tsrc = src; \
  int _n = 1; \
  while(1) { \
    if(src==tmpsub(eo,_n)) break; \
    if(_n==NTMPSUB) { \
      _n = 1; \
      tsrc = tmpsub(eo,_n); \
      QDP_D_eq_D(tsrc, src, qdpsub(oppsub(eo),lat));	\
      break; \
    } \
    _n++; \
  } \
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
#define NC QDP_get_nc(flw->links[0])
  TRACE;
  check_setup(flw);
  TRACE;
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  QDP_Subset even = QDP_even_L(lat);
  QDP_Subset odd = QDP_odd_L(lat);

  if(eo_in==eo_out) {
    if(eo_out==QOP_EVENODD) {
      wilson_hop(flw, out, in, sign, QOP_EVENODD, lat);
      clov(flw, kappa, out, in, all, 1);
    } else if(eo_out==QOP_EVEN) {
      clov(flw, kappa, out, in, even, 0);
    } else {
      clov(flw, kappa, out, in, odd, 0);
    }
  } else {
    if(eo_out==QOP_EVEN || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_ODD) {
	wilson_hop(flw, out, in, sign, QOP_EVEN, lat);
      } else if(eo_in==QOP_EVEN) {
	clov(flw, kappa, out, in, even, 0);
      } else {
	wilson_hop(flw, out, in, sign, QOP_EVEN, lat);
	clov(flw, kappa, out, in, even, 1);
      }
    }
    if(eo_out==QOP_ODD || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_EVEN) {
	wilson_hop(flw, out, in, sign, QOP_ODD, lat);
      } else if(eo_in==QOP_ODD) {
	clov(flw, kappa, out, in, odd, 0);
      } else {
	wilson_hop(flw, out, in, sign, QOP_ODD, lat);
	clov(flw, kappa, out, in, odd, 1);
      }
    }
  }
#undef NC
}

void
QOP_wilson_diaginv(QOP_info_t *info,
		   QOP_FermionLinksWilson *flw,
		   REAL kappa,
		   QOP_DiracFermion *out,
		   QOP_DiracFermion *in,
		   QOP_evenodd_t eo)
{
  QOP_wilson_diaginv_qdp(info, flw, kappa, out->df, in->df, eo);
}

void
QOP_wilson_diaginv_qdp(QOP_info_t *info,
		       QOP_FermionLinksWilson *flw,
		       REAL kappa,
		       QDP_DiracFermion *out,
		       QDP_DiracFermion *in,
		       QOP_evenodd_t eo)
{
#define NC QDP_get_nc(flw->links[0])
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  if(flw->clov==NULL) {
    QLA_Real f = 2*kappa;
    QDP_D_eq_r_times_D(out, &f, in, qdpsub(eo,lat));
  } else {
    if( flw->clovinv==NULL || flw->clovinvkappa!=kappa ) {
      get_clovinv(flw, kappa);
    }
    QDP_D_eq_zero(out, qdpsub(eo,lat));
    apply_clov(flw->clovinv, 0, out, in, qdpsub(eo,lat));
  }
#undef NC
}

#define cmplx(x) (*((QLA_Complex *)(&(x))))

static void
apply_clov_qla(NCPROT REAL *clov, QLA_Real m4,
	       QLA_DiracFermion *restrict clov_out,
	       QLA_DiracFermion *restrict clov_in,
	       QDP_Subset subset, QDP_Subset even, QDP_Subset odd)
{
  int start, end;
  if(subset==odd) start = QDP_subset_len(even);
  else start = 0;
  end = start + QDP_subset_len(subset);
  int nc2 = 2*QLA_Nc;
#pragma omp parallel for
  for(int x=start; x<end; x++) {
    for(int b=0; b<2; b++) {
      int xb = nc2*nc2*(2*x+b);  // chiral block offset (in REALs)
      int xbo = xb+nc2;
#define clov_diag(i) clov[xb+i]
#define clov_offd(i) cmplx(clov[xbo+2*i])
#define clov_offd2(i) cmplx(clov[xbo+i])
#define src(i) QLA_elem_D(clov_in[x],i/2,2*b+(i&1))
#define dest(i) QLA_elem_D(clov_out[x],i/2,2*b+(i&1))
      //flops = 6*44 = 264; bytes = 4*(36+12+12) = 240
#if QLA_Colors != 3
      {
	for(int i=0; i<nc2; i++) {
	  QLA_c_eq_r_times_c(dest(i), m4+clov_diag(i), src(i));
	}
	for(int j=0; j<nc2; j++) {
	  for(int i=j+1; i<nc2; i++) {
	    //cl[2nc-2+(4nc-3-j)*j+2i] = c(ic,is,jc,js)
	    int k = (4*QLA_Nc-3-j)*j + 2*(i-1);
	    QLA_c_peq_c_times_c(dest(i), clov_offd2(k), src(j));
	    QLA_c_peq_ca_times_c(dest(j), clov_offd2(k), src(i));
	  }
	}
      }
#else
      {
	QLA_Complex z0, z1, z2, z3, z4, z5;

	QLA_c_eq_r_times_c(z0, m4+clov_diag(0), src(0));
	QLA_c_eq_r_times_c(z1, m4+clov_diag(1), src(1));
	QLA_c_eq_r_times_c(z2, m4+clov_diag(2), src(2));
	QLA_c_eq_r_times_c(z3, m4+clov_diag(3), src(3));
	QLA_c_eq_r_times_c(z4, m4+clov_diag(4), src(4));
	QLA_c_eq_r_times_c(z5, m4+clov_diag(5), src(5));

	QLA_c_peq_c_times_c(z1, clov_offd(0), src(0));
	QLA_c_peq_ca_times_c(z0, clov_offd(0), src(1));

	QLA_c_peq_c_times_c(z3, clov_offd(9), src(2));
	QLA_c_peq_ca_times_c(z2, clov_offd(9), src(3));

	QLA_c_peq_c_times_c(z5, clov_offd(14), src(4));
	QLA_c_peq_ca_times_c(z4, clov_offd(14), src(5));

	QLA_c_peq_c_times_c(z2, clov_offd(1), src(0));
	QLA_c_peq_ca_times_c(z0, clov_offd(1), src(2));

	QLA_c_peq_c_times_c(z4, clov_offd(12), src(3));
	QLA_c_peq_ca_times_c(z3, clov_offd(12), src(4));

	QLA_c_peq_c_times_c(z5, clov_offd(8), src(1));
	QLA_c_peq_ca_times_c(z1, clov_offd(8), src(5));

	QLA_c_peq_c_times_c(z3, clov_offd(2), src(0));
	QLA_c_peq_ca_times_c(z0, clov_offd(2), src(3));

	QLA_c_peq_c_times_c(z4, clov_offd(7), src(1));
	QLA_c_peq_ca_times_c(z1, clov_offd(7), src(4));

	QLA_c_peq_c_times_c(z5, clov_offd(11), src(2));
	QLA_c_peq_ca_times_c(z2, clov_offd(11), src(5));

	QLA_c_peq_c_times_c(z4, clov_offd(3), src(0));
	QLA_c_peq_ca_times_c(z0, clov_offd(3), src(4));

	QLA_c_peq_c_times_c(z5, clov_offd(13), src(3));
	QLA_c_peq_ca_times_c(z3, clov_offd(13), src(5));

	QLA_c_peq_c_times_c(z2, clov_offd(5), src(1));
	QLA_c_peq_ca_times_c(z1, clov_offd(5), src(2));

	QLA_c_peq_c_times_c(z5, clov_offd(4), src(0));
	QLA_c_peq_ca_times_c(z0, clov_offd(4), src(5));

	QLA_c_peq_c_times_c(z4, clov_offd(10), src(2));
	QLA_c_peq_ca_times_c(z2, clov_offd(10), src(4));

	QLA_c_peq_c_times_c(z3, clov_offd(6), src(1));
	QLA_c_peq_ca_times_c(z1, clov_offd(6), src(3));

	QLA_c_peq_c(dest(0), z0);
	QLA_c_peq_c(dest(1), z1);
	QLA_c_peq_c(dest(2), z2);
	QLA_c_peq_c(dest(3), z3);
	QLA_c_peq_c(dest(4), z4);
	QLA_c_peq_c(dest(5), z5);
      }
#endif // QLA_Colors == 3
    }
  }
}

static void
apply_clov(REAL *clov, QLA_Real m4, QDP_DiracFermion *out,
	   QDP_DiracFermion *in, QDP_Subset subset)
{
#define NC QDP_get_nc(in)
  QDP_Lattice *lat = QDP_get_lattice_D(in);
  QDP_Subset even = QDP_even_L(lat);
  QDP_Subset odd = QDP_odd_L(lat);
  QLA_DiracFermion *clov_out, *clov_in;
  clov_out = QDP_expose_D(out);
  clov_in = QDP_expose_D(in);
  apply_clov_qla(NCARG clov, m4, clov_out, clov_in, subset, even, odd);
  QDP_reset_D(in);
  QDP_reset_D(out);
#undef NC
}

static void
clov(QOP_FermionLinksWilson *flw, REAL kappa, QDP_DiracFermion *out,
     QDP_DiracFermion *in, QDP_Subset subset, int peq)
{
  QLA_Real m4 = 0.5/kappa;
  if(flw->clov==NULL) {
    if(peq) {
      QDP_D_peq_r_times_D(out, &m4, in, subset);
    } else {
      QDP_D_eq_r_times_D(out, &m4, in, subset);
    }
  } else {
    if(!peq) {
      QDP_D_eq_zero(out, subset);
    }
    apply_clov(flw->clov, m4, out, in, subset);
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
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  QDP_Subset even = QDP_even_L(lat);
  QDP_Subset odd = QDP_odd_L(lat);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);

  subset = qdpsub(eo,lat);
  othersubset = qdpsub(oppsub(eo),lat);
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

  if(subset==even) othersubset = odd;
  else if(subset==odd) othersubset = even;
  else othersubset = all;

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  printf0("dslash0\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_D_veq_sD\n");
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, neighbor+mu, fwd+mu, subset,
		   QOP_wilson_nsvec);
      printf0("end QDP_D_veq_sD\n");
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_D\n");
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
			 othersubset, QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, neighbor+mu, fwd+mu,
		   subset, QOP_wilson_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  printf0("dslash0 - back\n");
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
    }
  } else {
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      printf0("QDP_H_veq_spproj_Ma_times_D\n");
      QDP_H_veq_spproj_Ma_times_D(htemp[ntmp]+12+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+12+mu, neighbor+mu,
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
    for(mu=0; mu<4; mu+=QOP_wilson_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->links+mu, dtemp[ntmp]+mu,
				  dir+mu, sgn+mu, subset, QOP_wilson_nvec);
    }
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
  TRACE;
  int mu, ntmp;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];
  QDP_Subset subset, othersubset;
  QDP_Lattice *lat = QDP_get_lattice_D(src);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);

  subset = qdpsub(eo,lat);
  othersubset = qdpsub(oppsub(eo),lat);
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
    sh[2*mu] = neighbor[mu];
    sh[2*mu+1] = neighbor[mu];
    sd[2*mu] = QDP_forward;
    sd[2*mu+1] = QDP_backward;
  }

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */

  TRACE;
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

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add to dest */

  TRACE;
  QDP_D_eq_zero(dest, subset);
  TRACE;
  if(shiftd_style(QOP_wilson_style)) {
    //QOP_wilson_nvec = 1;
    //printf0("%p %p %p %p %p %p %i\n", vdest, flw->dbllinks, dtemp[ntmp], dir, sgn, subset, QOP_wilson_nvec);
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      //TRACE;
      //printf0("%i %p %p %p %i %i\n", mu, vdest[mu], flw->dbllinks[mu], dtemp[ntmp][mu], dir[mu], sgn[mu]);
      //QDP_D_eq_zero(vdest[mu], subset);
      //TRACE;
      //QDP_M_eq_zero(flw->dbllinks[mu], subset);
      //QDP_ColorMatrix *m = QDP_create_M_L(lat);
      //QDP_M_eq_M(m, flw->dbllinks[mu], subset);
      //QDP_destroy_M(m);
      //TRACE;
      //QDP_D_eq_D(vdest[mu], dtemp[ntmp][mu], subset);
      //TRACE;
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->dbllinks+mu, dtemp[ntmp]+mu,
				  dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      //QDP_D_vpeq_spproj_M_times_D(vdest+mu, flw->dbllinks+mu, dtemp[ntmp],
      //		          dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      //TRACE;
    }
  } else {
    for(mu=0; mu<8; mu+=QOP_wilson_nvec) {
      TRACE;
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, flw->dbllinks+mu, htemp[ntmp]+mu,
				   dir+mu, sgn+mu, subset, QOP_wilson_nvec);
      TRACE;
    }
  }

  TRACE;
  if(shiftd_style(QOP_wilson_style)) {
    for(mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for(mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
  TRACE;
}

#if QOP_Colors == 3

// ---------------------------------------------------------- :
// new fermilab action IFLA -- added by bugra --------------- :
// New name is Ok-action ------------------------------------ :
// ---------------------------------------------------------- :
// ---------------------------------------------------------- :
static void
wilson_okaction_full_testtadpole(QOP_FermionLinksWilson *flw,
                                 QDP_DiracFermion *dest,
                                 QDP_DiracFermion *src,
                                 int sign,
                                 QOP_wilson_ifla_coeffs_t *coeffs,
                                 QOP_evenodd_t eo, int n);
void
QOP_wilson_ifla_dslash(QOP_info_t *info,
                       QOP_FermionLinksWilson *flw,
                       REAL kappa,
                       int sign,
                       QOP_wilson_ifla_coeffs_t *coeffs,
                       QOP_DiracFermion *out,
                       QOP_DiracFermion *in,
                       QOP_evenodd_t eo_out,
                       QOP_evenodd_t eo_in)
{
  //printf("This is the function : (QOP_wilson_ifla_dslash)\n");
  //printf("Currently being tested **\n");
  QOP_wilson_ifla_dslash_qdp(info,flw,kappa,sign,coeffs,out->df,in->df,eo_out,eo_in);
};

void QOP_wilson_ifla_dslash_qdp(QOP_info_t *info,
                                QOP_FermionLinksWilson *flw,
                                REAL kappa,
                                int  sign,
                                QOP_wilson_ifla_coeffs_t *coeffs,
                                QDP_DiracFermion *out,
                                QDP_DiracFermion *in,
                                QOP_evenodd_t eo_out,
                                QOP_evenodd_t eo_in)
{
#if 0
  // test version
  wilson_ifla_full_tad(flw,out,in,sign,coeffs,2,2);
  //wilson_okaction_full_4(flw,out,in,sign,coeffs,2,2);
#else
  // new version correct c_5 term
  wilson_okaction_full_testtadpole(flw,out,in,sign,coeffs,2,2);
#endif
};

// ****************************************************************
// ****************************************************************
//    _____.   -----.  .      .-----.  -----. .    .  |
//  ||     |  |        |      |     | |       |    |  |     
//  ||     |--.-----.  |      |-----| .-----. |----|  .-----.
//  ||     |        |  |      |     |       | |    |  |     |
//  ||-----.  .-----.   .----- .     . .-----. .    .  .-----.
// ****************************************************************
// ****************************************************************

static void
wilson_okaction_full_testtadpole(QOP_FermionLinksWilson *flw,
				 QDP_DiracFermion *result,
				 QDP_DiracFermion *source,
				 int sign,
				 QOP_wilson_ifla_coeffs_t *coeffs,
				 QOP_evenodd_t eo,
				 int n)
{

  //printf("LATEST V.4 FULL OK_ACTION\n");
  /* Tadpole improved version   */
  /* Here is the tadpole factor */
  /* Hard wired for now         */
  /* -------------------------- */
  static int dir;
  /* -------------------------- */
  QLA_Real tpf      = 1.0/(coeffs->u0);
  //warning: unused variable 'sff' ??
  //QLA_Real sff      = 1e10; // comment out for now to avoid warnings

  //printf("Tadpole Factor u0 = %e\n",tpf);
  /* Coefficients ------------- */
  QLA_Real Kappa   = coeffs->kapifla;
  QLA_Real r_s     = coeffs->r_s;
  QLA_Real zeta    = coeffs->zeta;
  QLA_Real c_E     = coeffs->c_E;
  QLA_Real c_B     = coeffs->c_B;
  /* Improvement Coefficients */
  QLA_Real c_1     = coeffs->c_1;
  QLA_Real c_2     = coeffs->c_2;
  QLA_Real c_3     = coeffs->c_3;
  QLA_Real c_4     = coeffs->c_4;
  QLA_Real c_5     = coeffs->c_5;
  QLA_Real c_EE    = coeffs->c_EE;
  /* ------------------------ */


  QLA_Real slc1    = 0.5;
  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  QDP_Subset even = QDP_even_L(lat);
  QDP_Subset odd = QDP_odd_L(lat);
  QDP_Shift *neighbor = QDP_neighbor_L(lat);

#if 0
  printf("&&IFLA COEFFICIENTS : \n");
  printf("kappa_wilson = %e\t r_s = %e\t zeta = %e\n",Kappa,r_s,zeta);
  printf("c_E          = %e\t c_B = %e\t c_EE = %e\n",c_E,c_B,c_EE);
  printf("c_1          = %e\t c_2 = %e\t c_3  = %e\n",c_1,c_2,c_3);
  printf("c_4          = %e\t c_5 = %e\t          \n",c_4,c_5);
#endif


  /* Gamma Matrix Indices in QDP ============================ */

  //warning: variable 'bidx' set but not used ??
  QLA_Int gidx[4],aidx[4]/*,bidx[4]*/;  // comment out for now

  gidx[0] = 1;  /* gamma[X]   = QDP_Gamma[1]   */
  gidx[1] = 2;  /* gamma[Y]   = QDP_Gamma[2]   */
  gidx[2] = 4;  /* gamma[Z]   = QDP_Gamma[4]   */
  gidx[3] = 8;  /* gamma[0]   = QDP_Gamma[8]   */

  aidx[0] = 9;  /* alpha[0]   = -QDP_Gamma[9]  */
  aidx[1] = 10; /* alpha[1]   = -QDP_Gamma[10] */
  aidx[2] = 12; /* alpha[2]   = -QDP_Gamma[12] */

#if 0 // comment out unused variables
  bidx[0] = 3;  /* i*Sigma[0] = QDP_Gamma[3]   */
  bidx[1] = 5;  /* i*Sigma[1] = QDP_Gamma[5]   */
  bidx[2] = 6;  /* i*Simga[2] = QDP_Gamma[6]   */
#endif
  /* --------------------------------------------------------- */

  /* DUMMY FIELDS ============================================ */
  int mu;
  QDP_DiracFermion *tempD1,*tempD2,*tempD3,*tempD4;
  QDP_DiracFermion *tempM,*tempP;
  QDP_DiracFermion *gamD;
  QDP_ColorMatrix  *tempG1,*tempG2,*tempG3;
  QDP_DiracFermion *psi_up[4];
  QDP_DiracFermion *psi_dw[4];
  QDP_ColorMatrix *f12,*f02,*f01;
  

  tempD1       = QDP_create_D_L(lat);
  tempD2       = QDP_create_D_L(lat);
  tempD3       = QDP_create_D_L(lat);
  tempD4       = QDP_create_D_L(lat);
  tempP        = QDP_create_D_L(lat);
  tempM        = QDP_create_D_L(lat);
  gamD         = QDP_create_D_L(lat);
  tempG1       = QDP_create_M_L(lat);
  tempG2       = QDP_create_M_L(lat);
  tempG3       = QDP_create_M_L(lat);
  f01          = QDP_create_M_L(lat);
  f02          = QDP_create_M_L(lat);
  f12          = QDP_create_M_L(lat);

  /* --------------------------------------------------------- */
  
  /* Reading the gauge fields from FermionLinksWilson */
  QDP_ColorMatrix  *gauge[4];

  //  Rescaling the fields back
  static QLA_Real scale1 = -2.0;
  for(mu = 0;mu<4; mu++){
    gauge[mu] = QDP_create_M_L(lat);
    //gauge[mu] = flw->links[mu];
    QDP_M_eq_r_times_M(gauge[mu],&scale1,flw->links[mu],all);
    // Here I scale with the tadpole factor 
    /*
    QDP_M_eq_r_times_M(gauge[mu],&tpf,gauge[mu],all);
    */
  };

  //Kappa = Kappa*tpf;

  /* This is 1/[ 2*Kappa_Wilson ] */
  QLA_Real m4 = 0.5/Kappa;

  //printf("1/(2*Kappa) = %e\n",m4);

  /*^^ result = (1/2*Kappa)*Source  */
  QDP_D_eq_zero(result,all);
  QDP_D_eq_zero(gamD,all);

#if 0
  static QLA_DiracFermion *qladf;
  static int st,nc,ns;
  static QLA_Real c1,c2;
  qladf = QDP_expose_D(source);
  printf("Exposing Dirac Fermion :: source \n");  
  for(st=1;st<2;st++){
    for(nc=0;nc<3;nc++){
      for(ns=0;ns<1;ns++){
	c1=QLA_real(QLA_elem_D(qladf[st],nc,ns));
	c2=QLA_imag(QLA_elem_D(qladf[st],nc,ns));
	printf("site : %d  color : %d  spin  : %d = (%e,%e)\n",st,nc,ns,c1,c2);
      };
    };
    printf("***********************************\n");
  };
  QDP_reset_D(source);
#endif
  // ------------------------------------------------
  //QDP_D_eq_r_times_D(source,&sff,source,all);



  //QLA_Real kapifla = m4+18.0*c_4; definition!!
  QDP_D_eq_r_times_D(result,&m4,source,all);
  
  /*^^ result -= [ D-Slash ]*Source */
  for(dir=0;dir<4;dir++){
    psi_up[dir] = QDP_create_D_L(lat);
    psi_dw[dir] = QDP_create_D_L(lat);
  };

  QLA_Real kat_a = 0.5*r_s*zeta+4.0*c_4;
  QLA_Real kat_b = 0.5*zeta-(c_1)-(6.0*c_2);
  
  
  //printf("r^primes_s = %e\n", kat_a);
  //printf("zeta^prime = %e\n", kat_b);
  

  for(dir=0;dir<4;dir++){

    QDP_D_eq_M_times_sD(psi_up[dir],gauge[dir],source,neighbor[dir],QDP_forward,all);
    QDP_D_eq_Ma_times_D(tempD1,gauge[dir],source,all);
    QDP_D_eq_sD(psi_dw[dir],tempD1,neighbor[dir],QDP_backward,all);
    
    if(dir!=3){
      QDP_D_eq_D_plus_D(tempD1,psi_up[dir],psi_dw[dir],all);
      QDP_D_eq_r_times_D(tempD1,&kat_a,tempD1,all);
      
      QDP_D_eq_D_minus_D(tempD2,psi_up[dir],psi_dw[dir],all);
      QDP_D_eq_gamma_times_D(tempD3,tempD2,gidx[dir],all);
      //
      if(c_3!=0){
	QDP_D_peq_D(gamD,tempD3,all); // \Gamma.D.Psi(x) for c_3
      };
      //
      QDP_D_eq_r_times_D(tempD3,&kat_b,tempD3,all);
      QDP_D_meq_D(tempD1,tempD3,all);
    }
    else {
      /* This is the time component */
      QDP_D_eq_D_plus_D( tempD1,psi_up[dir],psi_dw[dir],all);
      QDP_D_eq_D_minus_D(tempD2,psi_up[dir],psi_dw[dir],all);
      QDP_D_eq_gamma_times_D(tempD3,tempD2,gidx[dir],all);
      QDP_D_meq_D(tempD1,tempD3,all);
      QDP_D_eq_r_times_D(tempD1,&slc1,tempD1,all);
    };
    QDP_D_meq_D(result,tempD1,all);
  };

  
  /* ---------- up to here is standard -------------- */
 
  /* ------------------------------------------------ */
  /* ------- start of c_1+c_2+c4 diagonal terms ----- */
  /* ------------------------------------------------ */
  if((c_1!=0)||(c_2!=0)||(c_4!=0)){

    QLA_Real qc     = (0.5*c_1+c_2)*tpf;
    
    QDP_DiracFermion *u0u0,*u1u1,*u2u2,*d0d0,*d1d1,*d2d2;
    
    u0u0 = QDP_create_D_L(lat);
    u1u1 = QDP_create_D_L(lat);
    u2u2 = QDP_create_D_L(lat);
    d0d0 = QDP_create_D_L(lat);
    d1d1 = QDP_create_D_L(lat);
    d2d2 = QDP_create_D_L(lat);
    
    if((c_1!=0)||(c_2!=0)){
      // Create U_i(x)U_j(x+i)Psi(x+i+j)
      QDP_D_eq_M_times_sD(u0u0,gauge[0],psi_up[0],neighbor[0],QDP_forward,all);
      QDP_D_eq_M_times_sD(u1u1,gauge[1],psi_up[1],neighbor[1],QDP_forward,all);
      QDP_D_eq_M_times_sD(u2u2,gauge[2],psi_up[2],neighbor[2],QDP_forward,all);
      // Create Ud_i(x-i)Ud_j(x-i-j)Psi(x-i-j)
      QDP_D_eq_Ma_times_D(tempD1,gauge[0],psi_dw[0],all);
      QDP_D_eq_sD(d0d0,tempD1,neighbor[0],QDP_backward,all);
      QDP_D_eq_Ma_times_D(tempD1,gauge[1],psi_dw[1],all);
      QDP_D_eq_sD(d1d1,tempD1,neighbor[1],QDP_backward,all);
      QDP_D_eq_Ma_times_D(tempD1,gauge[2],psi_dw[2],all);
      QDP_D_eq_sD(d2d2,tempD1,neighbor[2],QDP_backward,all);
      // Construct and Add :
      QDP_D_eq_D_minus_D(tempD1,u0u0,d0d0,all);
      QDP_D_eq_gamma_times_D(tempD2,tempD1,gidx[0],all);
      QDP_D_eq_D_minus_D(tempD1,u1u1,d1d1,all);
      QDP_D_eq_gamma_times_D(tempD3,tempD1,gidx[1],all);
      QDP_D_eq_D_plus_D(tempP,tempD2,tempD3,all);
      QDP_D_eq_D_minus_D(tempD1,u2u2,d2d2,all);
      QDP_D_eq_gamma_times_D(tempD2,tempD1,gidx[2],all);
      QDP_D_peq_D(tempP,tempD2,all);
      QDP_D_eq_r_times_D(tempP,&qc,tempP,all);
      QDP_D_peq_D(result,tempP,all);
    };
    //  
    if(c_4!=0){
      QLA_Real tadc4 = c_4*tpf;
      QDP_D_eq_D_plus_D(tempD1,u0u0,d0d0,all);
      QDP_D_eq_D_plus_D(tempD2,u1u1,d1d1,all);
      QDP_D_eq_D_plus_D(tempP,tempD1,tempD2,all);
      QDP_D_eq_D_plus_D(tempD1,u2u2,d2d2,all);
      QDP_D_peq_D(tempP,tempD1,all);
      QDP_D_eq_r_times_D(tempP,&tadc4,tempP,all);
      QDP_D_peq_D(result,tempP,all);
    };
    
    QDP_destroy_D(u0u0);
    QDP_destroy_D(u1u1);
    QDP_destroy_D(u2u2);
    QDP_destroy_D(d0d0);
    QDP_destroy_D(d1d1);
    QDP_destroy_D(d2d2);
    
  };
  /* ------------------------------------------------ */
  /* --------- end of c_1+c_2+c4 diagonal terms ----- */
  /* ------------------------------------------------ */


  /* ------------------------------------------------ */
  /* ------- start of c_2 non-diagonal terms -------- */
  /* ------------------------------------------------ */

  if(c_2!=0){
    QLA_Real c2half = (0.5*c_2)*tpf;
    //1
    QDP_DiracFermion *s0,*s1,*s2,*d0,*d1,*d2;
    s0 = QDP_create_D_L(lat);
    s1 = QDP_create_D_L(lat);
    s2 = QDP_create_D_L(lat);
    d0 = QDP_create_D_L(lat);
    d1 = QDP_create_D_L(lat);
    d2 = QDP_create_D_L(lat);
    // -------------------------------------------------
    QDP_D_eq_D_plus_D( s0,psi_up[0],psi_dw[0],all);
    QDP_D_eq_D_plus_D( s1,psi_up[1],psi_dw[1],all);
    QDP_D_eq_D_plus_D( s2,psi_up[2],psi_dw[2],all);
    QDP_D_eq_D_minus_D(d0,psi_up[0],psi_dw[0],all);
    QDP_D_eq_D_minus_D(d1,psi_up[1],psi_dw[1],all);
    QDP_D_eq_D_minus_D(d2,psi_up[2],psi_dw[2],all);
    // -------------------------------------------------
    // 1
    QDP_D_eq_D_plus_D(tempP,s1,s2,all);
    QDP_D_eq_M_times_sD(tempD1,gauge[0],tempP,neighbor[0],QDP_forward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[0],tempP,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[0],QDP_backward,all);
    QDP_D_eq_D_minus_D(tempP,tempD1,tempD2,all); //D_0(S_1+S_2)
    //---
    QDP_D_eq_M_times_sD(tempD1,gauge[1],d0,neighbor[1],QDP_forward,all);
    QDP_D_eq_M_times_sD(tempD2,gauge[2],d0,neighbor[2],QDP_forward,all);
    QDP_D_eq_D_plus_D(tempM,tempD1,tempD2,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[1],d0,all);
    QDP_D_eq_sD(tempD1,tempD3,neighbor[1],QDP_backward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[2],d0,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[2],QDP_backward,all);
    QDP_D_peq_D(tempD1,tempD2,all);
    QDP_D_peq_D(tempM,tempD1,all);
    QDP_D_peq_D(tempP,tempM,all);
    QDP_D_eq_r_times_D(tempP,&c2half,tempP,all);
    QDP_D_eq_gamma_times_D(tempM,tempP,gidx[0],all);
    QDP_D_peq_D(result,tempM,all);
    // 2
    QDP_D_eq_D_plus_D(tempP,s0,s2,all);
    QDP_D_eq_M_times_sD(tempD1,gauge[1],tempP,neighbor[1],QDP_forward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[1],tempP,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[1],QDP_backward,all);
    QDP_D_eq_D_minus_D(tempP,tempD1,tempD2,all); //D_1(S_0+S_2)
    //---
    QDP_D_eq_M_times_sD(tempD1,gauge[0],d1,neighbor[0],QDP_forward,all);
    QDP_D_eq_M_times_sD(tempD2,gauge[2],d1,neighbor[2],QDP_forward,all);
    QDP_D_eq_D_plus_D(tempM,tempD1,tempD2,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[0],d1,all);
    QDP_D_eq_sD(tempD1,tempD3,neighbor[0],QDP_backward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[2],d1,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[2],QDP_backward,all);
    QDP_D_peq_D(tempD1,tempD2,all);
    QDP_D_peq_D(tempM,tempD1,all);
    QDP_D_peq_D(tempP,tempM,all);
    QDP_D_eq_r_times_D(tempP,&c2half,tempP,all);
    QDP_D_eq_gamma_times_D(tempM,tempP,gidx[1],all);
    QDP_D_peq_D(result,tempM,all);
    // 3
    QDP_D_eq_D_plus_D(tempP,s0,s1,all);
    QDP_D_eq_M_times_sD(tempD1,gauge[2],tempP,neighbor[2],QDP_forward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[2],tempP,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[2],QDP_backward,all);
    QDP_D_eq_D_minus_D(tempP,tempD1,tempD2,all); //D_2(S_0+S_1)
    //---
    QDP_D_eq_M_times_sD(tempD1,gauge[0],d2,neighbor[0],QDP_forward,all);
    QDP_D_eq_M_times_sD(tempD2,gauge[1],d2,neighbor[1],QDP_forward,all);
    QDP_D_eq_D_plus_D(tempM,tempD1,tempD2,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[0],d2,all);
    QDP_D_eq_sD(tempD1,tempD3,neighbor[0],QDP_backward,all);
    QDP_D_eq_Ma_times_D(tempD3,gauge[1],d2,all);
    QDP_D_eq_sD(tempD2,tempD3,neighbor[1],QDP_backward,all);
    QDP_D_peq_D(tempD1,tempD2,all);
    QDP_D_peq_D(tempM,tempD1,all);
    QDP_D_peq_D(tempP,tempM,all);
    QDP_D_eq_r_times_D(tempP,&c2half,tempP,all);
    QDP_D_eq_gamma_times_D(tempM,tempP,gidx[2],all);
    QDP_D_peq_D(result,tempM,all);
    

    // ---------------------------------------------------
    QDP_destroy_D(s0);
    QDP_destroy_D(s1);
    QDP_destroy_D(s2);
    QDP_destroy_D(d0);
    QDP_destroy_D(d1);
    QDP_destroy_D(d2);
    
  };
  
  
/* ------------------------------------------------ */
/* ------- end of c_2 non-diagonal terms ---------- */
/* ------------------------------------------------ */
 
/* ------------------------------------------------ */
/* ------- start of c_E and c_EE terms --- -------- */
/* ------------------------------------------------ */
 
 if((c_E!=0)||(c_EE!=0)){
   
   QLA_Real cEhalf   = (0.5*c_E*zeta)*tpf*tpf*tpf;
   QLA_Real cE2half  = (0.5*c_EE)*tpf*tpf*tpf*tpf;
   
   QDP_ColorMatrix  *e0,*e1,*e2;
   
   e0 = QDP_create_M_L(lat);
   e1 = QDP_create_M_L(lat);
   e2 = QDP_create_M_L(lat);

   f_mu_nu(e0,1.0,gauge,3,0); // E_0 = F_{30} :: Alpha[0] = -Qamma[9]
   f_mu_nu(e1,1.0,gauge,3,1); // E_1 = F_{31} :: Alpha[1] = -Qamma[10]
   f_mu_nu(e2,1.0,gauge,3,2); // E_2 = F_{32} :: Alpha[2] = -Qamma[12]
   // tempD4 = [U_3(x).Psi(x+3)-Ud_3(x-3).Psi(x-3)]
   QDP_D_eq_D_minus_D(tempD4,psi_up[3],psi_dw[3],all);
   // --------------------------------------------------- 
   /* C_E TERM */
   QDP_D_eq_M_times_D(tempD1,e0,source,all); //E_0 Psi(x)
   QDP_D_eq_gamma_times_D(tempP,tempD1,aidx[0],all);// +Qop[9]E_0Psi(x)  
   QDP_D_eq_M_times_D(tempD2,e1,source,all); //E_1 Psi(x)
   QDP_D_eq_gamma_times_D(tempM,tempD2,aidx[1],all);// +Qop[10]E_1Psi(x)
   QDP_D_peq_D(tempP,tempM,all); // (Qop[9]E_0+Qop[10]E_1)Psi(x)
   QDP_D_eq_M_times_D(tempD3,e2,source,all); // E_2 Psi(x)
   QDP_D_eq_gamma_times_D(tempM,tempD3,aidx[2],all); // +Qop[12]E_2Psi(x)
   QDP_D_peq_D(tempP,tempM,all); // (Qop[9]E_0+Qop[10]E_1+Qop[12]E_2)Psi(x)
   // \alpha.E.Psi(x) = -tempP
   QDP_D_eq_r_times_D(tempM,&cEhalf,tempP,all);
   QDP_D_peq_D(result,tempM,all); // plus becaue of - sign from the op //
   // --------------------------------------------------
   /* c_EE TERM */
   if(c_EE!=0){
     QDP_D_eq_M_times_sD(tempD1,gauge[3],tempP,neighbor[3],QDP_forward,all);
     QDP_D_eq_Ma_times_D(tempD3,gauge[3],tempP,all);
     QDP_D_eq_sD(tempD2,tempD3,neighbor[3],QDP_backward,all);
     QDP_D_meq_D(tempD1,tempD2,all);
     QDP_D_eq_gamma_times_D(tempM,tempD1,8,all);//-Gam[4]D_4\alphaE
     // --------------------------------------
     QDP_D_eq_gamma_times_D(tempD2,tempD4,gidx[0],all);
     QDP_D_eq_M_times_D(tempD3,e0,tempD2,all);
     QDP_D_peq_D(tempM,tempD3,all);
     //
     QDP_D_eq_gamma_times_D(tempD2,tempD4,gidx[1],all);
     QDP_D_eq_M_times_D(tempD3,e1,tempD2,all);
     QDP_D_peq_D(tempM,tempD3,all);
     //
     QDP_D_eq_gamma_times_D(tempD2,tempD4,gidx[2],all);
     QDP_D_eq_M_times_D(tempD3,e2,tempD2,all);
     QDP_D_peq_D(tempM,tempD3,all);
     // ----
     QDP_D_eq_r_times_D(tempM,&cE2half,tempM,all);
     QDP_D_meq_D(result,tempM,all);
     
   };
   QDP_destroy_M(e0);
   QDP_destroy_M(e1);
   QDP_destroy_M(e2);
  
 };
 /* ------------------------------------------------ */
 /* --------- end of c_E and c_EE terms ------------ */
 /* ------------------------------------------------ */


 
 
 /* ------------------------------------------------ */
 /* ------- start of c_B+c_3+c+5 ------ ------------ */
 /* ------------------------------------------------ */ 
 
 f_mu_nu(f12,1.0,gauge,1,2);
 f_mu_nu(f02,1.0,gauge,0,2); //*****
 f_mu_nu(f01,1.0,gauge,0,1); 

 // LOOP NO : 8
 //******if((c_B!=0)||(c_3!=0)){
   
   QLA_Real cBhalf = (0.5*c_B*zeta+8.0*c_5)*tpf*tpf*tpf;
   QLA_Real c3half = (0.5*c_3)*tpf*tpf*tpf*tpf;
   
   //printf("cBhalf = %e\n",cBhalf);
   //printf("c3half = %e\n",c3half);

   
   QDP_DiracFermion *b0psi;
   QDP_DiracFermion *b1psi;
   QDP_DiracFermion *b2psi;
   QDP_DiracFermion *sigB;
   b0psi = QDP_create_D_L(lat);
   b1psi = QDP_create_D_L(lat);
   b2psi = QDP_create_D_L(lat);
   sigB  = QDP_create_D_L(lat);

   // tempG1 = F_{12} = B_0
   // b0psi  = B_0 Psi(x)
   // tempD1 = Gamma[6].B_0.Psi(x) = (+)iSigma[0].B_0.Psi(x)
   /* f_mu_nu(f12,1.0,gauge,1,2); */
   QDP_D_eq_M_times_D(b0psi,f12,source,all);
   QDP_D_eq_gamma_times_D(tempD1,b0psi,6,all);

   // tempG2 = F_{20} = B_1
   // b1psi  = B_1 Psi(x)
   // tempD2 = Gamma[5].B_1.Psi(x) = (-)iSigma[1].B_1.Psi(x)
   /* f_mu_nu(f02,1.0,gauge,0,2); */ //*****
   QDP_D_eq_M_times_D(b1psi,f02,source,all);
   QDP_D_eq_gamma_times_D(tempD2,b1psi,5,all);
   
   // tempG3 = F_{01} = B_2
   // b2psi  = B_2 Psi(x)
   // tempD3 = Gamma[3].B_2.Psi(x) = (+)iSigma[2].B_2.Psi(x)
   /* f_mu_nu(f01,1.0,gauge,0,1); */
   QDP_D_eq_M_times_D(b2psi,f01,source,all);
   QDP_D_eq_gamma_times_D(tempD3,b2psi,3,all);
   
   // tempM = tempD1 - tempD2 Because : iSigma[1] = -Qamma[5]
   // sigB  = tempM  + tempD3
   QDP_D_eq_D_plus_D(tempM,tempD1,tempD2,all); //**** minus->plus
   QDP_D_eq_D_plus_D(sigB,tempM,tempD3,all);
   //
   QDP_D_eq_r_times_D(tempP,&cBhalf,sigB,all);
   // minus sign comes from the -c_B/2\Sigma.B term 
   QDP_D_meq_D(result,tempP,all);
    
   QDP_destroy_D(b0psi);
   QDP_destroy_D(b1psi);
   QDP_destroy_D(b2psi);
  
   // Now the c_3 part 
   //*****if(c3half!=0){
     // iSigma.B Gamma.D Psi(x)
     // 1 : tempG1 = F_{12}
     QDP_D_eq_M_times_D(tempD1,f12,gamD,all);
     QDP_D_eq_gamma_times_D(tempP,tempD1,6,all);
     // 2 : tempG2 = F_{20}
     QDP_D_eq_M_times_D(tempD1,f02,gamD,all);
     QDP_D_eq_gamma_times_D(tempD2,tempD1,5,all);
     QDP_D_peq_D(tempP,tempD2,all); //meq->peq
     //3
     QDP_D_eq_M_times_D(tempD1,f01,gamD,all);
     QDP_D_eq_gamma_times_D(tempD2,tempD1,3,all);
     QDP_D_peq_D(tempP,tempD2,all);
     //
     // gamma.D Sigma.B Psi(x)
     for(dir=0;dir<3;dir++){
       QDP_D_eq_M_times_sD(tempD1,gauge[dir],sigB,neighbor[dir],QDP_forward,all);
       QDP_D_eq_Ma_times_D(tempD3,gauge[dir],sigB,all);
       QDP_D_eq_sD(tempD2,tempD3,neighbor[dir],QDP_backward,all);
       QDP_D_meq_D(tempD1,tempD2,all);
       QDP_D_eq_gamma_times_D(tempD2,tempD1,gidx[dir],all);
       QDP_D_peq_D(tempP,tempD2,all);
     };
     //
     QDP_D_eq_r_times_D(tempP,&c3half,tempP,all);
     QDP_D_peq_D(result,tempP,all);
     //
     //*****}; // end of c3half
   
     //*****}; // end of LOOP 8
 
 
 // I no longer need tempP and tempM
 // empty memory
 //QDP_destroy_D(tempP);
     QDP_destroy_D(sigB);
     QDP_destroy_D(tempM);
     QDP_destroy_D(gamD);

 
 // Below is the new c_5 term where the three and five
 // link terms are seperated and multiplied with the 
 // appropriate power of tadpole factor. Keeping it separete for now.
 // Factor of 2 coming from \Delta is observed in c_B calculation above.
 // ::::::
 QLA_Real oneeight = +0.125;
 QLA_Real two      = +2.0;
 
 QLA_Real c5u3     = c_5*tpf*tpf;
 QLA_Real c5u5     = c_5*tpf*tpf*tpf*tpf;

 QDP_ColorMatrix *mat;
 mat = QDP_create_M_L(lat);
 
 //1
#if 1
 //f_mu_nu(f12,1,gauge,1,2);
 
 // (I) i\Sigma[0][ B_0 (u1+d1+u2+d2) + (u1+d1+u2+d2) B_0 ] psi(x)
 // i\Sigma[0] = QOP_Gamma[6] and B_0 = F_{12} 
 // ##############################################################
 // ##############################################################
 // ##############################################################
 // LINKS-SEPARATED VERSION 
 // ###########################################################
 // M = -U_2(x)U_1(x+2)Ud_2(x+1)+Ud_2(x-2)U_1(x-2)U_2(x+1-2)
 QDP_M_eq_sM(tempG1,gauge[1],neighbor[2],QDP_forward,all);
 QDP_M_eq_sM(tempG2,gauge[2],neighbor[1],QDP_forward,all);
 QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
 QDP_M_eq_M_times_M(tempG1,gauge[2],tempG3,all);
 //---
 QDP_M_eq_Ma_times_M(tempG2,gauge[2],gauge[1],all);
 QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[2],neighbor[1],QDP_forward,all);
 QDP_M_eq_sM(mat,tempG3,neighbor[2],QDP_backward,all);
 //---
 QDP_M_meq_M(mat,tempG1,all);
 QDP_M_eq_r_times_M(mat,&oneeight,mat,all); // times 2
 // ###########################################################
 // u1
 QDP_M_eq_M_times_Ma(tempG1,mat,gauge[1],all);
 QDP_M_eq_M_minus_M(tempG2,f12,tempG1,all);
 // R :QDP_D_eq_M_times_sD(tempD1,gauge[1],source,neighbor[1],QDP_forward,all);
 // W : tempD1 --> psi_up[1]
 QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[1],all); // YES
 // d1
 QDP_M_eq_Ma_times_M(tempG1,mat,gauge[1],all);
 QDP_M_eq_sM(tempG2,tempG1,neighbor[1],QDP_backward,all);
 QDP_M_eq_M_plus_M(tempG1,f12,tempG2,all);
 
 //R:QDP_D_eq_Ma_times_D(tempD1,gauge[1],source,all);
 //R:QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[1],QDP_backward,all);
 QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[1],all); // YES
 QDP_D_peq_D(tempD2,tempD3,all);
  //u1
 QDP_M_eq_Ma_times_M(tempG1,gauge[1],mat,all);
 QDP_M_eq_sM(tempG2,tempG1,neighbor[1],QDP_backward,all);
 QDP_M_eq_M_minus_M(tempG3,f12,tempG2,all);
 QDP_M_eq_M_times_sM(tempG1,gauge[1],tempG3,neighbor[1],QDP_forward,all);
 QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[1],QDP_forward,all);
 QDP_D_peq_D(tempD2,tempD3,all);
  //d1
 QDP_M_eq_M_times_Ma(tempG1,gauge[1],mat,all);
 QDP_M_eq_M_plus_M(tempG2,f12,tempG1,all);
 QDP_M_eq_Ma_times_M(tempG1,gauge[1],tempG2,all);
 QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
 QDP_D_eq_sD(tempD3,tempD1,neighbor[1],QDP_backward,all);
 //
 QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
  //
 QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
 QDP_D_eq_gamma_times_D(tempD4,tempD2,6,all);
 QDP_D_peq_D(result,tempD4,all); // all 5-links added to result
 //
 // u1 M(x)Psi(x+1)
 QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[1],QDP_forward,all);
 //d1
 QDP_D_eq_Ma_times_D(tempD1,mat,source,all); // Md(x)Psi(x)
 // Md(x-1)Psi(x-1)
 QDP_D_eq_sD(tempD3,tempD1,neighbor[1],QDP_backward,all);
 QDP_D_meq_D(tempD4,tempD3,all);
 QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
 //
 QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
 QDP_D_eq_gamma_times_D(tempD3,tempD4,6,all);
 QDP_D_peq_D(result,tempD3,all); // all 3-links added to result
 
 // ###########################################################
 QDP_M_eq_sM(tempG1,gauge[2],neighbor[1],QDP_forward,all);
 QDP_M_eq_sM(tempG2,gauge[1],neighbor[2],QDP_forward,all);
 QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
 QDP_M_eq_M_times_M(mat,gauge[1],tempG3,all);
 //---
 QDP_M_eq_Ma_times_M(tempG2,gauge[1],gauge[2],all);
 QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[1],neighbor[2],QDP_forward,all);
 QDP_M_eq_sM(tempG1,tempG3,neighbor[1],QDP_backward,all);
 //---
 QDP_M_meq_M(mat,tempG1,all);
 QDP_M_eq_r_times_M(mat,&oneeight,mat,all);
 // ###########################################################
 // u2
 QDP_M_eq_M_times_Ma(tempG1,mat,gauge[2],all);
 QDP_M_eq_M_minus_M(tempG2,f12,tempG1,all);
 // R : QDP_D_eq_M_times_sD(tempD1,gauge[2],source,neighbor[2],QDP_forward,all);
 // W : tempD1 = psi_up[2]-----------V
 QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[2],all); // YES
 // d2
 QDP_M_eq_Ma_times_M(tempG1,mat,gauge[2],all);
 QDP_M_eq_sM(tempG2,tempG1,neighbor[2],QDP_backward,all);
 QDP_M_eq_M_plus_M(tempG1,f12,tempG2,all);

 //R: QDP_D_eq_Ma_times_D(tempD1,gauge[2],source,all);
 //R: QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[2],QDP_backward,all);
 // W : V
 QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[2],all); // YES
 QDP_D_peq_D(tempD2,tempD3,all);
 //u2
 QDP_M_eq_Ma_times_M(tempG1,gauge[2],mat,all);
 QDP_M_eq_sM(tempG2,tempG1,neighbor[2],QDP_backward,all);
 QDP_M_eq_M_minus_M(tempG3,f12,tempG2,all);
 QDP_M_eq_M_times_sM(tempG1,gauge[2],tempG3,neighbor[2],QDP_forward,all);
 QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[2],QDP_forward,all);
 QDP_D_peq_D(tempD2,tempD3,all);
 //d2
 QDP_M_eq_M_times_Ma(tempG1,gauge[2],mat,all);
 QDP_M_eq_M_plus_M(tempG2,f12,tempG1,all);
 QDP_M_eq_Ma_times_M(tempG1,gauge[2],tempG2,all);
 QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
 QDP_D_eq_sD(tempD3,tempD1,neighbor[2],QDP_backward,all);
 //
 QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
 QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
 QDP_D_eq_gamma_times_D(tempD3,tempD2,6,all);
 QDP_D_peq_D(result,tempD3,all); // all 5-links added to result
 //
 // u1
 QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[2],QDP_forward,all);
 //d1
 QDP_D_eq_Ma_times_D(tempD1,mat,source,all);
 QDP_D_eq_sD(tempD3,tempD1,neighbor[2],QDP_backward,all);
 QDP_D_meq_D(tempD4,tempD3,all);
 QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
 //
 QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
 QDP_D_eq_gamma_times_D(tempD2,tempD4,6,all);
 QDP_D_peq_D(result,tempD2,all); // all 3-links added to result

 // I no longer need F12
 QDP_destroy_M(f12);
#endif

 //2
#if 1
 // ##############################################################
 // ##############################################################
 // ##############################################################
 // (I) i\Sigma[1][ B_1 (u0+d0+u2+d2) + (u0+d0+u2+d2) B_0 ] psi(x)
 // i\Sigma[1] = -QOP_Gamma[5] and B_1 = -F_{02} 
 // ##############################################################
 // ##############################################################
 // ##############################################################
 // LINKS-SEPARATED VERSION 
 // ###########################################################
 // M = Ud_2(x-2)U_0(x-0)U_2(x+0-2)-U_2(x)U_0(x+2)Ud_2(x+0)

 
 //f_mu_nu(f02,1,gauge,0,2);
 
 QDP_M_eq_sM(tempG1,gauge[0],neighbor[2],QDP_forward,all);
 QDP_M_eq_sM(tempG2,gauge[2],neighbor[0],QDP_forward,all);
 QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
 QDP_M_eq_M_times_M(tempG1,gauge[2],tempG3,all);
 //---
 QDP_M_eq_Ma_times_M(tempG2,gauge[2],gauge[0],all);
 QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[2],neighbor[0],QDP_forward,all);
 QDP_M_eq_sM(mat,tempG3,neighbor[2],QDP_backward,all);
 //---
 QDP_M_meq_M(mat,tempG1,all);
 QDP_M_eq_r_times_M(mat,&oneeight,mat,all);
 // ###########################################################
 // u1
 QDP_M_eq_M_times_Ma(tempG1,mat,gauge[0],all);
 QDP_M_eq_M_minus_M(tempG2,f02,tempG1,all);
 //R: QDP_D_eq_M_times_sD(tempD1,gauge[0],source,neighbor[0],QDP_forward,all);
 //W: tempD1 = psi_up[0]-----------V
 QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[0],all); //YES
  // d1
  QDP_M_eq_Ma_times_M(tempG1,mat,gauge[0],all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[0],QDP_backward,all);
  QDP_M_eq_M_plus_M(tempG1,f02,tempG2,all);
  
  //R:QDP_D_eq_Ma_times_D(tempD1,gauge[0],source,all);
  //R:QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[0],QDP_backward,all);
  //W: tempD1 = psi_dw[0]---^
  QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[0],all); //YES
  QDP_D_peq_D(tempD2,tempD3,all);
  //u1
  QDP_M_eq_Ma_times_M(tempG1,gauge[0],mat,all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[0],QDP_backward,all);
  QDP_M_eq_M_minus_M(tempG3,f02,tempG2,all);
  QDP_M_eq_M_times_sM(tempG1,gauge[0],tempG3,neighbor[0],QDP_forward,all);
  QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[0],QDP_forward,all);
  QDP_D_peq_D(tempD2,tempD3,all);
  //d1
  QDP_M_eq_M_times_Ma(tempG1,gauge[0],mat,all);
  QDP_M_eq_M_plus_M(tempG2,f02,tempG1,all);
  QDP_M_eq_Ma_times_M(tempG1,gauge[0],tempG2,all);
  QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[0],QDP_backward,all);
  //
  QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
  //
  QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
  QDP_D_eq_gamma_times_D(tempD4,tempD2,5,all);
  QDP_D_peq_D(result,tempD4,all); // 5-links added to result
  //
  // u1
  QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[0],QDP_forward,all);
  //d1
  QDP_D_eq_Ma_times_D(tempD1,mat,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[0],QDP_backward,all);
  QDP_D_meq_D(tempD4,tempD3,all);
  QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
  //
  QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
 
  QDP_D_eq_gamma_times_D(tempD3,tempD4,5,all);
  QDP_D_peq_D(result,tempD3,all);
  
  // ###########################################################
  // ###########################################################
  // ###########################################################
  
  QDP_M_eq_sM(tempG1,gauge[2],neighbor[0],QDP_forward,all);
  QDP_M_eq_sM(tempG2,gauge[0],neighbor[2],QDP_forward,all);
  QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
  QDP_M_eq_M_times_M(mat,gauge[0],tempG3,all);
  //---
  QDP_M_eq_Ma_times_M(tempG2,gauge[0],gauge[2],all);
  QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[0],neighbor[2],QDP_forward,all);
  QDP_M_eq_sM(tempG1,tempG3,neighbor[0],QDP_backward,all);
  //---
  QDP_M_meq_M(mat,tempG1,all);
  QDP_M_eq_r_times_M(mat,&oneeight,mat,all);
  // ###########################################################
  // u2
  QDP_M_eq_M_times_Ma(tempG1,mat,gauge[2],all);
  QDP_M_eq_M_minus_M(tempG2,f02,tempG1,all);
  //R:QDP_D_eq_M_times_sD(tempD1,gauge[2],source,neighbor[2],QDP_forward,all);
  // W: tempD1 = psi_up[2]-----------V
  QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[2],all); //YES
  // d2
  QDP_M_eq_Ma_times_M(tempG1,mat,gauge[2],all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[2],QDP_backward,all);
  QDP_M_eq_M_plus_M(tempG1,f02,tempG2,all);
  //R:QDP_D_eq_Ma_times_D(tempD1,gauge[2],source,all);
  //R:QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[2],QDP_backward,all);
  // W:
  QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[2],all); //YES
  QDP_D_peq_D(tempD2,tempD3,all);
  //u2
  QDP_M_eq_Ma_times_M(tempG1,gauge[2],mat,all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[2],QDP_backward,all);
  QDP_M_eq_M_minus_M(tempG3,f02,tempG2,all);
  QDP_M_eq_M_times_sM(tempG1,gauge[2],tempG3,neighbor[2],QDP_forward,all);
  QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[2],QDP_forward,all);
  QDP_D_peq_D(tempD2,tempD3,all);
  //d2
  QDP_M_eq_M_times_Ma(tempG1,gauge[2],mat,all);
  QDP_M_eq_M_plus_M(tempG2,f02,tempG1,all);
  QDP_M_eq_Ma_times_M(tempG1,gauge[2],tempG2,all);
  QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[2],QDP_backward,all);
  //
  QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
  QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
  QDP_D_eq_gamma_times_D(tempD3,tempD2,5,all);
  QDP_D_peq_D(result,tempD3,all); // 5 links added to the result 
  //
  // u1
  QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[2],QDP_forward,all);
  //d1
  QDP_D_eq_Ma_times_D(tempD1,mat,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[2],QDP_backward,all);
  QDP_D_meq_D(tempD4,tempD3,all);
  QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
  //
  QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
  QDP_D_eq_gamma_times_D(tempD2,tempD4,5,all);
  QDP_D_peq_D(result,tempD2,all);
  
  // I no longer need f_02
  QDP_destroy_M(f02);

#endif

  //3
#if 1

  // ##############################################################
  // ##############################################################
  // ##############################################################
  // (I) i\Sigma[2][ B_2 (u0+d0+u1+d1) + (u0+d0+u1+d1) B_2 ] psi(x)
  // i\Sigma[2] = +QOP_Gamma[3] and B_1 = +F_{01} 
  // ##############################################################
  // ##############################################################
  // ##############################################################
  // LINKS-SEPARATED VERSION 
  // ###########################################################
  // M = 

  //f_mu_nu(f01,1,gauge,0,1);

  QDP_M_eq_sM(tempG1,gauge[0],neighbor[1],QDP_forward,all);
  QDP_M_eq_sM(tempG2,gauge[1],neighbor[0],QDP_forward,all);
  QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
  QDP_M_eq_M_times_M(tempG1,gauge[1],tempG3,all);
  //---
  QDP_M_eq_Ma_times_M(tempG2,gauge[1],gauge[0],all);
  QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[1],neighbor[0],QDP_forward,all);
  QDP_M_eq_sM(mat,tempG3,neighbor[1],QDP_backward,all);
  //---
  QDP_M_meq_M(mat,tempG1,all);
  QDP_M_eq_r_times_M(mat,&oneeight,mat,all);
  // ###########################################################
  // u1
  QDP_M_eq_M_times_Ma(tempG1,mat,gauge[0],all);
  QDP_M_eq_M_minus_M(tempG2,f01,tempG1,all);
  //R:QDP_D_eq_M_times_sD(tempD1,gauge[0],source,neighbor[0],QDP_forward,all);
  //W:tempD1 = psi_up[0]------------V
  QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[0],all); //YES
  // d1
  QDP_M_eq_Ma_times_M(tempG1,mat,gauge[0],all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[0],QDP_backward,all);
  QDP_M_eq_M_plus_M(tempG1,f01,tempG2,all);
  //R:QDP_D_eq_Ma_times_D(tempD1,gauge[0],source,all);
  //R:QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[0],QDP_backward,all);
  //W: tempD1 = psi_dw[0];
  QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[0],all); //YES
 QDP_D_peq_D(tempD2,tempD3,all);
  //u1
  QDP_M_eq_Ma_times_M(tempG1,gauge[0],mat,all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[0],QDP_backward,all);
  QDP_M_eq_M_minus_M(tempG3,f01,tempG2,all);
  QDP_M_eq_M_times_sM(tempG1,gauge[0],tempG3,neighbor[0],QDP_forward,all);
  QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[0],QDP_forward,all);
  QDP_D_peq_D(tempD2,tempD3,all);
  //d1
  QDP_M_eq_M_times_Ma(tempG1,gauge[0],mat,all);
  QDP_M_eq_M_plus_M(tempG2,f01,tempG1,all);
  QDP_M_eq_Ma_times_M(tempG1,gauge[0],tempG2,all);
  QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[0],QDP_backward,all);
  //
  QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
  //
  QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
  QDP_D_eq_gamma_times_D(tempD4,tempD2,3,all);
  QDP_D_peq_D(result,tempD4,all); // 5-links added to the result
  //
  // u1
  QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[0],QDP_forward,all);
  //d1
  QDP_D_eq_Ma_times_D(tempD1,mat,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[0],QDP_backward,all);
  QDP_D_meq_D(tempD4,tempD3,all);
  QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
  //
  QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
  QDP_D_eq_gamma_times_D(tempD3,tempD4,3,all);
  QDP_D_peq_D(result,tempD3,all); // 3-links added to the result
  
  
  // ###########################################################
  QDP_M_eq_sM(tempG1,gauge[1],neighbor[0],QDP_forward,all);
  QDP_M_eq_sM(tempG2,gauge[0],neighbor[1],QDP_forward,all);
  QDP_M_eq_M_times_Ma(tempG3,tempG1,tempG2,all);
  QDP_M_eq_M_times_M(mat,gauge[0],tempG3,all);
  //---
  QDP_M_eq_Ma_times_M(tempG2,gauge[0],gauge[1],all);
  QDP_M_eq_M_times_sM(tempG3,tempG2,gauge[0],neighbor[1],QDP_forward,all);
  QDP_M_eq_sM(tempG1,tempG3,neighbor[0],QDP_backward,all);
  //---
  QDP_M_meq_M(mat,tempG1,all);
  QDP_M_eq_r_times_M(mat,&oneeight,mat,all);
  // ###########################################################
  // u2
  QDP_M_eq_M_times_Ma(tempG1,mat,gauge[1],all);
  QDP_M_eq_M_minus_M(tempG2,f01,tempG1,all);
  //R:QDP_D_eq_M_times_sD(tempD1,gauge[1],source,neighbor[1],QDP_forward,all);
  // W: tempD1 = psi_up[1]-----------V
  QDP_D_eq_M_times_D(tempD2,tempG2,psi_up[1],all); // YES
  // d2
  QDP_M_eq_Ma_times_M(tempG1,mat,gauge[1],all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[1],QDP_backward,all);
  QDP_M_eq_M_plus_M(tempG1,f01,tempG2,all);
  //R:QDP_D_eq_Ma_times_D(tempD1,gauge[1],source,all);
  //R:QDP_D_eq_M_times_sD(tempD3,tempG1,tempD1,neighbor[1],QDP_backward,all);
  //W : tempD1 =psi_dw[1]
  QDP_D_eq_M_times_D(tempD3,tempG1,psi_dw[1],all); // YES
  QDP_D_peq_D(tempD2,tempD3,all);
  //u2
  QDP_M_eq_Ma_times_M(tempG1,gauge[1],mat,all);
  QDP_M_eq_sM(tempG2,tempG1,neighbor[1],QDP_backward,all);
  QDP_M_eq_M_minus_M(tempG3,f01,tempG2,all);
  QDP_M_eq_M_times_sM(tempG1,gauge[1],tempG3,neighbor[1],QDP_forward,all);
  QDP_D_eq_M_times_sD(tempD3,tempG1,source,neighbor[1],QDP_forward,all);
  QDP_D_peq_D(tempD2,tempD3,all);
  //d2
  QDP_M_eq_M_times_Ma(tempG1,gauge[1],mat,all);
  QDP_M_eq_M_plus_M(tempG2,f01,tempG1,all);
  QDP_M_eq_Ma_times_M(tempG1,gauge[1],tempG2,all);
  QDP_D_eq_M_times_D(tempD1,tempG1,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[1],QDP_backward,all);
  //
  QDP_D_peq_D(tempD2,tempD3,all); // all 5-links
  QDP_D_eq_r_times_D(tempD2,&c5u5,tempD2,all);
  QDP_D_eq_gamma_times_D(tempD3,tempD2,3,all);
  QDP_D_peq_D(result,tempD3,all); // 5-links added
  //
  // u1
  QDP_D_eq_M_times_sD(tempD4,mat,source,neighbor[1],QDP_forward,all);
  //d1
  QDP_D_eq_Ma_times_D(tempD1,mat,source,all);
  QDP_D_eq_sD(tempD3,tempD1,neighbor[1],QDP_backward,all);
  QDP_D_meq_D(tempD4,tempD3,all);
  QDP_D_eq_r_times_D(tempD4,&two,tempD4,all);
  //
  QDP_D_eq_r_times_D(tempD4,&c5u3,tempD4,all);
  QDP_D_eq_gamma_times_D(tempD2,tempD4,3,all);
  QDP_D_peq_D(result,tempD2,all);

  QDP_destroy_M(f01);
#endif

#if 0
  //static QLA_DiracFermion *qladf;
  //static int st,nc,ns;
  //static QLA_Real c1,c2;
  qladf = QDP_expose_D(result);
  printf("Exposing Dirac Fermion :: result \n");  
  for(st=1;st<2;st++){
    for(nc=0;nc<3;nc++){
      for(ns=0;ns<1;ns++){
	c1=QLA_real(QLA_elem_D(qladf[st],nc,ns));
	c2=QLA_imag(QLA_elem_D(qladf[st],nc,ns));
	printf("site : %d  color : %d  spin  : %d = (%e,%e)\n",st,nc,ns,c1,c2);
      };
    };
    printf("***********************************\n");
  };
  QDP_reset_D(result);
#endif
  // ------------------------------------------------
  //  printf("Tracing my steps..... *** D-slash completed\n");

  
  


 // DESTROYING FIELDS ----------------------------------
 for(dir=0;dir<4;dir++){
 QDP_destroy_D(psi_up[dir]);
 QDP_destroy_D(psi_dw[dir]);
 QDP_destroy_M(gauge[dir]);
 };
 QDP_destroy_M(tempG1);
 QDP_destroy_M(tempG2);
 QDP_destroy_M(tempG3);
 QDP_destroy_D(tempD1);
 QDP_destroy_D(tempD2);
 QDP_destroy_D(tempD3);
 QDP_destroy_D(tempD4);
 QDP_destroy_M(mat);
 
 QDP_destroy_D(tempP);
};

#endif /* QOP_Colors == 3 */
