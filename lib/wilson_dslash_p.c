/* adapted from MILC version 6 */

/* conventions for raw packed clover term

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
c(0,1, 0,0)
c(1,0, 0,0)
c(1,1, 0,0)
c(2,0, 0,0)
c(2,1, 0,0)
c(1,0, 0,1)
c(1,1, 0,1)
c(2,0, 0,1)
c(2,1, 0,1)
c(1,1, 1,0)
c(2,0, 1,0)
c(2,1, 1,0)
c(2,0, 1,1)
c(2,1, 1,1)
c(2,1, 2,0)

and likewise for the second block by adding 2 to all is,js values above.
*/

#include <string.h>
#include <math.h>
#include <qop_internal.h>

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

#define CLOV_REALS (2*6*6) // 2 packed 6x6 Hermitian matrices
#define CLOV_SIZE (CLOV_REALS*sizeof(REAL)) 

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

/* ---------------------------------------------------------------- */
/* This part is added by Bugra :::::------------------------------- */
/* F_mu_nu calculates the $F_{\mu\nu}$ Field Strength Tensor  ----- */
/* --------------------------------------------------------------   */
static void
f_mu_nu(QDP_ColorMatrix *fmn, QLA_Real scale, QDP_ColorMatrix *link[], int mu, int nu)
{
  int order_flag;

  QDP_ColorMatrix *temp1,*temp2,*temp3,*temp4,*tmat4;
  QDP_ColorMatrix *pqt0,*pqt1,*pqt2,*pqt3;
  QDP_ColorMatrix *pqt4;

  temp1 = QDP_create_M();
  temp2 = QDP_create_M();
  temp3 = QDP_create_M();
  temp4 = QDP_create_M();
  tmat4 = QDP_create_M();
  pqt0  = QDP_create_M();
  pqt1  = QDP_create_M();
  pqt2  = QDP_create_M();
  pqt3  = QDP_create_M();
  pqt4  = QDP_create_M();

  /* chech if mu <nu */
  if(mu>nu){
    int i = mu;
    mu = nu;
    nu = i;
    order_flag=1;
  }
  else{
    order_flag=0;
  }

  /* Get pqt0 = U_nu(x+mu) : U_nu(x) from mu direction    */
  /* Get pqt1 = U_mu(x+nu) : U_mu(x) from nu direction    */
  QDP_M_eq_sM(pqt0,link[nu],QDP_neighbor[mu],QDP_forward,QDP_all);
  QDP_M_eq_sM(pqt1,link[mu],QDP_neighbor[nu],QDP_forward,QDP_all);

  /* creating a corner : temp1 = U^d_nu(x) U_mu(x) */
  QDP_M_eq_Ma_times_M(temp1,link[nu],link[mu],QDP_all);

  /* ------------------------------------------------- */
  /* creating a corner : fmn = [pqt0]*[pqt1^d]         */
  /*                     fmn = U_nu(x+mu) U^d_mu(x+nu) */
  QDP_M_eq_M_times_Ma(fmn,pqt0,pqt1,QDP_all);

  /* Create the following loops */
  /* U^d_nu(x) U_mu(x) U_nu(x+mu) U^d_mu(x+nu)  */
  QDP_M_eq_M_times_M(temp2,temp1,fmn,QDP_all);

  /* U_nu(x+mu) U^d_mu(x+nu) U^d_nu(x) U_mu(x) */
  QDP_M_eq_M_times_M(temp3,fmn,temp1,QDP_all);

  /* Creating the following +mu -nu plaquette */
  QDP_M_eq_sM(pqt2,temp2,QDP_neighbor[nu],QDP_backward,QDP_all);

  /* Creating the following -mu +nu plaquette */
  QDP_M_eq_sM(pqt3,temp3,QDP_neighbor[mu],QDP_backward,QDP_all);
  //QDP_discard_M(temp3); /* data in temp3 is no longer needed */

  /* creating +mu +nu plaquette and pit it in fmn      */
  /* tmat4 = U_mu(x) U_nu(x+mu) U^d_mu(x+nu)           */
  /* fmn   = U_mu(x) U_nu(x+mu) U^d_mu(x+nu) U^d_nu(x) */ 
  QDP_M_eq_M_times_M(tmat4,link[mu],fmn,QDP_all);
  QDP_M_eq_M_times_Ma(fmn,tmat4,link[nu],QDP_all);

  /* What is left is +mu -nu plaquette and adding them up */
  /* Right hand side of the clover field */
  QDP_M_eq_M_plus_M(fmn,pqt2,fmn,QDP_all);

  /* tmat4 = [(pqt[1]^d)   ] * [(tempmat1 )            ]  */
  /*       = [U^d_mu(x+nu) ] * [U^dagger_nu(x) *U_mu(x)]  */
  /* temp2 = [tmat4]*[pqt0 ]                            */
  /*       = [U^d_mu(x+nu) U^d_nu(x) U_mu(x)*[U_nu(x+mu)] */
  QDP_M_eq_Ma_times_M(tmat4,pqt1,temp1,QDP_all);
  QDP_M_eq_M_times_M(temp2,tmat4,pqt0,QDP_all);
  /* temp1 is a result of a shift and won't be needed   */
  QDP_discard_M(temp1);

  /* temp2 is now a plaquette -mu -nu and must be gathered */
  /* with displacement -mu-nu */
  /* pqt4 = U^d_mu(x-nu)U^d(x-mu-nu)U_mu(x-mu-nu)U_nu(x-nu) */
  QDP_M_eq_sM(temp4,temp2,QDP_neighbor[mu],QDP_backward,QDP_all);
  QDP_M_eq_sM(pqt4,temp4,QDP_neighbor[nu],QDP_backward,QDP_all);
  QDP_discard_M(temp2);

  /* Now gather -mu +nu  plaquette and add to fmn    */
  /* f_mn was the right hand side of the clover term */
  /* I add the third plaquette : f_mn = f_mn+pqt3    */
  /* U_nu(x) U^d(x-mu+nu) U^d_nu(x-mu) U_mu(x-mu)    */
  QDP_M_eq_M_plus_M(fmn,fmn,pqt3,QDP_all);

  /* finally add the last plaquette        */
  /* fmn = fmn+ pqt4                       */
  /* pqt4 is the last plaquette missing    */
  /* This completes the 4-plaquette        */
  QDP_M_eq_M_plus_M(fmn,fmn,pqt4,QDP_all);

  /* F_munu  is now 1/8 of f_mn-f_mn^dagger */
  /* QDP_T_eqop_Ta(Type *r, Type *a, subset) : r=conjugate(a) */
  /* tmat4  =Hermitian(fmn) */
  QDP_M_eq_Ma(tmat4,fmn,QDP_all);

  if(order_flag ==0){
    QDP_M_eq_M_minus_M(tmat4,fmn,tmat4,QDP_all);
  }
  else {
    QDP_M_eq_M_minus_M(tmat4,tmat4,fmn,QDP_all);
  }

  /* F_mn = 1/8 *tmat4 */
  /* QDP_T_eq_r_times_T(Type *r, QLA_real *a, Type *b, subset); */
  scale *= 0.125;
  QDP_M_eq_r_times_M(fmn, &scale, tmat4, QDP_all);

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
}

/* -------------------------------------------------------- */
/* -- Calculation of the clover term :: added by Bugra ---- */
/* -------------------------------------------------------- */
static void
get_clov(QLA_Real *clov, QDP_ColorMatrix *link[], QLA_Real cs, QLA_Real ct)
{
  /* Fist I create A[0],A[1],B[0] and B[1] matrices */
  QDP_ColorMatrix *a[2],*b[2];
  QDP_ColorMatrix *f_mn;
  QDP_ColorMatrix *xm,*ym;

  f_mn  = QDP_create_M();
  a[0]  = QDP_create_M();
  a[1]  = QDP_create_M();
  xm = QDP_create_M();
  ym = QDP_create_M();

  f_mu_nu(f_mn, cs, link, 0, 1);         /* f_mn = cs*F_{01}        */
  QDP_M_eq_M(xm, f_mn, QDP_all);         /* xm = cs*F_{01}        */
  QDP_M_eq_M(ym, f_mn, QDP_all);         /* ym = cs*F_{01}        */

  f_mu_nu(f_mn, ct, link, 2, 3);         /* f_mn = ct*F_{23}        */ 
  QDP_M_meq_M(xm, f_mn, QDP_all);        /* xm = cs*F_{01}-ct*F_{23} */
  QDP_M_peq_M(ym, f_mn, QDP_all);        /* ym = cs*F_{01}+ct*F_{23} */

  // PS: One cannot QDP_M_eq_i_M(a,b,subset) with a=b 
  QDP_M_eq_i_M(a[0], xm, QDP_all);       /* a[0] = i(cs*F_{01}-ct*F_{23}) */ 
  QDP_M_eq_i_M(a[1], ym, QDP_all);       /* a[1] = i(cs*F_{01}+ct*F_{23}) */

  /* PART 1 and PART 2 */
  QLA_ColorMatrix *A[2];
  A[0] = QDP_expose_M(a[0]);
  A[1] = QDP_expose_M(a[1]);

  for(int i=0; i<QDP_sites_on_node; i++) {
    for(int j=0; j<3; j++) {
      /* diagoal elements numbered from 00 to 05 for the matrix X */
      /* c(0,0, 0,0) = clov[0]  c(1,1, 1,1) = clov[3]             */
      /* c(0,1, 0,1) = clov[1]  c(2,0, 2,0) = clov[4]             */
      /* c(1,0, 1,0) = clov[2]  c(2,1, 2,1) = clov[5]             */
      //clov[(72*i)+(j   )] = QLA_real(QLA_elem_M(A[0][i],j,j));
      //clov[(72*i)+(j+3 )] = -QLA_real(QLA_elem_M(A[0][i],j,j));
      clov[(72*i)+(2*j  )] = QLA_real(QLA_elem_M(A[0][i],j,j));
      clov[(72*i)+(2*j+1)] = -QLA_real(QLA_elem_M(A[0][i],j,j));
      /* diagoal elements numbered from 36 to 41 for the matrix Y */
      /* c(0,2, 0,2) = clov[36]  c(1,3, 1,3) = clov[39]           */
      /* c(0,3, 0,3) = clov[37]  c(2,2, 2,2) = clov[40]           */
      /* c(1,2, 1,2) = clov[38]  c(2,3, 2,3) = clov[41]           */
      clov[(72*i)+(2*j+36)] = QLA_real(QLA_elem_M(A[1][i],j,j));
      clov[(72*i)+(2*j+37)] = -QLA_real(QLA_elem_M(A[1][i],j,j));
      /* -------------------------------------------------------- */
      /* 12 real number are assigned to the array clov            */

      /* Below the triangular which have A[0] and A[1] only       */
      for(int k=0; k<j; k++) { 
	/* c(1,0, 0,0) = clov[ 8]+i*clov[ 9] */
	/* c(2,0, 0,0) = clov[12]+i*clov[13] */
	/* c(2,0, 1,0) = clov[26]+i*clov[27] */
	int jk = 4 + (4*j) + (14*k);
	clov[(72*i)+(jk  ) ] = QLA_real(QLA_elem_M(A[0][i],j,k));
	clov[(72*i)+(jk+1) ] = QLA_imag(QLA_elem_M(A[0][i],j,k));
	clov[(72*i)+(jk+36)] = QLA_real(QLA_elem_M(A[1][i],j,k));
	clov[(72*i)+(jk+37)] = QLA_imag(QLA_elem_M(A[1][i],j,k));

	/* c(1,1, 0,1) = clov[18]+i*clov[19] */
	/* c(2,1, 0,1) = clov[22]+i*clov[23] */
	/* c(2,1, 1,1) = clov[32]+i*clov[33] */
	jk = 14 + (4*j) + (10*k);
	clov[(72*i)+(jk   )] = -QLA_real(QLA_elem_M(A[0][i],j,k));
	clov[(72*i)+(jk+1 )] = -QLA_imag(QLA_elem_M(A[0][i],j,k));
	clov[(72*i)+(jk+36)] = -QLA_real(QLA_elem_M(A[1][i],j,k));
	clov[(72*i)+(jk+37)] = -QLA_imag(QLA_elem_M(A[1][i],j,k));
      }
    }
  }

  QDP_reset_M(a[0]);
  QDP_reset_M(a[1]);

  QDP_destroy_M(a[0]);
  QDP_destroy_M(a[1]);

  /* Creating B[0] and B[1] : 3x3 matrices for X and Y                     */
  /* ------------------ NOTATION ------------------------------------------*/
  /* For the 3x3 part of X   : B[0] above                                  */
  /* c(1,1, 0,0) = clov[10]+i*clov[11]; c(1,1, 1,0) = clov[24]+i*clov[25]; */
  /* c(2,0, 0,0) = clov[12]+i*clov[13]; c(2,0, 1,0) = clov[26]+i*clov[27]; */
  /* c(2,1, 0,0) = clov[14]+i*clov[15]; c(2,1, 1,0) = clov[28]+i*clov[29]; */
  /* c(1,1, 0,1) = clov[18]+i*clov[19]; */
  /* c(2,0, 0,1) = clov[20]+i*clov[21]; */
  /* c(2,1, 0,1) = clov[22]+i*clov[23]; */
  /* For the 3x3 part of Y   : B[1] above                                  */
  /* c(1,3, 0,2) = clov[46]+i*clov[47]; c(1,3, 1,2) = clov[60]+i*clov[61]; */
  /* c(2,2, 0,2) = clov[48]+i*clov[49]; c(2,2, 1,2) = clov[62]+i*clov[63]; */
  /* c(2,3, 0,2) = clov[50]+i*clov[51]; c(2,3, 1,2) = clov[64]+i*clov[65]; */
  /* c(1,3, 0,3) = clov[54]+i*clov[55]; */
  /* c(2,2, 0,3) = clov[56]+i*clov[57]; */
  /* c(2,3, 0,3) = clov[58]+i*clov[59]; */

  b[0] = QDP_create_M();
  b[1] = QDP_create_M();

  f_mu_nu(xm, cs, link, 1, 2);               /* xm   = cs*F_{12}         */  
  f_mu_nu(f_mn, ct, link, 0, 3);             /* f_mn = ct*F_{03}         */

  QDP_M_eq_M_plus_M(ym, xm, f_mn, QDP_all);  /* ym = cs*F_{12}+ct*F_{03}    */
  QDP_M_meq_M(xm, f_mn, QDP_all);            /* xm = cs*F_{12}-ct*F_{03}    */

  QDP_M_eq_i_M(b[0], xm, QDP_all);           /* b[0] = i(cs*F_{12}-ct*F_{03}) */
  QDP_M_eq_i_M(b[1], ym, QDP_all);           /* b[1] = i(cs*F_{12}+ct*F_{03}) */

  f_mu_nu(f_mn, cs, link, 0, 2);             /* f_mn = cs*F_{02}           */
  QDP_M_meq_M(b[0], f_mn, QDP_all);          /* b[0] = b[0]-cs*F_{02}      */
  QDP_M_meq_M(b[1], f_mn, QDP_all);          /* b[1] = b[1]-cs*F_{02}      */

  f_mu_nu(f_mn, ct, link, 1, 3);             /* f_mn = ct*F_{13}           */
  QDP_M_meq_M(b[0], f_mn, QDP_all);          /* b[0] = b[0]-ct*F_{13}      */
  QDP_M_peq_M(b[1], f_mn, QDP_all);          /* b[1] = b[1]+ct*F_{13}      */

  /*  b[0] = i(cs*F_{12}-ct*F_{03})-(cs*F_{02}+ct*F_{13})    */
  /*  b[1] = i(cs*F_{12}+ct*F_{03})-(cs*F_{02}-ct*F_{13})    */

  QLA_ColorMatrix *B[2];
  B[0] = QDP_expose_M(b[0]);
  B[1] = QDP_expose_M(b[1]);
  for(int i=0; i<QDP_sites_on_node; i++) {
#if 1
    for(int j=0; j<3; j++) {
      for(int k=0; k<3; k++) {
	if(j<k) {
	  /* c(0,1, 1,0)* = clov[16]+i*clov[17] */
	  /* c(0,1, 2,0)* = clov[20]+i*clov[21] */
	  /* c(1,1, 2,0)* = clov[30]+i*clov[31] */
	  int jk = 12 + 10*j + 4*k;
	  clov[72*i+(jk   )] =  QLA_real(QLA_elem_M(B[0][i],j,k));
	  clov[72*i+(jk+1 )] = -QLA_imag(QLA_elem_M(B[0][i],j,k));
	  clov[72*i+(jk+36)] =  QLA_real(QLA_elem_M(B[1][i],j,k));
	  clov[72*i+(jk+37)] = -QLA_imag(QLA_elem_M(B[1][i],j,k));
	} else {
	  /* c(0,1, 0,0) = clov[ 6]+i*clov[ 7] */
	  /* c(1,1, 0,0) = clov[10]+i*clov[11] */
	  /* c(2,1, 0,0) = clov[14]+i*clov[15] */
	  /* c(1,1, 1,0) = clov[24]+i*clov[25] */
	  /* c(2,1, 1,0) = clov[28]+i*clov[29] */
	  /* c(2,1, 2,0) = clov[34]+i*clov[35] */
	  int jk = 6 + 4*j + k*(18-4*k);
	  clov[72*i+(jk   )] = QLA_real(QLA_elem_M(B[0][i],j,k));
	  clov[72*i+(jk+1 )] = QLA_imag(QLA_elem_M(B[0][i],j,k));
	  clov[72*i+(jk+36)] = QLA_real(QLA_elem_M(B[1][i],j,k));
	  clov[72*i+(jk+37)] = QLA_imag(QLA_elem_M(B[1][i],j,k));
	}
      }
    }
#else
    for(int tr=0; tr<2; tr++) {
      clov[(72*i)+6+(36*tr)]  = +QLA_real(QLA_elem_M(B[tr][i],0,0));
      clov[(72*i)+7+(36*tr)]  = +QLA_imag(QLA_elem_M(B[tr][i],0,0));
      clov[(72*i)+10+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],1,0));
      clov[(72*i)+11+(36*tr)] = +QLA_imag(QLA_elem_M(B[tr][i],1,0));
      clov[(72*i)+14+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],2,0));
      clov[(72*i)+15+(36*tr)] = +QLA_imag(QLA_elem_M(B[tr][i],2,0));
      clov[(72*i)+16+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],0,1));
      clov[(72*i)+17+(36*tr)] = -QLA_imag(QLA_elem_M(B[tr][i],0,1));
      clov[(72*i)+24+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],1,1));
      clov[(72*i)+25+(36*tr)] = +QLA_imag(QLA_elem_M(B[tr][i],1,1));
      clov[(72*i)+28+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],2,1));
      clov[(72*i)+29+(36*tr)] = +QLA_imag(QLA_elem_M(B[tr][i],2,1));
      clov[(72*i)+20+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],0,2));
      clov[(72*i)+21+(36*tr)] = -QLA_imag(QLA_elem_M(B[tr][i],0,2));
      clov[(72*i)+30+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],1,2));
      clov[(72*i)+31+(36*tr)] = -QLA_imag(QLA_elem_M(B[tr][i],1,2));
      clov[(72*i)+34+(36*tr)] = +QLA_real(QLA_elem_M(B[tr][i],2,2));
      clov[(72*i)+35+(36*tr)] = +QLA_imag(QLA_elem_M(B[tr][i],2,2));
    }
#endif
  }
  QDP_destroy_M(b[0]);
  QDP_destroy_M(b[1]);

  QDP_destroy_M(f_mn);
  QDP_destroy_M(xm);
  QDP_destroy_M(ym);

  /* 2*(3x3)=2*(18 real) = 36 real numbers are assiged to array clov  */
  /* ---------------------------------------------------------------- */
}

/* ---------------------------------------------------------------------- */
/* -------- End of the calculation of Clover Coefficients  -------------  */
/* ---------------------------------------------------------------------- */



// ************************************************************************ 

static void
clov_unpack(QLA_Complex u[6][6], REAL *p)
{
  int i, j, k;
  for(i=0; i<6; i++) {
    QLA_c_eq_r(u[i][i], p[i]);
  }
  k = 6;
  for(i=0; i<6; i++) {
    for(j=i+1; j<6; j++) {
      QLA_c_eq_ca(u[i][j], *((QLA_Complex *)(p+k)));
      QLA_c_eq_c(u[j][i], *((QLA_Complex *)(p+k)));
      k += 2;
    }
  }
}

static void
clov_pack(REAL *p, QLA_Complex u[6][6])
{
  int i, j, k;
  for(i=0; i<6; i++) {
    p[i] = QLA_real(u[i][i]);
  }
  k = 6;
  for(i=0; i<6; i++) {
    for(j=i+1; j<6; j++) {
      QLA_Complex z;
      QLA_c_eq_ca(z, u[i][j]);
      QLA_c_peq_c(z, u[j][i]);
      QLA_c_eq_r_times_c(*((QLA_Complex *)(p+k)), 0.5, z);
      k += 2;
    }
  }
}

static void
clov_invert(QLA_Complex ci[6][6], QLA_Complex c[6][6])
{
  int i, j, k;
  for(i=0; i<6; i++) {
    for(j=0; j<6; j++) {
      QLA_c_eq_r(ci[i][j], 0);
    }
    QLA_c_eq_r(ci[i][i], 1);
  }
  for(k=0; k<6; k++) {
    QLA_Complex s;
    //QLA_c_eq_r_div_c(s, 1, c[k][k]);
    {
      QLA_Real r;
      r = 1/QLA_norm2_c(c[k][k]);
      QLA_c_eq_r_times_c(s, r, c[k][k]);
      QLA_c_eq_ca(s, s);
    }
    for(j=0; j<6; j++) {
      QLA_Complex t;
      QLA_c_eq_c_times_c(t, s, c[k][j]);
      QLA_c_eq_c(c[k][j], t);
    }
    for(j=0; j<6; j++) {
      QLA_Complex t;
      QLA_c_eq_c_times_c(t, s, ci[k][j]);
      QLA_c_eq_c(ci[k][j], t);
    }
    for(i=0; i<6; i++) {
      if(i==k) continue;
      QLA_c_eq_c(s, c[i][k]);
      for(j=0; j<6; j++) {
	QLA_c_meq_c_times_c(c[i][j], s, c[k][j]);
      }
      for(j=0; j<6; j++) {
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
      printf("%f %f ", QLA_real(m[i*nc+j]), QLA_imag(m[i*nc+j]));
    }
    printf("\n");
  }
}
*/

static void
get_clovinv(QOP_FermionLinksWilson *flw, REAL kappa)
{
  QLA_Real m4 = 0.5/kappa;
  int i, j;
  if(flw->clovinv==NULL)
    QOP_malloc(flw->clovinv, REAL, QDP_sites_on_node*CLOV_REALS);
  for(i=0; i<2*QDP_sites_on_node; i++) {
    QLA_Complex cu[6][6], ciu[6][6];
    clov_unpack(cu, flw->clov+(CLOV_REALS/2)*i);
    for(j=0; j<6; j++) QLA_c_peq_r(cu[j][j], m4);
    //if(i==0) printcm(flw->clov, 3, 6);
    //if(i==0) printcm(cu, 6, 6);
    clov_invert(ciu, cu);
    clov_pack(flw->clovinv+(CLOV_REALS/2)*i, ciu);
    //if(i==0) printcm(cu, 6, 6);
    //if(i==0) printcm(ciu, 6, 6);
    //if(i==0) printcm(flw->clovinv, 3, 6);
  }
  flw->clovinvkappa = kappa;
}

QOP_FermionLinksWilson *
QOP_wilson_create_L_from_raw(REAL *links[], REAL *clov, QOP_evenodd_t evenodd)
{
  QOP_FermionLinksWilson *flw;
  QOP_GaugeField *gf;

  WILSON_INVERT_BEGIN;

  gf = QOP_create_G_from_raw(links, evenodd);
  flw = QOP_wilson_convert_L_from_qdp(gf->links, NULL);

  if(clov!=NULL) {
    QOP_malloc(flw->clov, REAL, QDP_sites_on_node*CLOV_REALS);
    memcpy(flw->clov, clov, QDP_sites_on_node*CLOV_SIZE);
  } else {
    flw->clov = NULL;
  }
  flw->clovinv = NULL;
  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qdpclov = NULL;
  flw->qopgf = gf;

  WILSON_INVERT_END;
  return flw;
}
/* --------------------------------------------------- */
/* This part is added by Bugra ----------------------- */

static QOP_FermionLinksWilson *
wilson_initialize_gauge_L()
{
  QOP_FermionLinksWilson *flw;

  WILSON_INVERT_BEGIN;

  QOP_malloc(flw          ,QOPPC(FermionLinksWilson),1);
  QOP_malloc(flw->links   ,QDPPC(ColorMatrix) *     ,4);
  QOP_malloc(flw->bcklinks,QDPPC(ColorMatrix) *     ,4);
  QOP_malloc(flw->dbllinks,QDPPC(ColorMatrix) *,     8);

  flw->dblstored = 0;
  flw->clov      = NULL;
  flw->clovinv   = NULL;
  flw->rawlinks  = NULL;
  flw->qopgf     = NULL;
  flw->qdpclov   = NULL;
  flw->eigcg.u   = NULL;

  WILSON_INVERT_END;

  return flw;
}


/* ----------------This part is added by Bugra ------------- */
/* --------------------------------------------------------- */
QOP_FermionLinksWilson *
QOP_wilson_create_L_from_G(QOP_info_t *info, 
			   QOP_wilson_coeffs_t *coeffs,
			   QOP_GaugeField *gauge)
{ 
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix        *newlinks[4];
  int                    i;

  WILSON_INVERT_BEGIN;

  /* initialize FermionLinksWilson and allocate memory---- */
  flw = wilson_initialize_gauge_L();

  /* First create QDP Color Matrices */
  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], gauge->links[i], QDP_all);
  }

  /* get the clover coefficients and put them in flw->clow */
  /* Usage : get_clov(QLA_Real *clov, QDP_ColorMatrix *link[], QLA_Real csw) */
  if(coeffs->clov_s != 0 || coeffs->clov_t != 0) {
    int nreals = QDP_sites_on_node*CLOV_REALS;
    QOP_malloc(flw->clov, REAL, nreals);
    get_clov(flw->clov, newlinks, 0.5*coeffs->clov_s, 0.5*coeffs->clov_t);
  }

  /* Check the anisotropy -------------------------------  */
  if(coeffs->aniso != 0.) {
    for(i=0; i<3; i++) {
      QLA_Real f = coeffs->aniso;
      QDP_M_eq_r_times_M(newlinks[i], &f, newlinks[i], QDP_all);
    }
  }

  /* Scale the links ------------------------------------- */
  for(i=0; i<4; i++) {
    QLA_Real f    = -0.5;
    QDP_M_eq_r_times_M(newlinks[i], &f, newlinks[i], QDP_all);
  }

  /* newlinks go to flw->links --------------------------- */
  for(i=0; i<4; i++) {
    flw->links[i] = newlinks[i];
  }

  WILSON_INVERT_END;
  return flw;
}

/* --------------------------------------------------------- */
/* --------------------------------------------------------- */

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
  int i;

  WILSON_INVERT_BEGIN;

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
  if(flw->clov) free(flw->clov);
  if(flw->clovinv) free(flw->clovinv);
  free(flw->bcklinks);
  free(flw->dbllinks);
  if(flw->eigcg.u) {
    for(i=0; i<flw->eigcg.numax; i++) {
      QDP_destroy_D(flw->eigcg.u[i]);
    }
    free(flw->eigcg.u);
    free(flw->eigcg.l);
  }
  free(flw);
  WILSON_INVERT_END;
}

QOP_FermionLinksWilson *
QOP_wilson_convert_L_from_raw(REAL *links[], REAL *clov,
			      QOP_evenodd_t evenodd)
{
  WILSON_INVERT_BEGIN;
  QOP_error("QOP_wilson_convert_L_from_raw unimplemented");
  WILSON_INVERT_END;
  return NULL;
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
  QOP_FermionLinksWilson *flw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  WILSON_INVERT_BEGIN;

  for(i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  flw = QOP_wilson_convert_L_from_qdp(newlinks, clov);
  flw->qdpclov = NULL;

  WILSON_INVERT_END;
  return flw;
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
  QOP_FermionLinksWilson *flw;
  int i;

  WILSON_INVERT_BEGIN;

  QOP_malloc(flw, QOPPC(FermionLinksWilson), 1);
  QOP_malloc(flw->links, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->bcklinks, QDPPC(ColorMatrix) *, 4);
  QOP_malloc(flw->dbllinks, QDPPC(ColorMatrix) *, 8);
  if(clov!=NULL) {
    int size = QDP_sites_on_node*CLOV_REALS;
    QOP_malloc(flw->clov, REAL, size);
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
	      flw->clov[k++] = -QLA_imag(z2); // - since we now store lower tri
	    }
	  }
	}
      }
      QDP_reset_P(clov);
    }
  } else {
    flw->clov = NULL;
  }
  flw->clovinv = NULL;

  flw->dblstored = 0;
  for(i=0; i<4; i++) {
    flw->links[i] = links[i];
  }
  // scale links
  for(i=0; i<4; i++) {
    QLA_Real f = -0.5;
    QDP_M_eq_r_times_M(flw->links[i], &f, flw->links[i], QDP_all);
  }

  //check_setup(flw);

  flw->rawlinks = NULL;
  flw->rawclov = NULL;
  flw->qopgf = NULL;
  flw->qdpclov = clov;
  flw->eigcg.u = NULL;

  WILSON_INVERT_END;
  return flw;
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
      clov(flw, kappa, out, in, QDP_all, 1);
      //printf("testd23\n");
    } else if(eo_out==QOP_EVEN) {
      clov(flw, kappa, out, in, QDP_even, 0);
    } else {
      clov(flw, kappa, out, in, QDP_odd, 0);
    }
  } else {
    if(eo_out==QOP_EVEN || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_ODD) {
	wilson_hop(flw, out, in, sign, QOP_EVEN);
      } else if(eo_in==QOP_EVEN) {
	clov(flw, kappa, out, in, QDP_even, 0);
      } else {
	wilson_hop(flw, out, in, sign, QOP_EVEN);
	clov(flw, kappa, out, in, QDP_even, 1);
      }
    }
    if(eo_out==QOP_ODD || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_EVEN) {
	wilson_hop(flw, out, in, sign, QOP_ODD);
      } else if(eo_in==QOP_ODD) {
	clov(flw, kappa, out, in, QDP_odd, 0);
      } else {
	wilson_hop(flw, out, in, sign, QOP_ODD);
	clov(flw, kappa, out, in, QDP_odd, 1);
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
  if(flw->clov==NULL) {
    QLA_Real f = 2*kappa;
    QDP_D_eq_r_times_D(out, &f, in, qdpsub(eo));
  } else {
    if( flw->clovinv==NULL || flw->clovinvkappa!=kappa ) {
      get_clovinv(flw, kappa);
    }
    QDP_D_eq_zero(out, qdpsub(eo));
    apply_clov(flw->clovinv, 0, out, in, qdpsub(eo));
  }
}

#define cmplx(x) (*((QLA_Complex *)(&(x))))

static void
apply_clov_qla(REAL *clov, QLA_Real m4, QLA_DiracFermion *restrict clov_out,
	       QLA_DiracFermion *restrict clov_in, QDP_Subset subset)
{
  //QLA_DiracFermion *clov_out, *clov_in;
  //clov_out = QDP_expose_D(out);
  //clov_in = QDP_expose_D(in);
  //{
    int x, start, end;
    if(subset==QDP_odd) start = QDP_subset_len(QDP_even);
    else start = 0;
    end = start + QDP_subset_len(subset);
    for(x=start; x<end; x++) {
      int b;
      for(b=0; b<2; b++) {
	int xb;
	xb = 36*(2*x+b);  // chiral block offset (in REALs)

#define clov_diag(i) clov[xb+i]
#define clov_offd(i) cmplx(clov[xb+6+2*i])
#define src(i) QLA_elem_D(clov_in[x],i/2,2*b+i%2)
#define dest(i) QLA_elem_D(clov_out[x],i/2,2*b+i%2)
	//flops = 6*44 = 264; bytes = 4*(36+12+12) = 240
#if 0
	{
	  QLA_Complex z0, z1, z2, z3, z4, z5;

	  QLA_c_eq_r_times_c(z0, m4+clov_diag(0), src(0));
	  QLA_c_eq_c_times_c(z1, clov_offd(0), src(0));
	  QLA_c_eq_c_times_c(z2, clov_offd(1), src(0));
	  QLA_c_eq_c_times_c(z3, clov_offd(2), src(0));
	  QLA_c_eq_c_times_c(z4, clov_offd(3), src(0));
	  QLA_c_eq_c_times_c(z5, clov_offd(4), src(0));

	  QLA_c_peq_ca_times_c(z0, clov_offd(0), src(1));
	  QLA_c_peq_r_times_c(z1, m4+clov_diag(1), src(1));
	  QLA_c_peq_c_times_c(z2, clov_offd(5), src(1));
	  QLA_c_peq_c_times_c(z3, clov_offd(6), src(1));
	  QLA_c_peq_c_times_c(z4, clov_offd(7), src(1));
	  QLA_c_peq_c_times_c(z5, clov_offd(8), src(1));

	  QLA_c_peq_ca_times_c(z0, clov_offd(1), src(2));
 	  QLA_c_peq_ca_times_c(z1, clov_offd(5), src(2));
	  QLA_c_peq_r_times_c(z2, m4+clov_diag(2), src(2));
	  QLA_c_peq_c_times_c(z3, clov_offd(9), src(2));
	  QLA_c_peq_c_times_c(z4, clov_offd(10), src(2));
	  QLA_c_peq_c_times_c(z5, clov_offd(11), src(2));

	  QLA_c_peq_ca_times_c(z0, clov_offd(2), src(3));
	  QLA_c_peq_ca_times_c(z1, clov_offd(6), src(3));
	  QLA_c_peq_ca_times_c(z2, clov_offd(9), src(3));
	  QLA_c_peq_r_times_c(z3, m4+clov_diag(3), src(3));
	  QLA_c_peq_c_times_c(z4, clov_offd(12), src(3));
	  QLA_c_peq_c_times_c(z5, clov_offd(13), src(3));

	  QLA_c_peq_ca_times_c(z0, clov_offd(3), src(4));
	  QLA_c_peq_ca_times_c(z1, clov_offd(7), src(4));
	  QLA_c_peq_ca_times_c(z2, clov_offd(10), src(4));
	  QLA_c_peq_ca_times_c(z3, clov_offd(12), src(4));
	  QLA_c_peq_r_times_c(z4, m4+clov_diag(4), src(4));
	  QLA_c_peq_c_times_c(z5, clov_offd(14), src(4));

	  QLA_c_peq_ca_times_c(z0, clov_offd(4), src(5));
	  QLA_c_peq_ca_times_c(z1, clov_offd(8), src(5));
	  QLA_c_peq_ca_times_c(z2, clov_offd(11), src(5));
	  QLA_c_peq_ca_times_c(z3, clov_offd(13), src(5));
	  QLA_c_peq_ca_times_c(z4, clov_offd(14), src(5));
	  QLA_c_peq_r_times_c(z5, m4+clov_diag(5), src(5));

	  QLA_c_peq_c(dest(0), z0);
	  QLA_c_peq_c(dest(1), z1);
	  QLA_c_peq_c(dest(2), z2);
	  QLA_c_peq_c(dest(3), z3);
	  QLA_c_peq_c(dest(4), z4);
	  QLA_c_peq_c(dest(5), z5);
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
#endif
      }
    }
    //}
    //QDP_reset_D(in);
    //QDP_reset_D(out);
}

static void
apply_clov(REAL *clov, QLA_Real m4, QDP_DiracFermion *out,
	   QDP_DiracFermion *in, QDP_Subset subset)
{
  QLA_DiracFermion *clov_out, *clov_in;
  clov_out = QDP_expose_D(out);
  clov_in = QDP_expose_D(in);
  apply_clov_qla(clov, m4, clov_out, clov_in, subset);
  QDP_reset_D(in);
  QDP_reset_D(out);
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
    for(mu=0; mu<4; mu+=QOP_wilson_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, flw->links+mu, vsrc+mu,
				  dir+mu, msgn+mu, othersubset,
				  QOP_wilson_nsvec);
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, QDP_neighbor+mu,
		   bck+mu, subset, QOP_wilson_nsvec);
    }
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
