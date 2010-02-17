#include <string.h>
#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define LU

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_dw_initQ;
extern int QOP_dw_style;
extern int QOP_dw_nsvec;
extern int QOP_dw_nvec;

static int old_style=-1;

static void
(*wilson_dslash)(QOP_FermionLinksDW *fldw,
                 QDP_DiracFermion *dest, QDP_DiracFermion *src,
                 int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash0(QOP_FermionLinksDW *fldw,
               QDP_DiracFermion *dest, QDP_DiracFermion *src,
               int sign, QDP_Subset subset, int ntmp);

static void
wilson_dslash1(QOP_FermionLinksDW *fldw,
               QDP_DiracFermion *dest, QDP_DiracFermion *src,
               int sign, QDP_Subset subset, int ntmp);

// -----------------------------------------------------------------
// Temporary variable management

// How many temporaries are we going to use?
#define NTMPSUB 2
#define NTMP (3*NTMPSUB)
#define NHTMP 16
#define NDTMP 12
static int dslash_initQ = 0; // Is dslash set up?
// The temporaries are born here
static QDP_HalfFermion *htemp[NTMP][NHTMP];
static QDP_DiracFermion *dtemp[NTMP][NDTMP];
static QDP_DiracFermion *tin[NTMP];
static QDP_HalfFermion **hf0=NULL, **hf1=NULL, *tth=NULL;
static QDP_DiracFermion *tt1=NULL, *tt2=NULL, **ttv=NULL;
// Some kind of automatic indexing; do we need this?
#define tmpnum(eo,n) ((eo)+3*((n)-1))
#define tmpsub(eo,n) tin[tmpnum(eo,n)]

#define check_setup(fldw,ls) \
{ \
  if( (!dslash_initQ) ) { \
    reset_temps(ls); \
  } \
  if( fldw->dblstored != dblstore_style(QOP_dw_style) ) { \
    double_store(fldw); \
  } \
}

static void
free_temps(int ls)
{
  if (dslash_initQ) {
    int i, j;

    QDP_destroy_H(tth);
    QDP_destroy_D(tt1);
    QDP_destroy_D(tt2);
    
    for (int s=0; s<ls; s++) {
      QDP_destroy_H(hf0[s]);
      QDP_destroy_H(hf1[s]);
      QDP_destroy_D(ttv[s]);
    }
    free(hf0);
    free(hf1);
    free(ttv);
    
    for (i=0; i<NTMP; i++)
      QDP_destroy_D(tin[i]);

    if (shiftd_style(old_style)) {
      for (i=0; i<NTMP; i++)
        for (j=0; j<NDTMP; j++)
          QDP_destroy_D(dtemp[i][j]);
    } else {
      for (i=0; i<NTMP; i++)
        for (j=0; j<NHTMP; j++)
          QDP_destroy_H(htemp[i][j]);
    }
  }
  dslash_initQ = 0;
}

static void
reset_temps(int ls)
{
  int i, j;

  free_temps(ls);
  
  if (dblstore_style(QOP_dw_style)) {
    wilson_dslash = wilson_dslash1;
  } else {
    wilson_dslash = wilson_dslash0;
  }
  
  tth = QDP_create_H();
  tt1 = QDP_create_D();
  tt2 = QDP_create_D();
  
  hf0 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  hf1 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  ttv = (QDP_DiracFermion**) malloc(ls*sizeof(QDP_DiracFermion*));
  for (int s=0; s<ls; s++) {
    hf0[s] = QDP_create_H();
    hf1[s] = QDP_create_H();
    ttv[s] = QDP_create_D();
  }
  
  for (i=0; i<NTMP; i++)
    tin[i] = QDP_create_D();

  if (shiftd_style(QOP_dw_style)) {
    for (i=0; i<NTMP; i++)
      for (j=0; j<NDTMP; j++)
        dtemp[i][j] = QDP_create_D();
  } else {
    for (i=0; i<NTMP; i++)
      for (j=0; j<NHTMP; j++)
        htemp[i][j] = QDP_create_H();
  }
  dslash_initQ = 1;
  old_style = QOP_dw_style;
}

static void
double_store( QOP_FermionLinksDW *fldw )
{
  int i;

  if (fldw->dblstored) {
    for (i=0; i<4; i++)
      QDP_destroy_M(fldw->bcklinks[i]);
    fldw->dblstored = 0;
  }

  if (dblstore_style(QOP_dw_style)) {
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<4; i++) {
      fldw->bcklinks[i] = QDP_create_M();
      fldw->dbllinks[2*i] = fldw->links[i];
      fldw->dbllinks[2*i+1] = fldw->bcklinks[i];
      QDP_M_eq_sM(m, fldw->links[i], QDP_neighbor[i],
                  QDP_backward, QDP_all);
      QDP_M_eq_Ma(fldw->bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    fldw->dblstored = dblstore_style(QOP_dw_style);
  }
}

// -----------------------------------------------------------------
// Gauge field management

QOP_FermionLinksDW *
QOP_dw_create_L_from_raw( REAL *links[], QOP_evenodd_t evenodd )
{
  QOP_FermionLinksDW *fldw;
  QOP_GaugeField *gf;

  DW_INVERT_BEGIN;

  gf = QOP_create_G_from_raw(links, evenodd);
  fldw = QOP_dw_convert_L_from_qdp(gf->links);

  fldw->qopgf = gf;

  DW_INVERT_END;
  return fldw;
}

static QOP_FermionLinksDW *
dw_initialize_gauge_L()
{
  QOP_FermionLinksDW *fldw;

  DW_INVERT_BEGIN;

  QOP_malloc(fldw          ,QOPPC(FermionLinksDW),1);
  QOP_malloc(fldw->links   ,QDPPC(ColorMatrix) * ,4);
  QOP_malloc(fldw->bcklinks,QDPPC(ColorMatrix) * ,4);
  QOP_malloc(fldw->dbllinks,QDPPC(ColorMatrix) * ,8);

  fldw->dblstored = 0;
  fldw->qopgf     = NULL;

  DW_INVERT_END;

  return fldw;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_G( QOP_info_t *info, 
                   QOP_dw_coeffs_t *coeffs,
                    QOP_GaugeField *gauge )
{ 
  QOP_FermionLinksDW *fldw;
  QDP_ColorMatrix    *newlinks[4];
  int                i;

  DW_INVERT_BEGIN;

  fldw = dw_initialize_gauge_L();
  
  for (i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], gauge->links[i], QDP_all);
  }

  DW_INVERT_END;
  return fldw;
}

void
QOP_dw_extract_L_to_raw( REAL *links[],
           QOP_FermionLinksDW *src,
                QOP_evenodd_t evenodd )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_extract_L_to_raw unimplemented.");
  DW_INVERT_END;
}

void
QOP_dw_destroy_L( QOP_FermionLinksDW *fldw )
{
  int i;

  DW_INVERT_BEGIN;

  if (fldw->qopgf) {
    QOP_destroy_G(fldw->qopgf);
  } else {
    for (i=0; i<4; i++) QDP_destroy_M(fldw->links[i]);
  }
  free(fldw->links);
  if (fldw->dblstored) {
    for (i=0; i<4; i++) QDP_destroy_M(fldw->bcklinks[i]);
  }
  free(fldw->bcklinks);
  free(fldw->dbllinks);
  /*if (fldw->eigcg.u) {
    for (i=0; i<fldw->eigcg.numax; i++)
      QDP_destroy_D(fldw->eigcg.u[i]);
    free(fldw->eigcg.u);
    free(fldw->eigcg.l);
  }*/
  free(fldw);
  DW_INVERT_END;
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_raw( REAL *links[],
                  QOP_evenodd_t evenodd )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_convert_L_from_raw unimplemented");
  DW_INVERT_END;
  return NULL;
}

void
QOP_dw_convert_L_to_raw( REAL ***links,
           QOP_FermionLinksDW *src,
                QOP_evenodd_t evenodd )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_convert_L_to_raw unimplemented");
  DW_INVERT_END;
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_G( QOP_info_t *info,
                    QOP_dw_coeffs_t *coeffs,
                     QOP_GaugeField *gauge )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_convert_L_from_G unimplemented");
  DW_INVERT_END;
  return NULL;
}

QOP_GaugeField *
QOP_dw_convert_L_to_G( QOP_FermionLinksDW *links )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_convert_L_to_G unimplemented");
  DW_INVERT_END;
  return NULL;
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_qdp( QDP_ColorMatrix *links[] )
{
  QOP_FermionLinksDW *fldw;
  QDP_ColorMatrix *newlinks[4];
  int i;

  DW_INVERT_BEGIN;

  for (i=0; i<4; i++) {
    newlinks[i] = QDP_create_M();
    QDP_M_eq_M(newlinks[i], links[i], QDP_all);
  }

  fldw = QOP_dw_convert_L_from_qdp(newlinks);

  DW_INVERT_END;
  return fldw;
}

void
QOP_dw_extract_L_to_qdp( QDP_ColorMatrix *links[],
                      QOP_FermionLinksDW *src )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_extract_L_to_qdp unimplemented");
  DW_INVERT_END;
}

QOP_FermionLinksDW *
QOP_dw_convert_L_from_qdp( QDP_ColorMatrix *links[] )
{
  QOP_FermionLinksDW *fldw;
  int i;

  DW_INVERT_BEGIN;

  QOP_malloc(fldw,           QOPPC(FermionLinksDW), 1);
  QOP_malloc(fldw->links,    QDPPC(ColorMatrix)*,   4);
  QOP_malloc(fldw->bcklinks, QDPPC(ColorMatrix)*,   4);
  QOP_malloc(fldw->dbllinks, QDPPC(ColorMatrix)*,   8);

  fldw->dblstored = 0;
  for (i=0; i<4; i++) {
    fldw->links[i] = links[i];
  }

  fldw->qopgf = NULL;

  DW_INVERT_END;
  return fldw;
}

void
QOP_dw_convert_L_to_qdp( QDP_ColorMatrix ***links,
                      QOP_FermionLinksDW *src )
{
  DW_INVERT_BEGIN;
  QOP_error("QOP_dw_convert_L_to_qdp unimplemented");
  DW_INVERT_END;
}


// -----------------------------------------------------------------
// Wilson operators

static void set_out(void);

static void
wilson_dslash0( QOP_FermionLinksDW *fldw,
                  QDP_DiracFermion *dest,
                  QDP_DiracFermion *src,
                               int sign,
                        QDP_Subset subset,
                               int ntmp )
{
  //printf("In dslash0\n");
  int mu;
  QDP_DiracFermion *vsrc[4];
  QDP_DiracFermion *vdest[4];
  QDP_ShiftDir fwd[4], bck[4];
  int dir[4], sgn[4], msgn[4];
  QDP_Subset othersubset;

  sign = -sign;

  for (mu=0; mu<4; mu++) {
    vsrc[mu] = src;
    vdest[mu] = dest;
    fwd[mu] = QDP_forward;
    bck[mu] = QDP_backward;
    dir[mu] = mu;
    sgn[mu] = sign;
    msgn[mu] = -sign;
  }

  if (subset==QDP_even) othersubset = QDP_odd;
  else if (subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  /* Take DW projection for src displaced in up direction, gather
     it to "our site" */

  printf0("dslash0\n");
  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_D_veq_sD\n");
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, QDP_neighbor+mu,
                   fwd+mu, subset, QOP_dw_nsvec);
      printf0("end QDP_D_veq_sD\n");
    }
  } else {
    for (mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_H_veq_spproj_D\n");
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
                         othersubset, QOP_dw_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, QDP_neighbor+mu,
                   fwd+mu, subset, QOP_dw_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Take DW projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */

  printf0("dslash0 - back\n");
  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<4; mu+=QOP_dw_nsvec) {
      QDP_D_veq_spproj_Ma_times_D(dtemp[ntmp]+8+mu, fldw->links+mu, vsrc+mu,
                                  dir+mu, msgn+mu, othersubset, 
                                  QOP_dw_nsvec);
      QDP_D_veq_sD(dtemp[ntmp]+4+mu, dtemp[ntmp]+8+mu, QDP_neighbor+mu,
                   bck+mu, subset, QOP_dw_nsvec);
    }
  } else {
    for (mu=0; mu<4; mu+=QOP_dw_nsvec) {
      printf0("QDP_H_veq_spproj_Ma_times_D\n");
      QDP_H_veq_spproj_Ma_times_D(htemp[ntmp]+12+mu, fldw->links+mu,
                                  vsrc+mu, dir+mu, msgn+mu, othersubset,
                                  QOP_dw_nsvec);
      printf0("QDP_H_veq_sH\n");
      QDP_H_veq_sH(htemp[ntmp]+4+mu, htemp[ntmp]+12+mu, QDP_neighbor+mu,
                   bck+mu, subset, QOP_dw_nsvec);
      printf0("end QDP_H_veq_sH\n");
    }
  }

  /* Set dest to zero */
  /* Take DW projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */

  printf0("dslash0 - fwd\n");
  set_out();

  if (shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, fldw->links+mu, dtemp[ntmp]+mu,
                                  dir+mu, sgn+mu, subset, QOP_dw_nvec);
    }
  } else {
    for (mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, fldw->links+mu, htemp[ntmp]+mu,
                                   dir+mu, sgn+mu, subset, QOP_dw_nvec);
    }
  }

  /* Take DW projection for src displaced in down direction,
     expand it, and add to dest */

  printf0("dslash0 - back\n");
  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_D(vdest+mu, dtemp[ntmp]+4+mu, subset, QOP_dw_nvec);
    }
  } else {
    for (mu=0; mu<4; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_H(vdest+mu, htemp[ntmp]+4+mu, dir+mu, msgn+mu,
                           subset, QOP_dw_nvec);
    }
  }

  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for (mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
}

static void
wilson_dslash1( QOP_FermionLinksDW *fldw,
                  QDP_DiracFermion *dest,
                  QDP_DiracFermion *src,
                               int sign,
                        QDP_Subset subset,
                               int ntmp )
{
  int mu;
  QDP_DiracFermion *vsrc[8];
  QDP_DiracFermion *vdest[8];
  QDP_Shift sh[8];
  QDP_ShiftDir sd[8];
  int dir[8], sgn[8];
  QDP_Subset othersubset=QDP_all;

  if (!shiftd_style(QOP_dw_style)) {
    if (subset==QDP_even) othersubset = QDP_odd;
    else if (subset==QDP_odd) othersubset = QDP_even;
    else othersubset = QDP_all;
  }

  sign = -sign;

  for (mu=0; mu<4; mu++) {
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

  /* Take DW projection for src displaced in up direction, gather
     it to "our site" */

  if (shiftd_style(QOP_dw_style)) {
    for(mu=0; mu<8; mu+=QOP_dw_nsvec) {
      QDP_D_veq_sD(dtemp[ntmp]+mu, vsrc+mu, sh+mu, sd+mu, subset,
                   QOP_dw_nsvec);
    }
  } else {
    for (mu=0; mu<8; mu+=QOP_dw_nsvec) {
      QDP_H_veq_spproj_D(htemp[ntmp]+8+mu, vsrc+mu, dir+mu, sgn+mu,
                         othersubset, QOP_dw_nsvec);
      QDP_H_veq_sH(htemp[ntmp]+mu, htemp[ntmp]+8+mu, sh+mu, sd+mu, subset,
                   QOP_dw_nsvec);
    }
  }

  /* Set dest to zero */
  /* Take DW projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add to dest */

  set_out();

  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<8; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_spproj_M_times_D(vdest+mu, fldw->dbllinks+mu,
                                  dtemp[ntmp]+mu, dir+mu, sgn+mu,
                                  subset, QOP_dw_nvec);
    }
  } else {
    for (mu=0; mu<8; mu+=QOP_dw_nvec) {
      QDP_D_vpeq_sprecon_M_times_H(vdest+mu, fldw->dbllinks+mu,
                                   htemp[ntmp]+mu, dir+mu, sgn+mu,
                                   subset, QOP_dw_nvec);
    }
  }

  if (shiftd_style(QOP_dw_style)) {
    for (mu=0; mu<8; mu++) {
      QDP_discard_D(dtemp[ntmp][mu]);
    }
  } else {
    for (mu=0; mu<8; mu++) {
      QDP_discard_H(htemp[ntmp][mu]);
    }
  }
}

// -----------------------------------------------------------------
// Even-odd operators

// Even-odd temporaries
static QLA_Real mw;
static int dirs[2], signs[2], lupass=0;
static QDP_DiracFermion *vout[2], *vin[2], *out0, *in0;
static QDP_Subset lusubset;

static void
set_out()
{
#ifdef LU
  if (lupass==0) {
    QDP_D_eq_zero(out0, lusubset);
  } else {
    QDP_D_eq_r_times_D(out0, &mw, in0, lusubset);
    QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, lusubset, 2);
  }
#else
  QDP_D_eq_r_times_D(out0, &mw, in0, QDP_all);
  QDP_D_vmeq_spproj_D(vout, vin, dirs, signs, QDP_all, 2);
#endif
}

void
QOPPC(Qxy)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL m0, REAL M, int ls, int sign,
           QDP_Subset subset, QDP_Subset osubset)
{
  QLA_Real half = -0.5;
  int i, ntmp = (subset==QDP_odd?0:1);
  lupass = 0;
  lusubset = subset;
  for (i=0; i<ls; i++) {
    out0 = out[i];
    QDP_D_eq_r_times_D(tt1, &half, in[i], osubset);
    wilson_dslash(fldw, out[i], tt1, sign, subset, ntmp);
  }
}

void
QOPPC(Qeo)(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
           QDP_DiracFermion *in[], REAL m0, REAL M, int ls)
{
  QOPPC(Qxy)(fldw, out, in, m0, M, ls, 1, QDP_even, QDP_odd);
}

void
QOPPC(Qoe)(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
           QDP_DiracFermion *in[], REAL m0, REAL M, int ls)
{
  QOPPC(Qxy)(fldw, out, in, m0, M, ls, 1, QDP_odd, QDP_even);
}

void
QOPPC(Seo)(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
           QDP_DiracFermion *in[], REAL m0, REAL M, int ls)
{
  QOPPC(Qxy)(fldw, out, in, m0, M, ls, -1, QDP_even, QDP_odd);
}

void
QOPPC(Soe)(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
           QDP_DiracFermion *in[], REAL m0, REAL M, int ls)
{
  QOPPC(Qxy)(fldw, out, in, m0, M, ls, -1, QDP_odd, QDP_even);
}

void
QOPPC(Qxxinv)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
       REAL m0, REAL M, int ls, int sign, QDP_Subset subset)
{
  QLA_Real fac, den, mba;
  QDP_DiracFermion *tin[2];
  QDP_HalfFermion *thf[2];
  int dirs[2]={4,4}, sgns[2]={1,-1};
  int i;

  sign = -sign;
  sgns[0] *= sign;
  sgns[1] *= sign;

  mba = 1.0/m0;
  den = 1.0/(1.0 + (M/m0)*pow(mba, ls-1));
  fac = -(M/m0)*den;
  QDP_H_eq_zero(hf0[ls-1], subset);
  for (i=0; i<ls; i++) {
    thf[0] = hf0[i];
    if (i==ls-1) thf[0] = tth;
    thf[1] = hf1[i];
    tin[0] = in[i];
    tin[1] = in[i];
    QDP_H_veq_spproj_D(thf, tin, dirs, sgns, subset, 2);

    // L_A^-1
    if (i<ls-1) {
      QDP_H_peq_r_times_H(hf0[ls-1], &fac, hf0[i], subset);
      fac *= mba;
    } else {
      QDP_H_peq_r_times_H(hf0[ls-1], &den, tth, subset);
    }

    // L_B^-1
    if(i==0) {
      QDP_H_eq_r_times_H(hf1[i], &mba, hf1[i], subset);
    } else {
      QDP_H_peq_H(hf1[i], hf1[i-1], subset);
      QDP_H_eq_r_times_H(hf1[i], &mba, hf1[i], subset);
    }
  }

  fac = -(M/m0)*pow(mba, ls-2);
  for (i=ls-1; i>=0; i--) {
    // R_A^-1
    if (i==ls-1) {
      QDP_H_eq_r_times_H(hf0[i], &mba, hf0[i], subset);
    } else {
      QDP_H_peq_H(hf0[i], hf0[i+1], subset);
      QDP_H_eq_r_times_H(hf0[i], &mba, hf0[i], subset);
    }

    // R_B^-1
    if (i==ls-1) {
      QDP_H_eq_r_times_H(hf1[i], &den, hf1[i], subset);
    } else {
      QDP_H_peq_r_times_H(hf1[i], &fac, hf1[ls-1], subset);
      fac *= m0;
    }

    QDP_D_eq_sprecon_H(out[i], hf0[i], dirs[0], sgns[0], subset);
    QDP_D_peq_sprecon_H(out[i], hf1[i], dirs[1], sgns[1], subset);
  }
}

void
QOPPC(Qooinv)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL m0, REAL M, int ls)
{
  QOPPC(Qxxinv)(out, in, m0, M, ls, 1, QDP_odd);
}

void
QOPPC(Qeeinv)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL m0, REAL M, int ls)
{
  QOPPC(Qxxinv)(out, in, m0, M, ls, 1, QDP_even);
}

void
QOPPC(Sooinv)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL m0, REAL M, int ls)
{
  QOPPC(Qxxinv)(out, in, m0, M, ls, -1, QDP_odd);
}

void
QOPPC(Seeinv)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL m0, REAL M, int ls)
{
  QOPPC(Qxxinv)(out, in, m0, M, ls, -1, QDP_even);
}

void
QOPPC(QoeQeeinv)(QOP_FermionLinksDW *fldw, QDP_DiracFermion *out[],
                 QDP_DiracFermion *in[], REAL m0, REAL M, int ls)
{
  QOPPC(Qeeinv)(out, in, m0, M, ls);
  QOPPC(Qoe)(fldw, out, out, m0, M, ls);
}

void
QOPPC(Qxx)(QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           QLA_Real m0, QLA_Real M, int ls, int sign,
           QDP_Subset subset, int add)
{
  int i;
  sign = -sign;
  for (i=0; i<ls; i++) {
    if (add) {
      QDP_D_peq_r_times_D(out[i], &m0, in[i], subset);
    } else {
      QDP_D_eq_r_times_D(out[i], &m0, in[i], subset);
    }

    if (i<ls-1) {
      QDP_D_meq_spproj_D(out[i], in[i+1], 4, sign, subset);
    } else {
      QDP_D_eq_r_times_D(tt1, &M, in[0], subset);
      QDP_D_peq_spproj_D(out[i], tt1, 4, sign, subset);
    }
    if (i==0) {
      QDP_D_eq_r_times_D(tt1, &M, in[ls-1], subset);
      QDP_D_peq_spproj_D(out[i], tt1, 4, -sign, subset);
    } else {
      QDP_D_meq_spproj_D(out[i], in[i-1], 4, -sign, subset);
    }
  }
}

// -----------------------------------------------------------------
// Domain-wall operators

void
QOPPC(dw_dslash1)(QOP_FermionLinksDW *fldw,
                  QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                  int sign, QDP_Subset subset, QDP_Subset othersubset,
                  QLA_Real m0, QLA_Real M, int ls)
{
#ifdef LU
  int i;

  if (sign>0) {
    QOPPC(Qeo)(fldw, out, in, m0, M, ls);
    QOPPC(Qeeinv)(out, out, m0, M, ls);
    QOPPC(Qoe)(fldw, out, out, m0, M, ls);
    QOPPC(Qooinv)(out, out, m0, M, ls);
  } else {
    QOPPC(Sooinv)(out, in, m0, M, ls);
    QOPPC(Seo)(fldw, out, out, m0, M, ls);
    QOPPC(Seeinv)(out, out, m0, M, ls);
    QOPPC(Soe)(fldw, out, out, m0, M, ls);
  }
  for (i=0; i<ls; i++) {
    QDP_D_eq_D_minus_D(out[i], in[i], out[i], subset);
  }

#else

  QLA_Real half = -0.5;
  QLA_Real mf = -M;
  int i;

  mw = m0;
  dirs[0] = dirs[1] = 4;
  signs[0] = -sign;
  signs[1] = sign;
  for (i=0; i<ls; i++) {
    vout[0] = vout[1] = out[i];
    vin[0] = in[(i+1)%ls];
    vin[1] = in[(i-1+ls)%ls];
    out0 = out[i];
    in0 = in[i];
    if (i==0) {
      QDP_D_eq_r_times_D(tt2, &mf, in[ls-1], QDP_all);
      vin[1] = tt2;
    } else if (i==ls-1) {
      QDP_D_eq_r_times_D(tt2, &mf, in[0], QDP_all);
      vin[0] = tt2;
    }
    QDP_D_eq_r_times_D(tt1, &half, in[i], QDP_all);
    wilson_dslash(fldw, out[i], tt1, sign, QDP_all, 0);
  }

#endif
}

void
QOPPC(dw_dslash2)(QOP_FermionLinksDW *fldw,
                  QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                  QDP_Subset subset, QDP_Subset othersubset,
                  QLA_Real m0, QLA_Real M, int ls)
{
  QOPPC(dw_dslash1)(fldw, ttv,  in,  1, subset, othersubset, m0, M, ls);
  QOPPC(dw_dslash1)(fldw, out, ttv, -1, subset, othersubset, m0, M, ls);
}

void
QOP_dw_dslash_qdp( QOP_info_t *info,
           QOP_FermionLinksDW *links,
                         REAL M5,
                         REAL m,
                          int sign,
             QDP_DiracFermion *out[],
             QDP_DiracFermion *in[],
                          int ls,
                QOP_evenodd_t eo_out,
                QOP_evenodd_t eo_in)
{
  check_setup(links,ls);
 
  if ( eo_out==QOP_EVEN || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_ODD ) {
      QOPPC(Qxy)(links, out, in, M5, m, ls, sign, QDP_even, QDP_odd);
    } else if ( eo_in==QOP_EVEN ) {
      QOPPC(Qxx)(out, in, M5, m, ls, sign, QDP_even, 0);
    } else {
      QOPPC(Qxy)(links, out, in, M5, m, ls, sign, QDP_even, QDP_odd);
      QOPPC(Qxx)(out, in, M5, m, ls, sign, QDP_even, 1);
    }
  }
 
  if ( eo_out==QOP_ODD || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_EVEN ) {
      QOPPC(Qxy)(links, out, in, M5, m, ls, sign, QDP_odd, QDP_even);
    } else if( eo_in==QOP_ODD ) {
      QOPPC(Qxx)(out, in, M5, m, ls, sign, QDP_odd, 0);
    } else {
      QOPPC(Qxy)(links, out, in, M5, m, ls, sign, QDP_odd, QDP_even);
      QOPPC(Qxx)(out, in, M5, m, ls, sign, QDP_odd, 1);
    }
  }
}

void
QOP_dw_dslash( QOP_info_t *info,
       QOP_FermionLinksDW *links,
                     REAL M5,
                     REAL m,
                      int sign,
         QOP_DiracFermion *out_pt[],
         QOP_DiracFermion *in_pt[],
                      int ls,
            QOP_evenodd_t eo_out,
            QOP_evenodd_t eo_in)
{
  QDP_DiracFermion *qo[ls], *qi[ls];
  for (int i=0; i<ls; i++) {
    qo[i] = out_pt[i]->df;
    qi[i] = in_pt[i]->df;
  }
  QOP_dw_dslash_qdp(info, links, M5, m, sign, qo, qi, ls, eo_out, eo_in);
}

void QOP_dw_diaginv_qdp( QOP_info_t *info,
                 QOP_FermionLinksDW *fldw,
                               REAL M5,
                               REAL m,
                   QDP_DiracFermion **out,
                   QDP_DiracFermion **in,
                                int ls,
                      QOP_evenodd_t eo)
{
  if (eo==QOP_EVEN || eo==QOP_EVENODD)
      QOPPC(Qxxinv)(out, in, M5, m, ls, 1, QDP_even);
  if (eo==QOP_ODD || eo==QOP_EVENODD)
      QOPPC(Qxxinv)(out, in, M5, m, ls, 1, QDP_odd);
}
