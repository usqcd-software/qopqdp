#include <string.h>
#include <math.h>
#include <qop_internal.h>

//#define printf0 QOP_printf0
#define printf0(...)

#define SDC_DEBUG 0

#define LU

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_dw_initQ;
extern int QOP_dw_style;
extern int QOP_dw_nsvec;
extern int QOP_dw_nvec;

static int old_style=-1;

// -----------------------------------------------------------------
// Temporary variable management

static int dslash_initQ = 0; // Is dslash set up?
// The temporaries are born here
static QDP_HalfFermion **thv0=NULL, **thv1=NULL, *tmph=NULL; // Used in Qxxinv
static QDP_DiracFermion *tmpd=NULL, // Used in Qxx
                        **tdv=NULL; // Used in dslash2

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

    QDP_destroy_H(tmph);
    QDP_destroy_D(tmpd);
    
    for (int s=0; s<ls; s++) {
      QDP_destroy_H(thv0[s]);
      QDP_destroy_H(thv1[s]);
      QDP_destroy_D(tdv[s]);
    }
    free(thv0);
    free(thv1);
    free(tdv);
  }
  dslash_initQ = 0;
}

static void
reset_temps(int ls)
{
  free_temps(ls);
  
  tmph = QDP_create_H();
  tmpd = QDP_create_D();
  
  thv0 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  thv1 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  tdv  = (QDP_DiracFermion**) malloc(ls*sizeof(QDP_DiracFermion*));
  for (int s=0; s<ls; s++) {
    thv0[s] = QDP_create_H();
    thv1[s] = QDP_create_H();
    tdv[s]  = QDP_create_D();
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
  
  fldw = QOP_dw_convert_L_from_qdp(newlinks);

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
// Even-odd operators
// Note that these operators are functions of M0 = M5-5, not M5
void
QOPPC(Qxy)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL M0, REAL mq, int ls, int sign,
           QDP_Subset subset, QDP_Subset osubset)
// Note the ordering of the subsets!
// It matches the ordering of xy, but not the standard ordering of QDP.
{
  check_setup(fldw,ls);
  int s;
  QLA_Real half = 0.5;
  QOP_FermionLinksWilson flw;
  flw.clovinvkappa = 0;//1/(2*M0+8);
  flw.dblstored = fldw->dblstored;
  flw.links = fldw->links;
  flw.bcklinks = fldw->bcklinks;
  flw.dbllinks = fldw->dbllinks;
  flw.qopgf = fldw->qopgf;
  flw.qdpclov = NULL;
  flw.clov = NULL;
  flw.clovinv = NULL;
  flw.rawlinks = fldw->raw;
  flw.rawclov = NULL;
  //flw.eigcg;
  for (s=0; s<ls; s++) {
#if SDC_DEBUG 
QLA_Real tmp; 
QDP_r_eq_norm2_D(&tmp, in[s], subset);
printf("Qxy %d: %g -> ",s,sqrt(tmp));
#endif
    QDP_D_eq_r_times_D(tmpd, &half, in[s], subset);
    QOP_wilson_dslash_qdp(NULL, &flw,
                          1.0/(2*M0+8), sign,
                		      out[s], tmpd,
                		      (osubset==QDP_even?QOP_EVEN:QOP_ODD),
                		      ( subset==QDP_even?QOP_EVEN:QOP_ODD));
#if SDC_DEBUG
QDP_r_eq_norm2_D(&tmp, out[s], osubset);
printf("%g\n",sqrt(tmp));
#endif
    fldw->dblstored = flw.dblstored;
    fldw->links = flw.links;
    fldw->bcklinks = flw.bcklinks;
    fldw->dbllinks = flw.dbllinks;
    fldw->qopgf = flw.qopgf;
    fldw->raw = flw.rawlinks;  
  }
}

void
QOPPC(Qeo)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL M0, REAL mq, int ls)
{
  QOPPC(Qxy)(fldw, out, in, M0, mq, ls, 1, QDP_even, QDP_odd);
}

void
QOPPC(Qoe)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL M0, REAL mq, int ls)
{
  QOPPC(Qxy)(fldw, out, in, M0, mq, ls, 1, QDP_odd, QDP_even);
}

void
QOPPC(Seo)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL M0, REAL mq, int ls)
{
  QOPPC(Qxy)(fldw, out, in, M0, mq, ls, -1, QDP_even, QDP_odd);
}

void
QOPPC(Soe)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           REAL M0, REAL mq, int ls)
{
  QOPPC(Qxy)(fldw, out, in, M0, mq, ls, -1, QDP_odd, QDP_even);
}

void
QOPPC(Qxxinv)(QOP_FermionLinksDW *fldw,
              QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL M0, REAL mq, int ls, int sign, QDP_Subset subset)
{
  QLA_Real fac, den, mba;
  QDP_DiracFermion *tin[2];
  QDP_HalfFermion *thf[2];
  int dirs[2]={4,4}, sgns[2]={-sign,sign};
  int s;
  
  check_setup(fldw,ls);

  mba = 1.0/M0;
  den = 1.0/(1.0 + (mq/M0)*pow(mba, ls-1));
  fac = -(mq/M0)*den;
  QDP_H_eq_zero(thv0[ls-1], subset);
  for (s=0; s<ls; s++) {
    thf[0] = (s<ls-1?thv0[s]:tmph);
    thf[1] = thv1[s];
    tin[0] = in[s];
    tin[1] = in[s];
    QDP_H_veq_spproj_D(thf, tin, dirs, sgns, subset, 2);

    // L_A^-1
    if (s<ls-1) {
      QDP_H_peq_r_times_H(thv0[ls-1], &fac, thv0[s], subset);
      fac *= mba;
    } else {
      QDP_H_peq_r_times_H(thv0[ls-1], &den, tmph, subset);
    }

    // L_B^-1
    if(s>0) {
      QDP_H_peq_H(thv1[s], thv1[s-1], subset);
      QDP_H_eq_r_times_H(thv1[s], &mba, thv1[s], subset);
    } else {
      QDP_H_eq_r_times_H(thv1[s], &mba, thv1[s], subset);
    }
  }

  fac = -(mq/M0)*pow(mba, ls-2);
  for (s=ls-1; s>=0; s--) {
    // R_A^-1
    if (s==ls-1) {
      QDP_H_eq_r_times_H(thv0[s], &mba, thv0[s], subset);
    } else {
      QDP_H_peq_H(thv0[s], thv0[s+1], subset);
      QDP_H_eq_r_times_H(thv0[s], &mba, thv0[s], subset);
    }

    // R_B^-1
    if (s==ls-1) {
      QDP_H_eq_r_times_H(thv1[s], &den, thv1[s], subset);
    } else {
      QDP_H_peq_r_times_H(thv1[s], &fac, thv1[ls-1], subset);
      fac *= M0;
    }

    QDP_D_eq_sprecon_H( out[s], thv0[s], dirs[0], sgns[0], subset);
    QDP_D_peq_sprecon_H(out[s], thv1[s], dirs[1], sgns[1], subset);
  }
}

void
QOPPC(Qooinv)(QOP_FermionLinksDW *fldw,
              QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL M0, REAL mq, int ls)
{
  QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, 1, QDP_odd);
}

void
QOPPC(Qeeinv)(QOP_FermionLinksDW *fldw,
              QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL M0, REAL mq, int ls)
{
  QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, 1, QDP_even);
}

void
QOPPC(Sooinv)(QOP_FermionLinksDW *fldw,
              QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL M0, REAL mq, int ls)
{
  QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, -1, QDP_odd);
}

void
QOPPC(Seeinv)(QOP_FermionLinksDW *fldw,
              QDP_DiracFermion *out[], QDP_DiracFermion *in[],
              REAL M0, REAL mq, int ls)
{
  QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, -1, QDP_even);
}

void
QOPPC(QoeQeeinv)(QOP_FermionLinksDW *fldw,
                 QDP_DiracFermion *out[], QDP_DiracFermion *in[],
                 REAL M0, REAL mq, int ls)
{
  QOPPC(Qeeinv)(fldw, out, in, M0, mq, ls);
  QOPPC(Qoe)(fldw, out, out, M0, mq, ls);
}

void
QOPPC(Qxx)(QOP_FermionLinksDW *fldw,
           QDP_DiracFermion *out[], QDP_DiracFermion *in[],
           QLA_Real M0, QLA_Real mq, int ls, int sign,
           QDP_Subset subset, int add)
{
  check_setup(fldw,ls);
  for (int s=0; s<ls; s++) {
    if (add) {
      QDP_D_peq_r_times_D(out[s], &M0, in[s], subset);
    } else {
      QDP_D_eq_r_times_D(out[s], &M0, in[s], subset);
    }
#if SDC_DEBUG
QLA_Real tmp;
QDP_r_eq_norm2_D(&tmp, in[s], subset);
printf("Qxx(M0) %d: %g -> ",s,sqrt(tmp));
QDP_r_eq_norm2_D(&tmp, out[s], subset);
printf("%g\n",sqrt(tmp));
#endif
    if (s<ls-1) {
      QDP_D_meq_spproj_D(out[s], in[s+1], 4, -sign, subset);
    } else {
      QDP_D_eq_r_times_D(tmpd, &mq, in[0], subset);
      QDP_D_peq_spproj_D(out[s], tmpd, 4, -sign, subset);
    }
    if (s>0) {
      QDP_D_meq_spproj_D(out[s], in[s-1], 4, sign, subset);
    } else {
      QDP_D_eq_r_times_D(tmpd, &mq, in[ls-1], subset);
      QDP_D_peq_spproj_D(out[s], tmpd, 4, sign, subset);
    }
#if SDC_DEBUG
QDP_r_eq_norm2_D(&tmp, out[s], subset);
printf("Qxx(mq) %d: %g\n",s,sqrt(tmp));
#endif
  }
}

// -----------------------------------------------------------------
// Domain-wall operators

void
QOPPC(dw_dslash_qdp)( QOP_info_t *info,
              QOP_FermionLinksDW *fldw,
                            REAL M5,
                            REAL mq,
                             int sign,
                QDP_DiracFermion *out[],
                QDP_DiracFermion *in[],
                             int ls,
                   QOP_evenodd_t eo_out,
                   QOP_evenodd_t eo_in)
{
  REAL M0 = M5-5;
  check_setup(fldw,ls);
 
  if ( eo_out==QOP_EVEN || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_ODD ) {
      QOPPC(Qxy)(fldw, out, in, M0, mq, ls, sign, QDP_odd, QDP_even); // Qoe
    } else if ( eo_in==QOP_EVEN ) {
      QOPPC(Qxx)(fldw, out, in, M0, mq, ls, sign, QDP_even, 0);       // Qee
    } else {
#if SDC_DEBUG
printf("Checkpoint OE/EE %g:%g %g %d %d\n",M5,M0, mq, ls, sign);
#endif
      QOPPC(Qxy)(fldw, out, in, M0, mq, ls, sign, QDP_odd, QDP_even); // Qoe
      QOPPC(Qxx)(fldw, out, in, M0, mq, ls, sign, QDP_even, 1);       // Qee
    }
  }

  if ( eo_out==QOP_ODD || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_EVEN ) {
      QOPPC(Qxy)(fldw, out, in, M0, mq, ls, sign, QDP_even, QDP_odd); // Qeo
    } else if( eo_in==QOP_ODD ) {
      QOPPC(Qxx)(fldw, out, in, M0, mq, ls, sign, QDP_odd, 0);        // Qoo
    } else {
#if SDC_DEBUG
printf("Checkpoint EO/OO %g:%g %g %d %d\n",M5,M0, mq, ls, sign);
#endif
      QOPPC(Qxy)(fldw, out, in, M0, mq, ls, sign, QDP_even, QDP_odd); // Qeo
      QOPPC(Qxx)(fldw, out, in, M0, mq, ls, sign, QDP_odd, 1);        // Qoo
    }
  }

}

void
QOPPC(dw_dslash)( QOP_info_t *info,
          QOP_FermionLinksDW *fldw,
                        REAL M5,
                        REAL mq,
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
  QOPPC(dw_dslash_qdp)(info, fldw, M5, mq, sign, qo, qi, ls, eo_out, eo_in);
}

void
QOPPC(dw_dslash2_qdp)( QOP_info_t *info,
               QOP_FermionLinksDW *fldw,
                             REAL M5,
                             REAL mq,
                 QDP_DiracFermion *out[],
                 QDP_DiracFermion *in[],
                              int ls,
                    QOP_evenodd_t eo_out,
                    QOP_evenodd_t eo_in) {
  check_setup(fldw,ls);
  QOPPC(dw_dslash_qdp)(info, fldw, M5, mq,  1, tdv, in, ls,
                       QOP_EVENODD, eo_in);
  QOPPC(dw_dslash_qdp)(info, fldw, M5, mq, -1, out, tdv, ls,
                       eo_out, QOP_EVENODD);
}

void
QOPPC(dw_dslash2)( QOP_info_t *info,
           QOP_FermionLinksDW *fldw,
                         REAL M5,
                         REAL mq,
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
  QOPPC(dw_dslash2_qdp)(info, fldw, M5, mq, qo, qi, ls, eo_out, eo_in);
}

void QOPPC(dw_diaginv_qdp)( QOP_info_t *info,
                    QOP_FermionLinksDW *fldw,
                                  REAL M5,
                                  REAL mq,
                      QDP_DiracFermion **out,
                      QDP_DiracFermion **in,
                                   int ls,
                         QOP_evenodd_t eo)
{
  REAL M0 = M5-5;
  check_setup(fldw,ls);
  if (eo==QOP_EVEN || eo==QOP_EVENODD)
      QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, 1, QDP_even);
  if (eo==QOP_ODD || eo==QOP_EVENODD)
      QOPPC(Qxxinv)(fldw, out, in, M0, mq, ls, 1, QDP_odd);
}
