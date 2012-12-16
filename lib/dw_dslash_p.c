//#define DO_TRACE
#include <math.h>
#include <qop_internal.h>

#define dblstore_style(x) ((x)&1)
#define shiftd_style(x) ((x)&2)

extern int QOP_dw_initQ;
extern int QOP_dw_style;

static int old_style=-1;

// -----------------------------------------------------------------
// Temporary variable management

static int dslash_initQ = 0; // Is dslash set up?
// The temporaries are born here
static QDP_HalfFermion **thv0=NULL, **thv1=NULL, *tmph=NULL; // Used in Qxxinv
static QDP_DiracFermion *tmpd=NULL, // Used in Qxx
                        **tdv=NULL, // Used in schur and reconstruct
                        **tdv2=NULL; // Used in dslash2, schur2 and reconstruct

#define check_setup(ls) \
{ \
  if( (!dslash_initQ) ) { \
    reset_temps(NCARG ls);  \
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
      QDP_destroy_D(tdv2[s]);
    }
    free(thv0);
    free(thv1);
    free(tdv);
    free(tdv2);
  }
  dslash_initQ = 0;
}

#define NC nc
static void
reset_temps(NCPROT int ls)
{
  free_temps(ls);
  
  tmph = QDP_create_H();
  tmpd = QDP_create_D();
  
  thv0 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  thv1 = (QDP_HalfFermion**) malloc(ls*sizeof(QDP_HalfFermion*));
  tdv  = (QDP_DiracFermion**) malloc(ls*sizeof(QDP_DiracFermion*));
  tdv2 = (QDP_DiracFermion**) malloc(ls*sizeof(QDP_DiracFermion*));
  for (int s=0; s<ls; s++) {
    thv0[s] = QDP_create_H();
    thv1[s] = QDP_create_H();
    tdv[s]  = QDP_create_D();
    tdv2[s]  = QDP_create_D();
  }

  dslash_initQ = 1;
  old_style = QOP_dw_style;
#undef NC
}

// -----------------------------------------------------------------
// Gauge field management

#define NC int nc
QOP_FermionLinksDW *
QOP_dw_create_L_from_raw(REAL *links[], QOP_evenodd_t evenodd )
#undef NC
{
#define NC nc
  QOP_FermionLinksDW *fldw;
  QOP_malloc(fldw, QOP_FermionLinksDW, 1);
  REAL noclov = 0.0; // No clover term

  fldw->flw = QOP_wilson_create_L_from_raw(links, &noclov, evenodd);

  return fldw;
#undef NC
}

QOP_FermionLinksDW *
QOP_dw_create_L_from_G( QOP_info_t *info, 
                   QOP_dw_coeffs_t *coeffs,
                    QOP_GaugeField *gauge )
{ 
  QOP_FermionLinksDW *fldw;
  QOP_malloc(fldw, QOP_FermionLinksDW, 1);

  QOP_wilson_coeffs_t wcoeffs;
  wcoeffs.clov_s = 0.0;
  wcoeffs.clov_t = 0.0;
  wcoeffs.aniso = 1.0;

  fldw->flw = QOP_wilson_create_L_from_G(info, &wcoeffs, gauge);

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
  QOP_wilson_destroy_L(fldw->flw);
  free(fldw);
}

#define NC int nc
QOP_FermionLinksDW *
QOP_dw_convert_L_from_raw( REAL *links[],
                  QOP_evenodd_t evenodd )
#undef NC
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
  QOP_malloc(fldw, QOP_FermionLinksDW, 1);

  fldw->flw = QOP_wilson_create_L_from_qdp(links, NULL);

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
  QOP_malloc(fldw, QOP_FermionLinksDW, 1);

  fldw->flw = QOP_wilson_convert_L_from_qdp(links, NULL);

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
// Note that these operators are functions of M0 = 5-M5, not M5.
// They will not be accessible to general users, who should use
// dslash, diaginv, schur, EO_project and EO_reconstruct instead.
static void
QOP_Qxy(QOP_FermionLinksDW *fldw,
	QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	QLA_Real M0, QLA_Real mq, int ls, int sign,
	QDP_Subset osubset, QDP_Subset subset)
{
#define NC QDP_get_nc(out[0])
  check_setup(ls);

  for (int s=0; s<ls; s++) {
    QOP_wilson_dslash_qdp(NULL, fldw->flw,
                          1.0/(2*M0+8) /*not used*/, sign,
			  out[s], in[s],
			  (osubset==QDP_even?QOP_EVEN:QOP_ODD),
			  ( subset==QDP_even?QOP_EVEN:QOP_ODD));
  }
#undef NC
}

static void
QOP_Qxx(QOP_FermionLinksDW *fldw,
	QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	QLA_Real M0, QLA_Real mq, int ls, int sign,
	QDP_Subset subset, int add)
{
#define NC QDP_get_nc(out[0])
  check_setup(ls);
  for (int s=0; s<ls; s++) {
    if (add) {
      QDP_D_peq_r_times_D(out[s], &M0, in[s], subset);
    } else {
      QDP_D_eq_r_times_D(out[s], &M0, in[s], subset);
    }
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
  }
#undef NC
}

static void
QOP_Qxxinv(QOP_FermionLinksDW *fldw,
	   QDP_DiracFermion *out[], QDP_DiracFermion *in[],
	   QLA_Real M0, QLA_Real mq, int ls, int sign, QDP_Subset subset)
{
#define NC QDP_get_nc(out[0])
  QLA_Real fac, den, mba;
  QDP_DiracFermion *tin[2];
  QDP_HalfFermion *thf[2];
  int dirs[2]={4,4}, sgns[2]={-sign,sign};
  int s;
  
  check_setup(ls);

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
#undef NC
}

/* For the standard case, where eo = e, project =
  [ 1   -Qeo 1/Qoo ]
  [ 0       1      ]
*/
void
QOP_dw_EO_project( QOP_FermionLinksDW *fldw,
		   QDP_DiracFermion *out[],
		   QDP_DiracFermion *in[],
		   QLA_Real M5,
		   QLA_Real mq,
		   int ls,
		   QOP_evenodd_t eo )
{
#define NC QDP_get_nc(out[0])
  check_setup(ls);
  QLA_Real M0 = 5-M5;
  QDP_Subset  eosub = (eo==QOP_EVEN?QDP_even:QDP_odd),
             oppsub = (eo==QOP_EVEN?QDP_odd:QDP_even);
  QOP_Qxxinv(fldw,  tdv,  in, M0, mq, ls, 1, oppsub);
  QOP_Qxy   (fldw, tdv2, tdv, M0, mq, ls, 1, eosub, oppsub);
  for (int s=0; s<ls; s++) QDP_D_eq_D_minus_D(out[s], in[s], tdv2[s], eosub);
  for (int s=0; s<ls; s++) QDP_D_eq_D(out[s], in[s], oppsub);
#undef NC
}

/* For the standard case, where eo = e, reconstruct =
  [       1/Qee          0   ]
  [ -1/Qoo Qoe 1/Qee   1/Qoo ]
*/
void
QOP_dw_EO_reconstruct( QOP_FermionLinksDW *fldw,
		       QDP_DiracFermion *out[],
		       QDP_DiracFermion *in[],
		       QLA_Real M5,
		       QLA_Real mq,
		       int ls,
		       QOP_evenodd_t eo )
{
#define NC QDP_get_nc(out[0])
  check_setup(ls);
  QLA_Real M0 = 5-M5;
  QDP_Subset  eosub = (eo==QOP_EVEN?QDP_even:QDP_odd),
             oppsub = (eo==QOP_EVEN?QDP_odd:QDP_even);
  QOP_Qxxinv(fldw, out, out, M0, mq, ls, 1, eosub);
  QOP_Qxy   (fldw, tdv, out, M0, mq, ls, 1, oppsub, eosub);
  for(int s=0; s<ls; s++) QDP_D_eq_D_minus_D(tdv[s], in[s], tdv[s], oppsub);
  QOP_Qxxinv(fldw, out, tdv, M0, mq, ls, 1, oppsub);

#if 0
  for (int s=0; s<ls; s++) QDP_D_eq_D(tdv[s], in[s], oppsub);
  QOP_Qxxinv(fldw, tdv,  in, M0, mq, ls, 1, eosub);
  QOP_Qxy   (fldw, out, tdv, M0, mq, ls, 1, oppsub, eosub);
  for (int s=0; s<ls; s++) QDP_D_eq_D_minus_D(tdv2[s], tdv[s], out[s], oppsub);
  QOP_Qxxinv(fldw, out, tdv2, M0, mq, ls, 1, oppsub);
  for (int s=0; s<ls; s++) QDP_D_eq_D(out[s], tdv[s], eosub);
#endif
#undef NC
}


// -----------------------------------------------------------------
// Domain-wall operators
// Note that in our convention, free-field M5 = 1
void
QOP_dw_dslash_qdp( QOP_info_t *info,
		   QOP_FermionLinksDW *fldw,
		   QLA_Real M5,
		   QLA_Real mq,
		   int sign,
		   QDP_DiracFermion *out[],
		   QDP_DiracFermion *in[],
		   int ls,
                   QOP_evenodd_t eo_out,
                   QOP_evenodd_t eo_in)
{
#define NC QDP_get_nc(out[0])
  QLA_Real M0 = 5-M5;
  check_setup(ls);
  
  if ( eo_out==QOP_EVEN || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_ODD ) {
      QOP_Qxy(fldw, out, in, M0, mq, ls, sign, QDP_even, QDP_odd); // Qeo
    } else if ( eo_in==QOP_EVEN ) {
      QOP_Qxx(fldw, out, in, M0, mq, ls, sign, QDP_even, 0);       // Qee
    } else {
      QOP_Qxy(fldw, out, in, M0, mq, ls, sign, QDP_even, QDP_odd); // Qeo
      QOP_Qxx(fldw, out, in, M0, mq, ls, sign, QDP_even, 1);       // Qee
    }
  }

  if ( eo_out==QOP_ODD || eo_out==QOP_EVENODD ) {
    if ( eo_in==QOP_EVEN ) {
      QOP_Qxy(fldw, out, in, M0, mq, ls, sign, QDP_odd, QDP_even); // Qoe
    } else if( eo_in==QOP_ODD ) {
      QOP_Qxx(fldw, out, in, M0, mq, ls, sign, QDP_odd, 0);        // Qoo
    } else {
      QOP_Qxy(fldw, out, in, M0, mq, ls, sign, QDP_odd, QDP_even); // Qoe
      QOP_Qxx(fldw, out, in, M0, mq, ls, sign, QDP_odd, 1);        // Qoo
    }
  }
#undef NC
}

void
QOP_dw_dslash( QOP_info_t *info,
	       QOP_FermionLinksDW *fldw,
	       QLA_Real M5,
	       QLA_Real mq,
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
  QOP_dw_dslash_qdp(info, fldw, M5, mq, sign, qo, qi, ls, eo_out, eo_in);
}

void
QOP_dw_dslash2_qdp( QOP_info_t *info,
		    QOP_FermionLinksDW *fldw,
		    QLA_Real M5,
		    QLA_Real mq,
		    QDP_DiracFermion *out[],
		    QDP_DiracFermion *in[],
		    int ls,
                    QOP_evenodd_t eo_out,
                    QOP_evenodd_t eo_in)
{
  QOP_dw_dslash_qdp(info, fldw, M5, mq,  1, tdv2, in, ls,
                       QOP_EVENODD, eo_in);
  QOP_dw_dslash_qdp(info, fldw, M5, mq, -1, out, tdv2, ls,
                       eo_out, QOP_EVENODD);
}

void
QOP_dw_dslash2( QOP_info_t *info,
		QOP_FermionLinksDW *fldw,
		QLA_Real M5,
		QLA_Real mq,
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
  QOP_dw_dslash2_qdp(info, fldw, M5, mq, qo, qi, ls, eo_out, eo_in);
}

void
QOP_dw_diaginv_qdp( QOP_info_t *info,
		    QOP_FermionLinksDW *fldw,
		    QLA_Real M5,
		    QLA_Real mq,
		    QDP_DiracFermion **out,
		    QDP_DiracFermion **in,
		    int ls,
                    QOP_evenodd_t eo)
{
  QLA_Real M0 = 5-M5;
  if (eo==QOP_EVEN || eo==QOP_EVENODD)
      QOP_Qxxinv(fldw, out, in, M0, mq, ls, 1, QDP_even);
  if (eo==QOP_ODD || eo==QOP_EVENODD)
      QOP_Qxxinv(fldw, out, in, M0, mq, ls, 1, QDP_odd);
}

/* For the standard case, where eo = e, schur =
  [ 1 - Qeo 1/Qoo Qoe 1/Qee   0 ]
  [           0               1 ]
*/
void
QOP_dw_schur_qdp( QOP_info_t *info,
		  QOP_FermionLinksDW *fldw,
		  QLA_Real M5,
		  QLA_Real mq,
		  int sign,
		  QDP_DiracFermion **out,
		  QDP_DiracFermion **in,
		  int ls,
                  QOP_evenodd_t eo)
{
  QLA_Real M0 = 5-M5;
  QDP_Subset  eosub = (eo==QOP_EVEN?QDP_even:QDP_odd),
             oppsub = (eo==QOP_EVEN?QDP_odd:QDP_even);
  if (sign>0) {
    QOP_Qxxinv(fldw, out,  in, M0, mq, ls, sign, eosub);
    QOP_Qxy   (fldw, tdv, out, M0, mq, ls, sign, oppsub, eosub);
    QOP_Qxxinv(fldw, out, tdv, M0, mq, ls, sign, oppsub);
    QOP_Qxy   (fldw, tdv, out, M0, mq, ls, sign, eosub, oppsub);
  } else {
    QOP_Qxy   (fldw, out,  in, M0, mq, ls, sign, oppsub, eosub);
    QOP_Qxxinv(fldw, tdv, out, M0, mq, ls, sign, oppsub);
    QOP_Qxy   (fldw, out, tdv, M0, mq, ls, sign, eosub, oppsub);
    QOP_Qxxinv(fldw, tdv, out, M0, mq, ls, sign, eosub);
  }
  for(int s=0; s<ls; s++)
    QDP_D_eq_D_minus_D(out[s], in[s], tdv[s], eosub);
}

void
QOP_dw_schur2_qdp( QOP_info_t *info,
		   QOP_FermionLinksDW *fldw,
		   QLA_Real M5,
		   QLA_Real mq,
		   QDP_DiracFermion *out[],
		   QDP_DiracFermion *in[],
		   int ls,
                   QOP_evenodd_t eo)
{
  QOP_dw_schur_qdp(info, fldw, M5, mq,  1, tdv2,   in, ls, eo);
  QOP_dw_schur_qdp(info, fldw, M5, mq, -1,  out, tdv2, ls, eo);
}
