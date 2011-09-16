//#define DO_TRACE
#include <qop_internal.h>

#define dblstore_style(x) ((x)&1)

extern QOP_asqtad_t QOP_asqtad;

static int old_style=-1;
static int old_optnum=-1;

#define NTMPSUB 4
#define NTMP (3*NTMPSUB)
#define NVTMP 24
static int dslash_setup = 0;
static QDP_ColorVector *vtemp[NTMP][NVTMP];
static QDP_ColorVector *tin[NTMP];
#define tmpnum(eo,n) ((eo)+3*((n)-1))
#define tmpsub(eo,n) tin[tmpnum(eo,n)]

#define check_setup(fla) \
{ \
  if( (!dslash_setup) || (QOP_asqtad.optnum != old_optnum) ) { \
    reset_temps(); \
  } \
  if( fla->dblstored != dblstore_style(QOP_asqtad.style) ) { \
    double_store(fla); \
  } \
}

static void
free_temps(void)
{
  if(dslash_setup) {
    int i, j;

    for(i=0; i<NTMP; i++) {
      QDP_destroy_V(tin[i]);
    }

    for(i=0; i<NTMP; i++) {
      for(j=0; j<NVTMP; j++) {
	QDP_destroy_V(vtemp[i][j]);
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
    tin[i] = QDP_create_V();
  }

  for(i=0; i<NTMP; i++) {
    for(j=0; j<NVTMP; j++) {
      vtemp[i][j] = QDP_create_V();
    }
  }

  dslash_setup = 1;
  old_style = QOP_asqtad.style;
  old_optnum = QOP_asqtad.optnum;
}

static void
double_store(QOP_FermionLinksAsqtad *fla)
{
  int i, nl2 = fla->nlinks/2;

  TRACE;
  if(fla->dblstored) {
    for(i=0; i<nl2; i++) {
      QDP_destroy_M(fla->bcklinks[i]);
    }
    fla->dblstored = 0;
  }
  TRACE;

  if(dblstore_style(QOP_asqtad.style)) {
    for(i=0; i<nl2; i++) {
      fla->bcklinks[i] = QDP_create_M();
    }
    if(nl2 == 8) {
      for(i=0; i<4; i++) {
	fla->dbllinks[4*i] = fla->fwdlinks[2*i];
	fla->dbllinks[4*i+1] = fla->fwdlinks[2*i+1];
	fla->dbllinks[4*i+2] = fla->bcklinks[2*i];
	fla->dbllinks[4*i+3] = fla->bcklinks[2*i+1];
      }
      for(i=0; i<16; i++) {
	fla->shifts_dbl[i] = QOP_asqtad.shifts_dbl[i];
	fla->shiftdirs_dbl[i] = QOP_asqtad.shiftdirs_dbl[i];
      }
    } else {
      for(i=0; i<4; i++) {
	fla->dbllinks[2*i] = fla->fwdlinks[i];
	fla->dbllinks[2*i+1] = fla->bcklinks[i];
      }
      for(i=0; i<4; i++) {
	fla->shifts_dbl[2*i] = QDP_neighbor[i];
	fla->shifts_dbl[2*i+1] = QDP_neighbor[i];
	fla->shiftdirs_dbl[2*i] = QDP_forward;
	fla->shiftdirs_dbl[2*i+1] = QDP_backward;
      }
    }
    TRACE;
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<nl2; i++) {
      QDP_M_eq_sM(m, fla->fwdlinks[i], fla->shifts[i],
                  QDP_backward, QDP_all);
      QDP_M_eqm_Ma(fla->bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    fla->dblstored = dblstore_style(QOP_asqtad.style);
  }
  TRACE;
}

QDP_ColorVector *
QOPPC(asqtad_dslash_get_tmp)(QOP_FermionLinksAsqtad *fla,
                             QOP_evenodd_t eo, int n)
{
  check_setup(fla);
  if(fla->nlinks==8) n += 2;
  if(n>=1 && n<=NTMPSUB) return tmpsub(eo,n);
  else return NULL;
}


/********************/
/*  link functions  */
/********************/

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_convert_L_from_qdp)(QDP_ColorMatrix *fatlinks[],
				 QDP_ColorMatrix *longlinks[])
{
  int i, nl, nl2;
  QOPPC(FermionLinksAsqtad) *fla;

  ASQTAD_DSLASH_BEGIN;
  if(!QOP_asqtad.inited) QOP_asqtad_invert_init();

  TRACE;
  nl = 8;
  if(longlinks) nl = 16;
  nl2 = nl/2;

  TRACE;
  QOP_malloc(fla, QOPPC(FermionLinksAsqtad), 1);
  QOP_malloc(fla->fatlinks, QDPPC(ColorMatrix) *, 4);
  if(longlinks) { QOP_malloc(fla->longlinks, QDPPC(ColorMatrix) *, 4); }
  else fla->longlinks = NULL;
  QOP_malloc(fla->fwdlinks, QDPPC(ColorMatrix) *, nl2);
  QOP_malloc(fla->bcklinks, QDPPC(ColorMatrix) *, nl2);
  QOP_malloc(fla->dbllinks, QDPPC(ColorMatrix) *, nl);
  TRACE;

  fla->dblstored = 0;
  for(i=0; i<4; i++) {
    fla->fatlinks[i] = fatlinks[i];
    if(longlinks) fla->longlinks[i] = longlinks[i];
  }
  if(longlinks) {
    for(i=0; i<4; i++) {
      fla->fwdlinks[2*i] = fatlinks[i];
      fla->fwdlinks[2*i+1] = longlinks[i];
    }
    for(i=0; i<8; i++) {
      fla->shifts[i] = QOP_asqtad.shifts[i];
    }
  } else {
    for(i=0; i<4; i++) {
      fla->fwdlinks[i] = fatlinks[i];
      fla->shifts[i] = QDP_neighbor[i];
    }
  }
  TRACE;
  // scale links
  for(i=0; i<nl2; i++) {
    QLA_Real f = 0.5;
    QDP_M_eq_r_times_M(fla->fwdlinks[i], &f, fla->fwdlinks[i], QDP_all);
  }

  //AB Need to set to NULL explicitly if unused
  if(!longlinks) fla->longlinks=NULL;

  fla->nlinks = nl;
  TRACE;
  check_setup(fla);
  TRACE;
  fla->eigcg.u = NULL;

  TRACE;
  ASQTAD_DSLASH_END;
  return fla;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_qdp)(QDP_ColorMatrix *fatlinks[],
				QDP_ColorMatrix *longlinks[])
{
  QOPPC(FermionLinksAsqtad) *fla;
  QDP_ColorMatrix *fl[4], *ll[4];
  int i;

  ASQTAD_DSLASH_BEGIN;

  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    QDP_M_eq_M(fl[i], fatlinks[i], QDP_all);
    if(longlinks) {
      ll[i] = QDP_create_M();
      QDP_M_eq_M(ll[i], longlinks[i], QDP_all);
    }
  }

  if(longlinks) {
    fla = QOPPC(asqtad_convert_L_from_qdp)(fl, ll);
  } else {
    fla = QOPPC(asqtad_convert_L_from_qdp)(fl, NULL);
  }

  ASQTAD_DSLASH_END;
  return fla;
}

void
QOPPC(asqtad_load_L_from_qdp)(QOPPC(FermionLinksAsqtad) *fla,
			      QDP_ColorMatrix *fatlinks[],
			      QDP_ColorMatrix *longlinks[])
{
  int i;

  ASQTAD_DSLASH_BEGIN;

  // copy and scale links
  for(i=0; i<4; i++) {
    QLA_Real f = 0.5;
    QDP_M_eq_r_times_M(fla->fatlinks[i], &f, fatlinks[i], QDP_all);
    if(longlinks) {
      // need to create naik links if not done & update pointers
      QDP_M_eq_r_times_M(fla->longlinks[i], &f, longlinks[i], QDP_all);
    } else {
      // need to free naik links & update pointers
    }
  }
  fla->dblstored = 0;
  check_setup(fla);
  ASQTAD_DSLASH_END;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_convert_L_from_raw)(REAL *fatlinks[], REAL *longlinks[],
				 QOP_evenodd_t evenodd)
{
  QDP_ColorMatrix *fl[4], *ll[4], **lp;
  int i;

  ASQTAD_DSLASH_BEGIN;

  lp = NULL;
  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    QDP_insert_M(fl[i], (QLA_ColorMatrix *)fatlinks[i], QDP_all);
    if(longlinks) {
      ll[i] = QDP_create_M();
      QDP_insert_M(ll[i], (QLA_ColorMatrix *)longlinks[i], QDP_all);
      lp = ll;
    }
  }

  ASQTAD_DSLASH_END;
  return QOPPC(asqtad_convert_L_from_qdp)(fl, lp);
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_raw)(REAL *fatlinks[], REAL *longlinks[],
				QOP_evenodd_t evenodd)
{
  QDP_ColorMatrix *fl[4], *ll[4], **lp;
  int i;

  ASQTAD_DSLASH_BEGIN;

  lp = NULL;
  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    QOP_qdp_eq_raw(M, fl[i], fatlinks[i], evenodd);
    if(longlinks) {
      ll[i] = QDP_create_M();
      QOP_qdp_eq_raw(M, ll[i], longlinks[i], evenodd);
      lp = ll;
    }
  }

  ASQTAD_DSLASH_END;
  return QOPPC(asqtad_convert_L_from_qdp)(fl, lp);
}

void
QOPPC(asqtad_load_L_from_raw)(QOPPC(FermionLinksAsqtad) *fla,
			      REAL *fatlinks[], REAL *longlinks[],
			      QOP_evenodd_t evenodd)
{
  int i, nl = 8;

  ASQTAD_DSLASH_BEGIN;

  for(i=0; i<4; i++) {
    QOP_qdp_eq_raw(M, fla->fatlinks[i], fatlinks[i], evenodd);
    if(longlinks) {
      // need to create naik links if not done & update pointers
      QOP_qdp_eq_raw(M, fla->longlinks[i], longlinks[i], evenodd);
    } else {
      nl = 4;
      // need to free naik links & update pointers
    }
  }
  fla->dblstored = 0;
  // scale links
  for(i=0; i<nl; i++) {
    QLA_Real f = 0.5;
    QDP_M_eq_r_times_M(fla->fwdlinks[i], &f, fla->fwdlinks[i], QDP_all);
  }

  check_setup(fla);
  ASQTAD_DSLASH_END;
}

/* Computes the staple :
                 mu
              +-------+
        nu    |       |
              |       |
              X       X
  Where the mu link can be any su3_matrix. The result is saved in staple.
  if staple==NULL then the result is not saved.
  It also adds the computed staple to the fatlink[mu] with weight coef.
*/
static void
compute_gen_staple(QDP_ColorMatrix *staple, int mu, int nu,
		   QDP_ColorMatrix *link, double dcoef,
		   QDP_ColorMatrix **gauge, QDP_ColorMatrix **fl,
		   QDP_ColorMatrix *ts0, QDP_ColorMatrix *ts1,
		   QDP_ColorMatrix *tmat1, QDP_ColorMatrix *tmat2,
		   QDP_ColorMatrix *ts2)
{
  QLA_Real coef = dcoef;
  //QDP_ColorMatrix *ts0, *ts1;
  //QDP_ColorMatrix *tmat1, *tmat2;

  //ts0 = QDP_create_M();
  //ts1 = QDP_create_M();
  //tmat1 = QDP_create_M();
  //tmat2 = QDP_create_M();

  /* Upper staple */
  if(link!=gauge[mu])
    QDP_M_eq_sM(ts0, link, QDP_neighbor[nu], QDP_forward, QDP_all);
  //QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);

  if(staple!=NULL) {  /* Save the staple */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(staple, gauge[nu], tmat1, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    QDP_M_eq_M_times_M(tmat2, gauge[nu], tmat1, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, tmat2, QDP_all);
  }

  /* lower staple */
  //QDP_M_eq_sM(ts1, gauge[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(tmat1, gauge[nu], link, QDP_all);
  QDP_M_eq_M_times_M(tmat2, tmat1, ts1, QDP_all);
  QDP_M_eq_sM(ts2, tmat2, QDP_neighbor[nu], QDP_backward, QDP_all);

  if(staple!=NULL) { /* Save the staple */
    QDP_M_peq_M(staple, ts2, QDP_all);
    QDP_M_peq_r_times_M(fl[mu], &coef, staple, QDP_all);
  } else {  /* No need to save the staple. Add it to the fatlinks */
    QDP_M_peq_r_times_M(fl[mu], &coef, ts2, QDP_all);
  }

  if(link!=gauge[mu]) QDP_discard_M(ts0);
  //QDP_discard_M(ts1);
  QDP_discard_M(ts2);
  //QDP_destroy_M(ts0);
  //QDP_destroy_M(ts1);
  //QDP_destroy_M(tmat1);
  //QDP_destroy_M(tmat2);
} /* compute_gen_staple */

static void
make_imp_links(QOP_info_t *info, QDP_ColorMatrix *fl[], QDP_ColorMatrix *ll[],
	       QOP_asqtad_coeffs_t *coeffs, QDP_ColorMatrix *gf[])
{
  QDP_ColorMatrix *staple=NULL, *tempmat1=NULL, *t1=NULL, *t2=NULL, *tsl[4], *tsg[4][4], *ts1[4], *ts2[4];
  double nflop = 0;
  double dtime;

#if 0
  {
#define p(s) QOP_printf0("%20s %g\n", #s, coeffs->s)
    p(one_link);
    p(three_staple);
    p(five_staple);
    p(seven_staple);
    p(lepage);
    p(naik);
#undef p
  }
#endif

  if(coeffs->three_staple || coeffs->lepage || coeffs->five_staple || coeffs->seven_staple || coeffs->naik) {
    nflop = 61632;
    staple = QDP_create_M();
    tempmat1 = QDP_create_M();
    t1 = QDP_create_M();
    t2 = QDP_create_M();
    for(int dir=0; dir<4; dir++) {
      tsl[dir] = QDP_create_M();
      ts1[dir] = QDP_create_M();
      ts2[dir] = QDP_create_M();
      for(int nu=0; nu<4; nu++) {
	tsg[dir][nu] = NULL;
	if(dir!=nu) {
	  tsg[dir][nu] = QDP_create_M();
	  QDP_M_eq_sM(tsg[dir][nu], gf[dir], QDP_neighbor[nu], QDP_forward, QDP_all);
	}
      }
    }
  }

  dtime = -QOP_time();

  /* to fix up the Lepage term, included by a trick below */
  QLA_Real one_link = coeffs->one_link - 6.0*coeffs->lepage;

  for(int dir=0; dir<4; dir++) {
    QDP_M_eq_r_times_M(fl[dir], &one_link, gf[dir], QDP_all);
    if(coeffs->three_staple || coeffs->lepage || coeffs->five_staple || coeffs->seven_staple) {
      for(int nu=0; nu<4; nu++) if(nu!=dir) {
	  if(coeffs->three_staple) {
	    compute_gen_staple(staple, dir, nu, gf[dir], coeffs->three_staple,
			       gf, fl, tsg[dir][nu], tsg[nu][dir], t1, t2, ts2[nu]);
	  }
	  if(coeffs->lepage) {
	    compute_gen_staple(NULL, dir, nu, staple, coeffs->lepage, gf, fl,
			       tsl[nu], tsg[nu][dir], t1, t2, ts2[nu]);
	  }
	  if(coeffs->five_staple || coeffs->seven_staple) {
	    for(int rho=0; rho<4; rho++) if((rho!=dir)&&(rho!=nu)) {
		compute_gen_staple(tempmat1, dir, rho, staple, coeffs->five_staple,
				   gf, fl, tsl[rho], tsg[rho][dir], t1, t2, ts2[rho]);
		if(coeffs->seven_staple) {
		  for(int sig=0; sig<4; sig++) {
		    if((sig!=dir)&&(sig!=nu)&&(sig!=rho)) {
		      compute_gen_staple(NULL, dir, sig, tempmat1, coeffs->seven_staple,
					 gf, fl, ts1[sig], tsg[sig][dir],t1,t2,ts2[sig]);
		    }
		  } /* sig */
		}
	      } /* rho */
	  }
	} /* nu */
    }
  } /* dir */

  /* long links */
  if(ll) {
    QLA_Real naik = coeffs->naik;
    for(int dir=0; dir<4; dir++) {
      QDP_M_eq_sM(staple, gf[dir], QDP_neighbor[dir], QDP_forward, QDP_all);
      QDP_M_eq_M_times_M(tempmat1, gf[dir], staple, QDP_all);
      QDP_M_eq_sM(ts1[dir], tempmat1, QDP_neighbor[dir], QDP_forward, QDP_all);
      QDP_M_eq_M_times_M(ll[dir], gf[dir], ts1[dir], QDP_all);
      QDP_M_eq_r_times_M(ll[dir], &naik, ll[dir], QDP_all);
    }
  }

  if(coeffs->three_staple || coeffs->lepage || coeffs->five_staple || coeffs->seven_staple || coeffs->naik) {
    QDP_destroy_M(staple);
    QDP_destroy_M(tempmat1);
    QDP_destroy_M(t1);
    QDP_destroy_M(t2);
    for(int dir=0; dir<4; dir++) {
      QDP_destroy_M(tsl[dir]);
      QDP_destroy_M(ts1[dir]);
      QDP_destroy_M(ts2[dir]);
      for(int nu=0; nu<4; nu++) {
	//if(dir!=nu) QDP_destroy_M(tsg[dir][nu]);
	if(tsg[dir][nu]!=NULL) QDP_destroy_M(tsg[dir][nu]);
      }
    }
  }

  dtime += QOP_time();
  //node0_printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,
  //       (Real)nflop*volume/(1e6*dtime*numnodes()) );

  info->final_sec = dtime;
  info->final_flop = nflop*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_G)(QOP_info_t *info,
			      QOP_asqtad_coeffs_t *coeffs,
			      QOPPC(GaugeField) *gauge)
{
  QOPPC(FermionLinksAsqtad) *fla;
  QDP_ColorMatrix *fl[4], *ll[4], **lp;
  int i;

  ASQTAD_DSLASH_BEGIN;

  TRACE;
  lp = ll;
  if(coeffs->naik==0) lp = NULL;
  for(i=0; i<4; i++) {
    fl[i] = QDP_create_M();
    if(lp) ll[i] = QDP_create_M();
  }
  TRACE;
  make_imp_links(info, fl, lp, coeffs, gauge->links);
  TRACE;
  fla = QOPPC(asqtad_convert_L_from_qdp)(fl, lp);
  TRACE;

  ASQTAD_DSLASH_END;
  return fla;
}

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_r_times_L)(QLA_Real *s,
			      QOPPC(FermionLinksAsqtad) *fla_src)
{
  QOPPC(FermionLinksAsqtad) *fla;
  QDP_ColorMatrix **fl_src = fla_src->fatlinks;
  QDP_ColorMatrix **ll_src = fla_src->longlinks;
  QDP_ColorMatrix *fl_dst[4], *ll_dst[4];
  int i;

  ASQTAD_DSLASH_BEGIN;

  /* Copy and multiply.  Correct for the 1/2 in source */

  for(i=0; i<4; i++) {
    QLA_Real two_s = 2.0 * (*s);
    fl_dst[i] = QDP_create_M();
    QDP_M_eq_r_times_M(fl_dst[i], &two_s, fl_src[i], QDP_all);
    if(ll_src) {
      ll_dst[i] = QDP_create_M();
      QDP_M_eq_r_times_M(ll_dst[i], &two_s, ll_src[i], QDP_all);
    }
  }

  /* Hand over the QDP fields to the links structure */

  if(ll_src)
    fla = QOPPC(asqtad_convert_L_from_qdp)(fl_dst, ll_dst);
  else
    fla = QOPPC(asqtad_convert_L_from_qdp)(fl_dst, NULL);


  ASQTAD_DSLASH_END;
  return fla;
}

/* Copy constructor */

QOPPC(FermionLinksAsqtad) *
QOPPC(asqtad_create_L_from_L)(QOPPC(FermionLinksAsqtad) *fla_src)
{
  QOPPC(FermionLinksAsqtad) *fla;
  QLA_Real one = 1.0;

  ASQTAD_DSLASH_BEGIN;

  /* Copy the source links */
  fla = QOPPC(asqtad_create_L_from_r_times_L)(&one, fla_src);

  ASQTAD_DSLASH_END;
  return fla;
}

void
QOPPC(asqtad_L_peq_L)(QOPPC(FermionLinksAsqtad) *fla,
		      QOPPC(FermionLinksAsqtad) *fla1)
{
  int i;

  ASQTAD_DSLASH_BEGIN;

  TRACE;
  for(i = 0; i < 4; i++){
    QDP_M_peq_M(fla->fatlinks[i], fla1->fatlinks[i], QDP_all);
    if(fla->longlinks && fla1->longlinks)
      QDP_M_peq_M(fla->longlinks[i], fla1->longlinks[i], QDP_all);
  }

  TRACE;

  double_store(fla);

  ASQTAD_DSLASH_END;

}

void
QOPPC(asqtad_load_L_from_G)(QOP_info_t *info,
			    QOPPC(FermionLinksAsqtad) *fla,
			    QOP_asqtad_coeffs_t *coeffs,
			    QOPPC(GaugeField) *gauge)
{
  int i, nl;

  ASQTAD_DSLASH_BEGIN;

  nl = fla->nlinks;
  make_imp_links(info, fla->fatlinks, fla->longlinks, coeffs, gauge->links);
  fla->dblstored = 0;
  // scale links
  for(i=0; i<nl/2; i++) {
    QLA_Real f = 0.5;
    QDP_M_eq_r_times_M(fla->fwdlinks[i], &f, fla->fwdlinks[i], QDP_all);
  }
  check_setup(fla);

  ASQTAD_DSLASH_END;
}

void
QOPPC(asqtad_extract_L_to_qdp)(QDP_ColorMatrix *fatlinks[],
			       QDP_ColorMatrix *longlinks[],
			       QOP_FermionLinksAsqtad *src)
{
  int i;

  ASQTAD_DSLASH_BEGIN;

  // copy and unscale links
  for(i=0; i<4; i++) {
    QLA_Real f = 2.0;
    QDP_M_eq_r_times_M(fatlinks[i], &f, src->fatlinks[i], QDP_all);
    // need to check if src->longlinks exists
    if(longlinks) QDP_M_eq_r_times_M(longlinks[i], &f, src->longlinks[i], QDP_all);
  }

  ASQTAD_DSLASH_END;
}

void
QOPPC(asqtad_extract_L_to_raw)(REAL *fatlinks[], REAL *longlinks[],
			       QOP_FermionLinksAsqtad *fla,
			       QOP_evenodd_t evenodd)
{
  int i, nl;

  ASQTAD_DSLASH_BEGIN;

  nl = fla->nlinks;
  // unscale links
  for(i=0; i<nl/2; i++) {
    QLA_Real f = 2.0;
    QDP_M_eq_r_times_M(fla->fwdlinks[i], &f, fla->fwdlinks[i], QDP_all);
  }
  for(i=0; i<4; i++) {
    QOP_raw_eq_qdp(M, fatlinks[i], fla->fatlinks[i], evenodd);
    // need to check if src->longlinks exists
    if(longlinks) { QOP_raw_eq_qdp(M, longlinks[i], fla->longlinks[i], evenodd); }
  }
  // scale links
  for(i=0; i<nl/2; i++) {
    QLA_Real f = 0.5;
    QDP_M_eq_r_times_M(fla->fwdlinks[i], &f, fla->fwdlinks[i], QDP_all);
  }

  ASQTAD_DSLASH_END;
}

void
QOPPC(asqtad_destroy_L)(QOPPC(FermionLinksAsqtad) *fla)
{
  int i, nl;

  ASQTAD_DSLASH_BEGIN;

  if(fla == NULL)return;

  nl = fla->nlinks;
  for(i=0; i<nl/2; i++) {
    QDP_destroy_M(fla->fwdlinks[i]);
  }
  if(fla->dblstored) {
    for(i=0; i<nl/2; i++) {
      QDP_destroy_M(fla->bcklinks[i]);
    }
  }
  free(fla->fatlinks);
  if(fla->longlinks) free(fla->longlinks);
  free(fla->fwdlinks);
  free(fla->bcklinks);
  free(fla->dbllinks);
  if(fla->eigcg.u) {
    for(i=0; i<fla->eigcg.numax; i++) {
      QDP_destroy_V(fla->eigcg.u[i]);
    }
    free(fla->eigcg.u);
    free(fla->eigcg.l);
  }
  free(fla);

  ASQTAD_DSLASH_END;
}


/********************/
/* Dslash functions */
/********************/

static void
asqtad_dslash0(QOP_FermionLinksAsqtad *fla,
	       QDP_ColorVector *dest, QDP_ColorVector *src,
	       QOP_evenodd_t eo, int n);

static void
asqtad_dslash1(QOP_FermionLinksAsqtad *fla,
	       QDP_ColorVector *dest, QDP_ColorVector *src,
	       QOP_evenodd_t eo, int n);

static void
asqtad_dslash0_dir(QOP_FermionLinksAsqtad *fla,
		   QDP_ColorVector *dest, QDP_ColorVector *src,
		   int dir, int fb, double wtfat, double wtlong,
		   QOP_evenodd_t eo, int n);

static void
asqtad_dslash1_dir(QOP_FermionLinksAsqtad *fla,
		   QDP_ColorVector *dest, QDP_ColorVector *src,
		   int dir, int fb, double wtfat, double wtlong,
		   QOP_evenodd_t eo, int n);

#define asqtad_hop(fla, dest, src, eo) \
{ \
  QDP_ColorVector *tsrc = src; \
  int _n = 1; \
  while(1) { \
    if(src==tmpsub(eo,_n)) break; \
    if(_n==NTMPSUB) { \
      _n = (fla->nlinks==16) ? 1 : 3;		\
      tsrc = tmpsub(eo,_n); \
      QDP_V_eq_V(tsrc, src, qdpsub(oppsub(eo))); \
      break; \
    } \
    _n++; \
  } \
  /*printf("%i %i\n", eo, _n);*/ \
  if(dblstore_style(QOP_asqtad.style)) { \
    asqtad_dslash1(fla, dest, tsrc, eo, _n); \
  } else { \
    asqtad_dslash0(fla, dest, tsrc, eo, _n); \
  } \
}

#define mass_term(fla, out, in, mass, sum, subset) { \
  QLA_Real qm = mass; \
  if(sum==NULL) { QDP_V_eq_r_times_V(out, &qm, in, subset); } \
  else { QDP_V_eq_r_times_V_plus_V(out, &qm, in, sum, subset); } \
}

void
QOPPC(asqtad_dslash)(QOP_info_t *info,
		     QOP_FermionLinksAsqtad *fla,
		     REAL mass,
		     QOP_ColorVector *out,
		     QOP_ColorVector *in,
		     QOP_evenodd_t eo_out,
		     QOP_evenodd_t eo_in)
{
  QOP_asqtad_dslash_qdp(info,fla,mass,out->cv,in->cv,eo_out,eo_in);
}

void
QOPPC(asqtad_diaginv)(QOP_info_t *info,
		      QOP_FermionLinksAsqtad *fla,
		      REAL mass,
		      QOP_ColorVector *out,
		      QOP_ColorVector *in,
		      QOP_evenodd_t eo)
{
  QOP_asqtad_diaginv_qdp(info,fla,mass,out->cv,in->cv,eo);
}

void
QOPPC(asqtad_dslash_qdp)(QOP_info_t *info,
			 QOP_FermionLinksAsqtad *fla,
			 REAL mass,
			 QDP_ColorVector *out,
			 QDP_ColorVector *in,
			 QOP_evenodd_t eo_out,
			 QOP_evenodd_t eo_in)
{
  check_setup(fla);

  if(eo_in==eo_out) {
    if(eo_out==QOP_EVENODD) {
      asqtad_hop(fla, out, in, QOP_EVENODD);
      mass_term(fla, out, in, mass, out, QDP_all);
    } else if(eo_out==QOP_EVEN) {
      mass_term(fla, out, in, mass, NULL, QDP_even);
    } else {
      mass_term(fla, out, in, mass, NULL, QDP_odd);
    }
  } else {
    if(eo_out==QOP_EVEN || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_ODD) {
        asqtad_hop(fla, out, in, QOP_EVEN);
      } else if(eo_in==QOP_EVEN) {
        mass_term(fla, out, in, mass, NULL, QDP_even);
      } else {
        asqtad_hop(fla, out, in, QOP_EVEN);
        mass_term(fla, out, in, mass, out, QDP_even);
      }
    }
    if(eo_out==QOP_ODD || eo_out==QOP_EVENODD) {
      if(eo_in==QOP_EVEN) {
        asqtad_hop(fla, out, in, QOP_ODD);
      } else if(eo_in==QOP_ODD) {
        mass_term(fla, out, in, mass, NULL, QDP_odd);
      } else {
        asqtad_hop(fla, out, in, QOP_ODD);
        mass_term(fla, out, in, mass, out, QDP_odd);
      }
    }
  }
}

void
QOPPC(asqtad_diaginv_qdp)(QOP_info_t *info,
			  QOP_FermionLinksAsqtad *asqtad,
			  REAL mass,
			  QDP_ColorVector *out,
			  QDP_ColorVector *in,
			  QOP_evenodd_t eo)
{
  QLA_Real f = 1.0/mass;
  QDP_V_eq_r_times_V(out, &f, in, qdpsub(eo));
}

/* shift (parallel transport a la dslash)

   (Used to construct vector interpolating operators.  Same as dslash,
   except always acts on all sites and works in only one direction.)
*/

#define asqtad_hop_dir(fla, dest, src, dir, fb, wtfat, wtlong, eo)	\
{ \
  QDP_ColorVector *tsrc = src; \
  int _n = 1; \
  while(1) { \
    if(src==tmpsub(eo,_n)) break; \
    if(_n==NTMPSUB) { \
      _n = (fla->nlinks==16) ? 1 : 3;		\
      tsrc = tmpsub(eo,_n); \
      QDP_V_eq_V(tsrc, src, qdpsub(oppsub(eo))); \
      break; \
    } \
    _n++; \
  } \
  /*printf("%i %i\n", eo, _n);*/ \
  if(dblstore_style(QOP_asqtad.style)) { \
    asqtad_dslash1_dir(fla, dest, tsrc, dir, fb, wtfat, wtlong, eo, _n);	\
  } else { \
    asqtad_dslash0_dir(fla, dest, tsrc, dir, fb, wtfat, wtlong, eo, _n);	\
  } \
}

void
QOPPC(asqtad_dslash_dir_qdp)(QOP_info_t *info,
			     QOP_FermionLinksAsqtad *fla,
			     int dir, int fb,
			     double wtfat, double wtlong,
			     QDP_ColorVector *out,
			     QDP_ColorVector *in,
			     QOP_evenodd_t eo_out)
{
  check_setup(fla);

  asqtad_hop_dir(fla, out, in, dir, fb, wtfat, wtlong, eo_out);
}


void
QOPPC(asqtad_dslash_dir)(QOP_info_t *info,
			 QOP_FermionLinksAsqtad *fla,
			 int dir, int fb,
			 double wtfat, double wtlong,
			 QOP_ColorVector *out,
			 QOP_ColorVector *in,
			 QOP_evenodd_t eo_out)
{
  QOP_asqtad_dslash_dir_qdp(info,fla,dir,fb,wtfat,wtlong,out->cv,in->cv,eo_out);
}

/* internal functions */

static void
asqtad_dslash0(QOP_FermionLinksAsqtad *fla,
	       QDP_ColorVector *dest, QDP_ColorVector *src,
	       QOP_evenodd_t eo, int n)
{
  int i, ntmp, nsv, nv, nl2=fla->nlinks/2;
  QDP_ColorVector *vsrc[nl2], *vdest[nl2], **temp;
  QDP_Subset dsubset, ssubset;
  dsubset = qdpsub(eo);
  ssubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);
  temp = vtemp[ntmp];

  for(i=0; i<nl2; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }
  nsv = QOP_asqtad.nsvec;
  if(nsv>nl2) nsv = nl2;
  nv = QOP_asqtad.nvec;
  if(nv>nl2) nv = nl2;

  /* Start gathers from positive directions */
  for(i=0; i<nl2; i+=nsv) {
    QDP_V_veq_sV(temp+i, vsrc+i, fla->shifts+i, QOP_common.shiftfwd+i,
		 dsubset, nsv);
  }

  /* Multiply by adjoint matrix at other sites */
  /* Start gathers from negative directions */
  for(i=0; i<nl2; i+=nsv) {
    QDP_V_veq_Ma_times_V(temp+16+i, fla->fwdlinks+i, vsrc+i, ssubset,
			 nsv);
    QDP_V_veq_sV(temp+8+i, temp+16+i, fla->shifts+i,
		 QOP_common.shiftbck+i, dsubset, nsv);
  }

  /* Multiply by matrix for positive directions and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<nl2; i+=nv) {
    QDP_V_vpeq_M_times_V(vdest+i, fla->fwdlinks+i, temp+i, dsubset, nv);
  }
  for(i=0; i<nl2; i++) QDP_discard_V(temp[i]);

  /* Wait gathers from negative directions, accumulate (negative) */
  for(i=0; i<nl2; i+=nv) {
    QDP_V_vmeq_V(vdest+i, temp+8+i, dsubset, nv);
  }
  for(i=0; i<nl2; ++i) QDP_discard_V(temp[i+8]);
}

static void
asqtad_dslash1(QOP_FermionLinksAsqtad *fla,
	       QDP_ColorVector *dest, QDP_ColorVector *src,
	       QOP_evenodd_t eo, int n)
{
  int i, ntmp, nsv, nv, nl=fla->nlinks;
  QDP_ColorVector *vsrc[nl], *vdest[nl], **temp;
  QDP_Subset dsubset, ssubset;
  dsubset = qdpsub(eo);
  ssubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);
  temp = vtemp[ntmp];

  for(i=0; i<nl; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }
  nsv = QOP_asqtad.nsvec;
  if(nsv>nl) nsv = nl;
  nv = QOP_asqtad.nvec;
  if(nv>nl) nv = nl;

  /* Start gathers from all directions */
  for(i=0; i<nl; i+=nsv) {
    QDP_V_veq_sV(temp+i, vsrc+i, fla->shifts_dbl+i,
		 fla->shiftdirs_dbl+i, dsubset, nsv);
  }

  /* Wait gathers from all directions, multiply by matrix and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<nl; i+=nv) {
    QDP_V_vpeq_M_times_V(vdest+i, fla->dbllinks+i, temp+i, dsubset, nv);
  }
  for(i=0; i<nl; ++i) QDP_discard_V(temp[i]);
}

static void
asqtad_dslash0_dir(QOP_FermionLinksAsqtad *fla,
		   QDP_ColorVector *dest, QDP_ColorVector *src,
		   int dir, int fb, double wtfat, double wtlong,
		   QOP_evenodd_t eo, int n)
{
  int i, ntmp, nsv, nv, nl2=fla->nlinks/2;
  QDP_ColorVector *vsrc[nl2], *vdest[nl2], **temp;
  QDP_Subset dsubset, ssubset;
  QLA_Real r[2];

  r[0] = wtfat;
  r[1] = wtlong;
  dsubset = qdpsub(eo);
  ssubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);
  temp = vtemp[ntmp];

  for(i=0; i<nl2; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  nsv = nv = 2;
  i = 2*dir;

  if(fb > 0){

    /* Shift from forward direction */

    /* Gather from positive directions and apply weights */
    QDP_V_veq_sV(temp+i, vsrc+i, fla->shifts+i, QOP_common.shiftfwd+i,
		 dsubset, nsv);

    QDP_V_veq_r_times_V(temp+i, r, temp+i, dsubset, nsv);

    /* Multiply by matrix for positive directions and accumulate */
    QDP_V_vpeq_M_times_V(vdest+i, fla->fwdlinks+i, temp+i, dsubset, nv);
    
  } else {

    /* Shift from backward direction */

    /* Multiply by adjoint matrix at other sites */
    /* Gather from negative directions, apply weights and accumulate */
  
    QDP_V_veq_Ma_times_V(temp+16+i, fla->fwdlinks+i, vsrc+i, ssubset,
			 nsv);
    QDP_V_veq_sV(temp+8+i, temp+16+i, fla->shifts+i,
		 QOP_common.shiftbck+i, dsubset, nsv);

    QDP_V_veq_r_times_V(temp+8+i, r, temp+8+i, dsubset, nsv);

    QDP_V_vpeq_V(vdest+i, temp+8+i, dsubset, nv);

  }

  for(i=0; i<nl2; i++){
    QDP_discard_V(temp[i]);
    QDP_discard_V(temp[i+nl2]);
  }
}

static void
asqtad_dslash1_dir(QOP_FermionLinksAsqtad *fla,
		   QDP_ColorVector *dest, QDP_ColorVector *src,
		   int dir, int fb, double wtfat, double wtlong,
		   QOP_evenodd_t eo, int n)
{
  int i, ntmp, nsv, nv, nl=fla->nlinks;
  QDP_ColorVector *vsrc[nl], *vdest[nl], **temp;
  QDP_Subset dsubset, ssubset;
  QLA_Real r[2];

  r[0] = wtfat;
  r[1] = wtlong;
  dsubset = qdpsub(eo);
  ssubset = qdpsub(oppsub(eo));
  ntmp = tmpnum(eo,n);
  temp = vtemp[ntmp];

  for(i=0; i<nl; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  nsv = nv = 2;

  if(fb > 0)
    i = 4*dir;
  else
    i = 4*dir + 2;

  /* Start gathers */
  QDP_V_veq_sV(temp+i, vsrc+i, fla->shifts_dbl+i,
	       fla->shiftdirs_dbl+i, dsubset, nsv);
  
  /* Wait gathers from all directions, apply weights, multiply by
     matrix and accumulate */
  QDP_V_veq_r_times_V(temp+i, r, temp+i, dsubset, nsv);

  QDP_V_vpeq_M_times_V(vdest+i, fla->dbllinks+i, temp+i, dsubset, nv);

  for(i=0; i<nl; ++i) QDP_discard_V(temp[i]);

}

/* rephase */

static int nd;
static int bc_dir;
static int bc_coord;
static int *bc_r0;
static QLA_Complex bc_phase;
static int staggered_sign_bits;

static void
set_staggered_phases(QLA_ColorMatrix *m, int coords[]){
  if(staggered_sign_bits) {
    int s=0;
    for(int i=0; i<nd; i++){
      int n = QDP_coord_size(i);
      int rshift = (coords[i] + n - bc_r0[i]) % n;
      s += ((staggered_sign_bits>>i)&1) * rshift;
    }
    if(s&1) {
      QLA_M_eqm_M(m, m);
    }
  }
}

static void
rephase_fat_bdry_func(QLA_ColorMatrix *m, int coords[])
{
  if(bc_dir>=0) {
    int rshift = (coords[bc_dir] + bc_coord - bc_r0[bc_dir]) % bc_coord;
    if(rshift == bc_coord-1) {
      QLA_ColorMatrix t;
      QLA_M_eq_c_times_M(&t, &bc_phase, m);
      QLA_M_eq_M(m, &t);
    }
  }

  set_staggered_phases(m, coords);
}

static void
rephase_long_bdry_func(QLA_ColorMatrix *m, int coords[])
{
  if(bc_dir>=0) {
    /* Apply the phase to three long links from nd-3 to nd-1 */
    int rshift = (coords[bc_dir] + bc_coord - bc_r0[bc_dir]) % bc_coord;
    if(rshift >= bc_coord-3) {
      QLA_ColorMatrix t;
      QLA_M_eq_c_times_M(&t, &bc_phase, m);
      QLA_M_eq_M(m, &t);
    }
  }

  set_staggered_phases(m, coords);
}

void
QOPPC(asqtad_rephase_L)(QOPPC(FermionLinksAsqtad) *fla,
	      int *r0,
	      QOP_bc_t *bc,
	      QOP_staggered_sign_t *sign)
{
  int nl2 = fla->nlinks/2;
  bc_r0 = r0;
  nd = QDP_ndim();
  for(int i=0; i<nd; i++) {
    bc_dir = -1;
    staggered_sign_bits = 0;
    if(bc && (bc->phase[i].re!=1. || bc->phase[i].im!=0.)) {
      bc_dir = i;
      bc_coord = QDP_coord_size(i);
      QLA_c_eq_r_plus_ir(bc_phase, bc->phase[i].re, bc->phase[i].im);
    }
    if(sign) {
      staggered_sign_bits = sign->signmask[i];
    }
    if(sign || bc_dir >= 0){
      if(nl2 == 8) {
	/* Forward fat links */
	QDP_M_eq_func(fla->fwdlinks[2*i], rephase_fat_bdry_func, QDP_all);
	/* Forward long links */
	QDP_M_eq_func(fla->fwdlinks[2*i+1], rephase_long_bdry_func, QDP_all);
      } else {
	/* Forward fat links only in this case */
	QDP_M_eq_func(fla->fwdlinks[i], rephase_fat_bdry_func, QDP_all);
      }
    }
  }

  /* Refresh double links if we are using them */
  double_store(fla);
}
