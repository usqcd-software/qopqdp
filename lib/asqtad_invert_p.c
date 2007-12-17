/**************************************************************************
Asqtad inverter conventions:

in even-odd blocks (where (a,b) is either (e,o) or (o,e))

D = ( m     D_ab )
    ( D_ba  m    )

D^-1 = ( m A^-1      -D_ab B^-1 ) = ( m A^-1      -A^-1 D_ab             )
       ( -D_ba A^-1  m B^-1     )   ( -D_ba A^-1  [1 + D_ba A^-1 D_ab]/m )

with A = (m^2 - D_ab D_ba) and B = (m^2 - D_ba D_ab)

with even-odd preconditioning we can write the solution as

D^-1 (x) = ( A^-1 z                )
     (y)   ( y/m - D_ba A^-1 z / m )

with z = m x - D_ab y.
***************************************************************************/
#include <qop_internal.h>

extern QOP_asqtad_t QOP_asqtad;
extern QDP_ColorVector *QOPPC(asqtad_dslash_get_tmp)
     (QOP_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);

/* inverter stuff */

static QOP_FermionLinksAsqtad *gl_fla;
static REAL gl_mass;
static QOP_evenodd_t gl_eo;
static QDP_ColorVector *gl_tmp, *gl_tmp2;

/* calculates out = mass*in - D*in on subset eo */
static void
project(QOP_FermionLinksAsqtad *fla, QLA_Real mass, QDP_ColorVector *out,
	QDP_ColorVector *in, QOP_evenodd_t eo)
{
  //QOP_asqtad_diaginv_qdp(NULL, fla, mass, gl_tmp2, in, oppsub(eo));
  //QOP_asqtad_dslash_qdp(NULL, fla, mass, out, gl_tmp2, eo, oppsub(eo));
  //QDP_V_eq_V_minus_V(out, in, out, qdpsub(eo));
  QOP_asqtad_dslash_qdp(NULL, fla, mass, out, in, eo, oppsub(eo));
  QDP_V_eq_r_times_V_minus_V(out, &mass, in, out, qdpsub(eo));
}

/* calculates out = (src - D*soln)/mass on subset eo */
static void
reconstruct(QOP_FermionLinksAsqtad *fla, QLA_Real mass, QDP_ColorVector *out,
	    QDP_ColorVector *soln, QDP_ColorVector *src, QOP_evenodd_t eo)
{
  QDP_V_eq_V(gl_tmp, soln, qdpsub(oppsub(eo)));
  QOP_asqtad_dslash_qdp(NULL, fla, mass, out, gl_tmp, eo, oppsub(eo));
  QDP_V_eq_V_minus_V(out, src, out, qdpsub(eo));
  QOP_asqtad_diaginv_qdp(NULL, fla, mass, out, out, eo);
}

void
QOPPC(asqtad_invert_d2)(QDP_ColorVector *out, QDP_ColorVector *in,
			QDP_Subset subset)
{
  QLA_Real m2 = gl_mass*gl_mass;
  //QLA_Real m2 = -1.0/(gl_mass*gl_mass);
  QOP_asqtad_dslash_qdp(NULL,gl_fla,gl_mass, gl_tmp2, in, oppsub(gl_eo),gl_eo);
  QOP_asqtad_dslash_qdp(NULL,gl_fla,gl_mass, out,gl_tmp2, gl_eo,oppsub(gl_eo));
  QDP_V_eq_r_times_V_minus_V(out, &m2, in, out, qdpsub(gl_eo));
  //QDP_V_eq_r_times_V_plus_V(out, &m2, out, in, qdpsub(gl_eo));
}

void
QOP_asqtad_invert(QOP_info_t *info,
		  QOP_FermionLinksAsqtad *fla,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  REAL mass,
		  QOP_ColorVector *out,
		  QOP_ColorVector *in)
{
  /* cg has 5 * 12 = 60 flops/site/it */
  /* MdagM -> 2*(66+72*15)+12 = 2304 flops/site */
  double dtime = 0;
  double nflop = 0.5 * 2364;
  double rsqminold;
  QLA_Real rsq, rsqstop, insq;
  QDP_ColorVector *qdpin, *qdpout;
  QDP_ColorVector *cgp, *cgr;
  QDP_Subset insub, cgsub;
  QOP_evenodd_t ineo, cgeo;
  int iter = 0;
  int max_iter_old = inv_arg->max_iter;
  int max_restarts_old = inv_arg->max_restarts;
  int nrestart = -1, max_restarts = inv_arg->max_restarts;
  if(max_restarts<=0) max_restarts = 5;

  if(QOP_asqtad.cgtype==1) {
    fla->eigcg.numax = QOP_asqtad.eigcg_numax;
    fla->eigcg.m = QOP_asqtad.eigcg_m;
    fla->eigcg.nev = QOP_asqtad.eigcg_nev;
  }

  ineo = inv_arg->evenodd;
  insub = qdpsub(ineo);

  qdpin = QDP_create_V();
  qdpout = QDP_create_V();

  gl_fla = fla;
  gl_mass = mass;

  cgeo = ineo;
  if(ineo==QOP_EVENODD) cgeo = QOP_EVEN;
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  cgp = QOPPC(asqtad_dslash_get_tmp)(fla, oppsub(cgeo), 1);
  cgr = QOPPC(asqtad_dslash_get_tmp)(fla, oppsub(cgeo), 2);
  gl_tmp = cgr;
  gl_tmp2 = QOPPC(asqtad_dslash_get_tmp)(fla, cgeo, 1);

  QDP_V_eq_zero(qdpin, QDP_all);
  if(ineo==cgeo) {
    QLA_Real m = mass;
    //QDP_V_eq_V(qdpin, in->cv, insub);
    QDP_V_eq_r_times_V(qdpin, &m, in->cv, insub);
  } else {
    QDP_V_eq_zero(cgp, QDP_all);
    QDP_V_eq_V(gl_tmp2, in->cv, insub);
    project(fla, mass, qdpin, gl_tmp2, cgeo);
    QDP_V_eq_V(qdpin, in->cv, qdpsub(oppsub(ineo)));
  }
  //QOP_asqtad_diaginv_qdp(NULL, fla, mass, cgp, qdpin, cgeo);
  //QDP_V_eq_V(qdpin, cgp, cgsub);
  //if(ineo!=QOP_EVENODD && ineo!=cgeo) {
  //QDP_V_eq_zero(qdpout, cgsub);
  //reconstruct(fla, mass, qdpout, out->cv, qdpout, oppsub(ineo));
  //}
  QDP_V_eq_V(qdpout, out->cv, cgsub);

  QDP_r_eq_norm2_V(&insq, in->cv, insub);
  rsqstop = insq * res_arg->rsqmin;
  VERB(LOW, "ASQTAD_INVERT: rsqstop = %g\n", rsqstop);
  rsqminold = res_arg->rsqmin;
  res_arg->rsqmin *= 0.9;
  inv_arg->max_restarts = 0;
  do {
    inv_arg->max_iter = max_iter_old - iter;

    dtime -= QOP_time();

    if(QOP_asqtad.cgtype==1) {
      QOPPC(invert_eigcg_V)(QOPPC(asqtad_invert_d2), inv_arg, res_arg,
			    qdpout, qdpin, cgp, cgsub, &fla->eigcg);
    } else {
      QOPPC(invert_cg_V)(QOPPC(asqtad_invert_d2), inv_arg, res_arg,
			 qdpout, qdpin, cgp, cgsub);
    }
    //printf("resid = %g\n", res_arg->final_rsq);

    dtime += QOP_time();

    reconstruct(fla, mass, qdpout, qdpout, qdpin, oppsub(cgeo));
    //QOP_asqtad_dslash_qdp(NULL, fla, mass, cgr, qdpout, oppsub(ineo), QOP_EVENODD);
    //QDP_r_eq_norm2_V(&rsq, cgr, qdpsub(oppsub(ineo)));
    //printf("0 ?= %g\n", rsq);

    // get final residual
    //QDP_r_eq_norm2_D(&rsq, qdpout, QDP_all);
    //printf("nrm = %g\n", rsq);
    QOP_asqtad_dslash_qdp(NULL, fla, mass, cgr, qdpout, ineo, QOP_EVENODD);
    QDP_V_meq_V(cgr, in->cv, insub);
    QDP_r_eq_norm2_V(&rsq, cgr, insub);
    //printf("%i %i rsq = %g\tprec rsq = %g\trsqstop = %g\n", nrestart,
    //res_arg->final_iter, rsq, res_arg->final_rsq, rsqstop);
    res_arg->rsqmin *= 0.9*rsqstop/rsq;
    iter += res_arg->final_iter;
    VERB(LOW, "ASQTAD_INVERT: iter %i rsq = %g\n", iter, rsq);
  } while((rsq>rsqstop)&&(++nrestart<max_restarts));

  QDP_V_eq_V(out->cv, qdpout, insub);

  QDP_destroy_V(qdpin);
  QDP_destroy_V(qdpout);

  inv_arg->max_iter = max_iter_old;
  inv_arg->max_restarts = max_restarts_old;
  res_arg->rsqmin = rsqminold;
  res_arg->final_iter = iter;
  res_arg->final_rsq = rsq/insq;
  res_arg->final_restart = nrestart;

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}

void
QOP_asqtad_invert_multi(QOP_info_t *info,
			QOP_FermionLinksAsqtad *fla,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t **res_arg[],
			REAL *masses[], int nmass[],
			QOP_ColorVector **out_pt[],
			QOP_ColorVector *in_pt[],
			int nsrc)
{
  /* cg has 5 * 12 = 60 flops/site/it */
  /* MdagM -> 2*(66+72*15)+12 = 2304 flops/site */
  double dtime;
  double nflop = 0.5 * 2364;
  double nflopm = 0.5 * 30; /* per extra mass */
  QDP_ColorVector *cgp, *cgr;
  QDP_Subset insub, cgsub;
  QOP_evenodd_t ineo, cgeo;
  int i, j;

  ineo = inv_arg->evenodd;
  insub = qdpsub(ineo);

  gl_fla = fla;

  cgeo = ineo;
  if(ineo==QOP_EVENODD) {
    cgeo = QOP_EVEN;
    printf("warning: QOP_asqtad_invert_multi on ALL subset not supported yet!\n");
  }
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  cgp = QOPPC(asqtad_dslash_get_tmp)(fla, oppsub(cgeo), 1);
  cgr = QOPPC(asqtad_dslash_get_tmp)(fla, oppsub(cgeo), 2);
  gl_tmp = cgr;
  gl_tmp2 = QOPPC(asqtad_dslash_get_tmp)(fla, cgeo, 1);

  info->final_flop = 0;
  dtime = -QDP_time();

#if 0  // fake version
  for(i=0; i<nsrc; i++) {
    for(j=0; j<nmass[i]; j++) {
      gl_mass = masses[i][j];

      QOPPC(invert_cg_V)(QOPPC(asqtad_invert_d2), inv_arg, res_arg[i][j],
			 out_pt[i][j]->cv, in_pt[i]->cv, cgp, cgsub);

      info->final_flop += nflop*res_arg[i][j]->final_iter*QDP_sites_on_node;
    }
  }
#else // real multimass
  if( (nsrc==2) && (nmass[0]==1) && (nmass[1]==1) &&
      (masses[0][0]!=masses[1][0]) ) {  /* two source version */
    QLA_Real shifts[2], st, rsq0, rsq1;
    QDP_ColorVector *incv[2], *outcv[2], *x0, *src;
    int imin=0, i;
    QOP_resid_arg_t *ra[2];
    double rsqsave[2];

    x0 = QDP_create_V();
    src = QDP_create_V();

    for(i=0; i<2; i++) {
      shifts[i] = masses[i][0]*masses[i][0];
      if(shifts[i]<shifts[imin]) imin = i;
      incv[i] = in_pt[i]->cv;
      outcv[i] = out_pt[i][0]->cv;
      ra[i] = res_arg[i][0];
    }
    gl_mass = masses[imin][0];
    for(i=0; i<2; i++) shifts[i] -= gl_mass*gl_mass;

    st = 1/(shifts[1]-shifts[0]);
    QDP_V_eq_V_minus_V(x0, incv[1], incv[0], cgsub);
    QDP_V_eq_r_times_V(x0, &st, x0, cgsub);
    QDP_V_eq_V(cgp, x0, cgsub);
    QOPPC(asqtad_invert_d2)(src, cgp, cgsub);
    QDP_V_eq_V_minus_V(src, incv[imin], src, cgsub);

    QDP_r_eq_norm2_V(&rsq1, src, cgsub);
    for(i=0; i<2; i++) {
      QDP_r_eq_norm2_V(&rsq0, incv[i], cgsub);
      rsqsave[i] = res_arg[i][0]->rsqmin;
      res_arg[i][0]->rsqmin *= rsq0/rsq1;
    }

    QOPPC(invert_cgms_V)(QOPPC(asqtad_invert_d2), inv_arg, ra, shifts, 2,
			 outcv, src, cgp, cgsub);

    info->final_flop += (nflop+nflopm)*ra[imin]->final_iter*QDP_sites_on_node;

    for(i=0; i<2; i++) {
      res_arg[i][0]->rsqmin = rsqsave[i];
      QDP_V_peq_V(outcv[i], x0, cgsub);
    }

    QDP_destroy_V(x0);
    QDP_destroy_V(src);
  } else { // regular multimass
    for(i=0; i<nsrc; i++) {
      //QLA_Real shifts[nmass[i]];
      //QDP_ColorVector *cv[nmass[i]];
      // work around bug in XLC
      QLA_Real *shifts;
      QDP_ColorVector **cv;
      int jmin=0;
      shifts = (QLA_Real *) malloc(nmass[i]*sizeof(QLA_Real));
      cv = (QDP_ColorVector **) malloc(nmass[i]*sizeof(QDP_ColorVector *));
      for(j=0; j<nmass[i]; j++) {
	shifts[j] = masses[i][j]*masses[i][j];
	if(shifts[j]<shifts[jmin]) jmin = j;
	cv[j] = out_pt[i][j]->cv;
      }
      gl_mass = masses[i][jmin];
      for(j=0; j<nmass[i]; j++) shifts[j] -= gl_mass*gl_mass;

      QOPPC(invert_cgms_V)(QOPPC(asqtad_invert_d2), inv_arg, res_arg[i],
			   shifts, nmass[i], cv, in_pt[i]->cv, cgp, cgsub);

      info->final_flop += (nflop+nflopm*(nmass[i]-1))
	* res_arg[i][0]->final_iter * QDP_sites_on_node;
      free(shifts);
      free(cv);
    }
  }
#endif  // real multimass

  dtime += QDP_time();

  for(i=0; i<nsrc; i++) {
    for(j=0; j<nmass[i]; j++) {
      //QLA_Real minv = 1.0/masses[i][j];
      //QDP_V_eq_r_times_V(out_pt[i][j]->cv, &minv, out_pt[i][j]->cv, cgsub);
      QLA_Real m = masses[i][j];
      QDP_V_eq_r_times_V(out_pt[i][j]->cv, &m, out_pt[i][j]->cv, cgsub);
    }
  }

  info->final_sec = dtime;
  info->status = QOP_SUCCESS;
}
