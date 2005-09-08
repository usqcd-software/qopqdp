#include <sys/time.h>
#include <qop.h>

#if QOP_Precision == 2
#define PREC(x) QOP_D_##x
#define REAL double
#else
#define PREC(x) QOP_F_##x
#define REAL float
#endif

#define printf0 if(QDP_this_node==0) printf

extern int QOP_asqtad_inited;
extern int QOP_style;
extern int QOP_nsvec;
extern int QOP_nvec;
extern QDP_Shift QOP_asqtad_shifts[8], QOP_neighbor3[4];
extern QDP_Shift QOP_asqtad_shifts_dbl[16];
extern QDP_ShiftDir QOP_shiftfwd[8], QOP_shiftbck[8];
extern QDP_ShiftDir QOP_asqtad_shiftdirs_dbl[16];

static int old_style=0;
static int old_nsvec=0;
static int old_nvec=0;

static int congrad_setup = 0;
static int dblstored = 0;
static QDP_ColorMatrix *fwdlinks[8];
static QDP_ColorMatrix *bcklinks[8], *dbllinks[16];
static QDP_ColorVector *ttt, *tttt, *r, *p;
static QDP_ColorVector *temp1[24], *temp2[24];

static int have_rawlinks = 0;
static QDP_ColorMatrix *raw_fatlinks[4], *raw_longlinks[4];

static void PREC(asqtad_mdslash2)(QDP_ColorVector *out, QDP_ColorVector *in,
				  QDP_Subset subset, QDP_Subset othersubset,
				  QLA_Real m2x4);

double
dclock(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6*tv.tv_usec;
}

static void
load_link(QDP_ColorMatrix *cm, REAL *lnk)
{
  QDP_insert_M(cm, (QLA_ColorMatrix *)lnk, QDP_all);
}

static void
set_V_from_real(QDP_ColorVector *q, REAL *r)
{
  QDP_insert_V(q, (QLA_ColorVector *)r, QDP_all);
}

static void
set_real_from_V(REAL *r, QDP_ColorVector *q)
{
  QDP_extract_V((QLA_ColorVector *)r, q, QDP_all);
}

static void
free_temps(void)
{
  if(congrad_setup) {
    int i, n;

    QDP_destroy_V(ttt);
    QDP_destroy_V(tttt);
    QDP_destroy_V(r);
    QDP_destroy_V(p);

    if(old_style==0) n = 24;
    else n = 16;
    for(i=0; i<n; i++) {
      QDP_destroy_V(temp1[i]);
      QDP_destroy_V(temp2[i]);
    }
  }
  congrad_setup = 0;
}

static void
double_store(void)
{
  if( (QOP_style==1) && (!dblstored) ) {
    int i;
    QDP_ColorMatrix *m = QDP_create_M();
    for(i=0; i<8; i++) {
      QDP_M_eq_sM(m, fwdlinks[i], QOP_asqtad_shifts[i], QDP_backward, QDP_all);
      QDP_M_eqm_Ma(bcklinks[i], m, QDP_all);
    }
    QDP_destroy_M(m);
    dblstored = 1;
  }
}

static void
reset_temps(void)
{
  int i, n;

  if(QOP_style!=old_style) {
    if(QOP_style==0) {
      if(congrad_setup) {
	for(i=0; i<8; i++) {
	  QDP_destroy_M(bcklinks[i]);
	}
      }
    } else {
      for(i=0; i<8; i++) {
	bcklinks[i] = QDP_create_M();
      }
      for(i=0; i<4; i++) {
	dbllinks[4*i] = fwdlinks[2*i];
	dbllinks[4*i+1] = fwdlinks[2*i+1];
	dbllinks[4*i+2] = bcklinks[2*i];
	dbllinks[4*i+3] = bcklinks[2*i+1];
      }
    }
    dblstored = 0;
  }
  double_store();

  free_temps();

  ttt = QDP_create_V();
  tttt = QDP_create_V();
  r = QDP_create_V();
  p = QDP_create_V();

  if(QOP_style==0) n = 24;
  else n = 16;
  for(i=0; i<n; i++) {
    temp1[i] = QDP_create_V();
    temp2[i] = QDP_create_V();
  }
  congrad_setup = 1;
}

QOP_status_t
PREC(asqtad_invert_load_links_raw)(REAL *fatlinks[], REAL *longlinks[])
{
  int i;

  if(!have_rawlinks) {
    have_rawlinks = 1;
    for(i=0; i<4; i++) {
      raw_fatlinks[i] = QDP_create_M();
      raw_longlinks[i] = QDP_create_M();
    }
  }
  for(i=0; i<4; i++) {
    load_link(raw_fatlinks[i], fatlinks[i]);
    load_link(raw_longlinks[i], longlinks[i]);
  }
  PREC(asqtad_invert_load_links_qdp)(raw_fatlinks, raw_longlinks);

  return QOP_SUCCESS;
}

QOP_status_t
PREC(asqtad_invert_load_links_qdp)(QDP_ColorMatrix *fatlnk[],
				   QDP_ColorMatrix *longlnk[])
{
  int i;

  dblstored = 0;
  for(i=0; i<4; i++) {
    fwdlinks[2*i] = fatlnk[i];
    fwdlinks[2*i+1] = longlnk[i];
  }

  if(!congrad_setup) {
    reset_temps();
  }

  if( (QOP_style != old_style) ||
      (QOP_nsvec != old_nsvec) ||
      (QOP_nvec != old_nvec) ) {
    reset_temps();
    old_style = QOP_style;
    old_nsvec = QOP_nsvec;
    old_nvec = QOP_nvec;
  }

  double_store();

  return QOP_SUCCESS;
}

int
PREC(asqtad_inv_raw)(QOP_invert_arg *inv_arg, REAL *out_pt, REAL *in_pt)
{
  QDP_ColorVector *in, *out;
  int nit;

  in = QDP_create_V();
  out = QDP_create_V();
  set_V_from_real(in, in_pt);
  set_V_from_real(out, out_pt);
  nit = PREC(asqtad_inv_qdp)(inv_arg, out, in);
  set_real_from_V(out_pt, out);

  return nit;
}

int
PREC(asqtad_inv_qdp)(QOP_invert_arg *inv_arg,
		     QDP_ColorVector *out, QDP_ColorVector *in)
{
  QLA_Real a, b;
  QLA_Real rsq, oldrsq, pkp;
  QLA_Real m2x4;
  QLA_Real insq;
  QLA_Real rsqstop;
  QDP_Subset subset, othersubset;
  int iteration=0;

  double dtimec;
  double nflop;
  nflop = 1187;

  if(inv_arg->evenodd==QOP_EVEN) {
    subset = QDP_even;
    othersubset = QDP_odd;
  } else if(inv_arg->evenodd==QOP_ODD) {
    subset = QDP_odd;
    othersubset = QDP_even;
  } else {
    subset = QDP_all;
    othersubset = QDP_all;
    nflop *= 2;
  }

  m2x4 = 4*inv_arg->mass*inv_arg->mass;

  if( (QOP_style != old_style) ||
      (QOP_nsvec != old_nsvec) ||
      (QOP_nvec != old_nvec) ) {
    reset_temps();
    old_style = QOP_style;
    old_nsvec = QOP_nsvec;
    old_nvec = QOP_nvec;
  }

  dtimec = -dclock();

  do {
    QDP_V_eq_V(p, out, subset);
    //printf0("first dirac op\n"); fflush(stdout);
    PREC(asqtad_mdslash2)(ttt, p, subset, othersubset, m2x4);
    //printf0("done first dirac op\n"); fflush(stdout);
    iteration++;

    QDP_V_eq_V_plus_V(r, in, ttt, subset);
    QDP_r_eq_norm2_V(&insq, in, subset);
    QDP_r_eq_norm2_V(&rsq, r, subset);
    oldrsq = rsq;
    QDP_V_eq_zero(p, subset);
    rsqstop = inv_arg->rsqmin * insq;

    while( (rsq>rsqstop) && (iteration%inv_arg->restart!=0) ) {
      //printf0("iteration = %-5i resid = %g\n", iteration, rsq);
      b = rsq / oldrsq;
      QDP_V_eq_r_times_V_plus_V(p, &b, p, r, subset);

      oldrsq = rsq;

      PREC(asqtad_mdslash2)(ttt, p, subset, othersubset, m2x4);
      iteration++;

      QDP_r_eq_re_V_dot_V(&pkp, p, ttt, subset);

      a = - rsq / pkp;

      QDP_V_peq_r_times_V(out, &a, p, subset);
      QDP_V_peq_r_times_V(r, &a, ttt, subset);
      QDP_r_eq_norm2_V(&rsq, r, subset);
    }

  } while( (rsq>rsqstop) && (iteration<inv_arg->max_iter) );

  dtimec += dclock();
  inv_arg->final_rsq = rsq;
  inv_arg->final_iter = iteration;
  inv_arg->final_sec = dtimec;
  inv_arg->final_flop = nflop*iteration*QDP_sites_on_node;
  if( rsq <= rsqstop ) {
    if(QDP_this_node==-1) {
      printf("CONGRAD5: time = %g iters = %i mflops = %g\n",
             dtimec, iteration,
             (double) ( nflop*iteration*QDP_sites_on_node / (1.0e6*dtimec) ));
      fflush(stdout);
    }
  } else {
    if(QDP_this_node==0) {
      printf("CG not converged after %d iterations, res. = %e wanted %e\n",
             iteration, rsq, rsqstop);
      fflush(stdout);
    }
  }

  return(iteration);
}

QOP_status_t
PREC(asqtad_invert_unload_links)(void)
{
  dblstored = 0;
  return QOP_SUCCESS;
}


/* internal functions */

static void
PREC(asqtad_dslash0)(QDP_ColorVector *dest, QDP_Subset dsubset,
		     QDP_ColorVector *src, QDP_Subset ssubset,
		     QDP_ColorVector *temp[])
{
  QDP_ColorVector *vsrc[8], *vdest[8];
  int i;

  for(i=0; i<8; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  /* Start gathers from positive directions */
  for(i=0; i<8; i+=QOP_nsvec) {
    QDP_V_veq_sV(temp+i, vsrc+i, QOP_asqtad_shifts+i, QOP_shiftfwd+i,
		 dsubset, QOP_nsvec);
  }

  /* Multiply by adjoint matrix at other sites */
  /* Start gathers from negative directions */
  for(i=0; i<8; i+=QOP_nsvec) {
    QDP_V_veq_Ma_times_V(temp+16+i, fwdlinks+i, vsrc+i, ssubset, QOP_nsvec);
    QDP_V_veq_sV(temp+8+i, temp+16+i, QOP_asqtad_shifts+i, QOP_shiftbck+i,
		 dsubset, QOP_nsvec);
  }

  /* Wait gathers from positive directions, multiply by matrix and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<8; i+=QOP_nvec) {
    QDP_V_vpeq_M_times_V(vdest+i, fwdlinks+i, temp+i, dsubset, QOP_nvec);
  }
  for(i=0; i<8; i++) QDP_discard_V(temp[i]);

  /* Wait gathers from negative directions, accumulate (negative) */
  for(i=0; i<8; i+=QOP_nvec) {
    QDP_V_vmeq_V(vdest+i, temp+8+i, dsubset, QOP_nvec);
  }
  for(i=0; i<8; ++i) QDP_discard_V(temp[i+8]);
}

static void
PREC(asqtad_dslash1)(QDP_ColorVector *dest, QDP_Subset dsubset,
		     QDP_ColorVector *src, QDP_Subset ssubset,
		     QDP_ColorVector *temp[])
{
  QDP_ColorVector *vsrc[16], *vdest[16];
  int i;

  for(i=0; i<16; i++) {
    vsrc[i] = src;
    vdest[i] = dest;
  }

  /* Start gathers from all directions */
  for(i=0; i<16; i+=QOP_nsvec) {
    QDP_V_veq_sV(temp+i, vsrc+i, QOP_asqtad_shifts_dbl+i,
		 QOP_asqtad_shiftdirs_dbl+i, dsubset, QOP_nsvec);
  }

  /* Wait gathers from all directions, multiply by matrix and accumulate */
  QDP_V_eq_zero(dest, dsubset);
  for(i=0; i<16; i+=QOP_nvec) {
    QDP_V_vpeq_M_times_V(vdest+i, dbllinks+i, temp+i, dsubset, QOP_nvec);
  }
  for(i=0; i<16; ++i) QDP_discard_V(temp[i]);
}

static void
PREC(asqtad_mdslash2)(QDP_ColorVector *out, QDP_ColorVector *in,
		      QDP_Subset subset, QDP_Subset othersubset, QLA_Real m2x4)
{
  if(QOP_style==0) {
    PREC(asqtad_dslash0)(tttt, othersubset, in, subset, temp1);
    PREC(asqtad_dslash0)(out, subset, tttt, othersubset, temp2);
  } else {
    PREC(asqtad_dslash1)(tttt, othersubset, in, subset, temp1);
    PREC(asqtad_dslash1)(out, subset, tttt, othersubset, temp2);
  }
  QDP_V_meq_r_times_V(out, &m2x4, in, subset);
}
