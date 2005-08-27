#include <sys/time.h>
#include <qop.h>

#if QOP_Precision == 2
#define PREC(x) QOP_D_##x
#define REAL double
#else
#define PREC(x) QOP_F_##x
#define REAL float
#endif

static int congrad_setup = 0;
static QDP_ColorVector *ttt, *tttt, *r, *p;
static int have_rawlinks = 0;
static QDP_ColorMatrix *raw_fatlinks[4], *raw_longlinks[4];
static QDP_ColorMatrix *fwdlinks[8];
static QDP_ColorVector *temp1[16], *temp2[16];
static QDP_ColorVector *fb_tempvec[8];
static int temps_inited = 0;
static int fb_tempvec_inited = 0;
static QDP_Shift asqtad_shifts[8], neighbor3[4];
static QDP_ShiftDir shiftfwd[8], shiftbck[8];

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

  for(i=0; i<4; i++) {
    fwdlinks[2*i] = fatlnk[i];
    fwdlinks[2*i+1] = longlnk[i];
  }

  if(!congrad_setup) {
    int disp[4]={0,0,0,0};
    congrad_setup = 1;
    ttt = QDP_create_V();
    tttt = QDP_create_V();
    r = QDP_create_V();
    p = QDP_create_V();
    for(i=0; i<4; i++) {
      asqtad_shifts[2*i] = QDP_neighbor[i];
      disp[i] = 3;
      neighbor3[i] = QDP_create_shift(disp);
      disp[i] = 0;
      asqtad_shifts[2*i+1] = neighbor3[i];
    }
    for(i=0; i<8; i++) {
      shiftfwd[i] = QDP_forward;
      shiftbck[i] = QDP_backward;
    }
  }

  if(!temps_inited) {
    temps_inited = 1;
    for(i=0; i<16; i++) {
      temp1[i] = QDP_create_V();
      temp2[i] = QDP_create_V();
    }
  }

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
  dtimec = -dclock();

  do {
    QDP_V_eq_V(p, out, subset);
    PREC(asqtad_mdslash2)(ttt, p, subset, othersubset, m2x4);
    iteration++;

    QDP_V_eq_V_plus_V(r, in, ttt, subset);
    QDP_r_eq_norm2_V(&insq, in, subset);
    QDP_r_eq_norm2_V(&rsq, r, subset);
    oldrsq = rsq;
    QDP_V_eq_zero(p, subset);
    rsqstop = inv_arg->rsqmin * insq;

    while( (rsq>rsqstop) && (iteration%inv_arg->restart!=0) ) {
      //printf("iteration = %-5i resid = %g\n", iteration, rsq);
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

  if( rsq <= rsqstop ) {
    dtimec += dclock();
    if(QDP_this_node==0) {
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
  return QOP_SUCCESS;
}


/* internal functions */

static void
PREC(asqtad_dslash)(QDP_ColorVector *dest, QDP_Subset dsubset,
		    QDP_ColorVector *src, QDP_Subset ssubset,
		    QDP_ColorVector *temp[])
{
  QDP_ColorVector *srcvec[8], *destvec[8];
  int i;

  for(i=0; i<8; i++) {
    srcvec[i] = src;
    destvec[i] = dest;
  }

  if(!fb_tempvec_inited) {
    fb_tempvec_inited = 1;
    for(i=0; i<8; i++) {
      fb_tempvec[i] = QDP_create_V();
    }
  }

  /* Start gathers from positive directions */
  QDP_V_veq_sV(temp, srcvec, asqtad_shifts, shiftfwd, dsubset, 8);

  /* Multiply by adjoint matrix at other sites */
  QDP_V_veq_Ma_times_V(fb_tempvec, fwdlinks, srcvec, ssubset, 8);

  /* Start gathers from negative directions */
  QDP_V_veq_sV(temp+8, fb_tempvec, asqtad_shifts, shiftbck, dsubset, 8);

  /* Wait gathers from positive directions, multiply by matrix and
     accumulate */
  QDP_V_eq_zero(dest, dsubset);
  QDP_V_vpeq_M_times_V(destvec, fwdlinks, temp, dsubset, 8);
  for(i=0; i<8; ++i) QDP_discard_V(temp[i]);

  /* Wait gathers from negative directions, accumulate (negative) */
  QDP_V_vmeq_V(destvec, temp+8, dsubset, 8);
  for(i=0; i<8; ++i) QDP_discard_V(temp[i+8]);
}

static void
PREC(asqtad_mdslash2)(QDP_ColorVector *out, QDP_ColorVector *in,
		      QDP_Subset subset, QDP_Subset othersubset, QLA_Real m2x4)
{
  PREC(asqtad_dslash)(tttt, othersubset, in, subset, temp1);
  PREC(asqtad_dslash)(out, subset, tttt, othersubset, temp2);
  QDP_V_meq_r_times_V(out, &m2x4, in, subset);
}
