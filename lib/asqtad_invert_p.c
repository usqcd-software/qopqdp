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

If we have the approximate preconditioned solution, v,
the preconditioned residual, s, is

s = z - A v

and the full residual is

r = (x) - D ( v              ) = ( s/m )
    (y)     ( y/m - D_ba v/m )   ( 0 )

the stopping criterion is

|r|^2/(|x|^2+|y|^2) < rsqmin  (outer)
|s|^2/|z|^2         < rsqmin' (inner)

equating these gives

rsqmin' = rsqmin (|x|^2+|y|^2) m^2 / |z|^2

***************************************************************************/
#define DO_TRACE
#include <qop_internal.h>
#include <math.h>

extern QOP_asqtad_t QOP_asqtad;
extern QDP_ColorVector *QOP_asqtad_dslash_get_tmp
     (QOP_FermionLinksAsqtad *fla, QOP_evenodd_t eo, int n);

/* inverter stuff */

static QOP_FermionLinksAsqtad *gl_fla;
static REAL gl_mass;
static QOP_evenodd_t gl_eo;
static QDP_ColorVector *gl_tmp, *gl_tmp2;

#if 0
static void
dumpvec(QDP_ColorVector *in, QDP_Subset sub)
{
  QLA_Real nrm2;
  QDP_r_eq_norm2_V(&nrm2, in, sub);
  QOP_printf0("norm2 = %g\n", nrm2);
  int i;
  QDP_loop_sites(i, sub, {
      int c[4];
      QDP_get_coords(c, 0, i);
      printf(" %i %i %i %i", c[0], c[1], c[2], c[3]);
      for(int j=0; j<QLA_Nc; j++) {
	QLA_ColorVector *v = QDP_site_ptr_readonly_V(in, i);
	printf(" %10g", QLA_real(QLA_elem_V(*v, j)));
	printf(" %10g", QLA_imag(QLA_elem_V(*v, j)));
      }
      printf("\n");
    });
  exit(0);
}
#endif

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

static void
QOP_asqtad_invert_d2(QDP_ColorVector *out, QDP_ColorVector *in,
		     QDP_Subset subset)
{
  QOP_asqtad_ddagd(NULL, gl_fla, gl_mass, out, in, gl_eo);
}

static QOP_Real
QOP_asqtad_invert_d2_norm2(QDP_ColorVector *out, QDP_ColorVector *in,
			   QDP_Subset subset)
{
  return QOP_asqtad_ddagd_norm2(NULL, gl_fla, gl_mass, out, in, gl_eo);
}

#if QOP_Precision == 'D'
static QOP_F_FermionLinksAsqtad *gl_flaf;
#if 0
static void
QOP_F_asqtad_invert_d2(QDP_F_ColorVector *out, QDP_F_ColorVector *in,
		       QDP_Subset subset)
{
  QOP_F_asqtad_ddagd(NULL, gl_flaf, gl_mass, out, in, gl_eo);
}
#endif
static QOP_F_Real
QOP_F_asqtad_invert_d2_norm2(QDP_F_ColorVector *out, QDP_F_ColorVector *in,
			     QDP_Subset subset)
{
  return QOP_F_asqtad_ddagd_norm2(NULL, gl_flaf, gl_mass, out, in, gl_eo);
}
#endif

void
QOP_asqtad_invert(QOP_info_t *info,
		  QOP_FermionLinksAsqtad *fla,
		  QOP_invert_arg_t *inv_arg,
		  QOP_resid_arg_t *res_arg,
		  REAL mass,
		  QOP_ColorVector *out,
		  QOP_ColorVector *in)
{
  QOP_asqtad_invert_qdp(info, fla, inv_arg, res_arg, mass, out->cv, in->cv);
}

#if 0
void
asqtad_invert_threadfunc_qdp(void *args0)
{
  void **args = args0;
  QOP_info_t *info = args[0];
  QOP_FermionLinksAsqtad *fla = args[1];
  QOP_invert_arg_t *inv_arg = args[2];
  QOP_resid_arg_t *res_arg = args[3];
  REAL mass = *(REAL*)args[4];
  QDP_ColorVector *out = args[5];
  QDP_ColorVector *in = args[6];
  QOP_asqtad_invert_qdp(info, fla, inv_arg, res_arg, mass, out, in);
  //QDP_thread_barrier();
}

void
QOP_asqtad_invert_threaded_qdp(QOP_info_t *info,
			       QOP_FermionLinksAsqtad *fla,
			       QOP_invert_arg_t *inv_arg,
			       QOP_resid_arg_t *res_arg,
			       REAL mass,
			       QDP_ColorVector *out,
			       QDP_ColorVector *in,
			       int nthreads)
{
  void *args[7];
  args[0] = info;
  args[1] = fla;
  args[2] = inv_arg;
  args[3] = res_arg;
  args[4] = &mass;
  args[5] = out;
  args[6] = in;
  //QDP_create_threads(nthreads, 1, asqtad_invert_threadfunc_qdp, args);
}

void
QOP_asqtad_invert_threaded(QOP_info_t *info,
			   QOP_FermionLinksAsqtad *fla,
			   QOP_invert_arg_t *inv_arg,
			   QOP_resid_arg_t *res_arg,
			   REAL mass,
			   QOP_ColorVector *out,
			   QOP_ColorVector *in,
			   int nthreads)
{
  QOP_asqtad_invert_threaded_qdp(info, fla, inv_arg, res_arg, mass, out->cv, in->cv, nthreads);
}
#endif

static void
QOP_asqtad_solve_multiA(QOP_info_t *info,
			QOP_FermionLinksAsqtad *fla,
			QOP_invert_arg_t *inv_arg,
			QOP_resid_arg_t *res_arg[],
			REAL masses[],
			QDP_ColorVector *out[],
			QDP_ColorVector *in[],
			int nsolve)
{
  for(int i=0; i<nsolve; i++) {
    QOP_asqtad_invert_qdp(info, fla, inv_arg, res_arg[i], masses[i],
			  out[i], in[i]);
  }
#if 0
  for(int i=0; i<nsolve; i++) {
    if(ineo==QOP_ODD) reconstruct(fla, mass[i], y[i], x[i], x[i], QOP_EVEN);
    else reconstruct(fla, mass[i], y[i], x[i], in[i], QOP_EVEN);
    if(ineo==QOP_EVEN) reconstruct(fla, mass[i], y[i], x[i], x[i], QOP_ODD);
    else reconstruct(fla, mass[i], y[i], x[i], in[i], QOP_ODD);
  }
  {
    // update x from y e/o

    for(int i=0; i<nsolve; i++) {
      QOP_asqtad_dslash_qdp(NULL, fla, mass[i], r[i], x[i], QOP_EVENODD, QOP_EVENODD);
      QDP_V_meq_V(r[i], in[i], insub);
    }
    if(ineo!=QOP_ODD) QDP_r_veq_norm2_V(innrm2e, in, QOP_EVEN, nsolve);
    else for(int i=0; i<nsolve; i++) innrm2e[i] = 0;
    if(ineo!=QOP_EVEN) QDP_r_veq_norm2_V(innrm2o, in, QOP_ODD, nsolve);
    else for(int i=0; i<nsolve; i++) innrm2o[i] = 0;
    for(int i=0; i<nsolve; i++) {
      innrm2et += innrm2
	innrm2[i] = innrm2e[i] + innrm2o[i];
    }
    // run preconditioned solve on larger input subset
    if(insqe>=insqo) cgeo = QOP_EVEN;
    else cgeo = QOP_ODD;
    cgsub = qdpsub(cgeo);
    gl_eo = cgeo;
    for(int i=0; i<nsolve; i++) {
      QOP_asqtad_dslash_qdp(NULL, fla, mass[i], r[i], x[i], QOP_EVENODD, QOP_EVENODD);
      QDP_V_meq_V(r[i], in[i], insub);
    }
  }
#endif
}

// r[i][opsub] = A out[i][insub]
static void
refine(QDP_ColorVector *out[], QDP_ColorVector *in[], QDP_ColorVector *r[],
       QLA_Real r2[], QLA_Real m2[], QLA_Complex scale[], int no,
       QOP_evenodd_t ineo, QDP_ColorVector *x[], QDP_ColorVector *Ax[], int nx,
       int oi[], int xi[], int ni, QDP_ColorVector *tv)
{
  QDP_Subset insub = qdpsub(ineo);
  QOP_evenodd_t opeo = oppsub(ineo);
  QDP_Subset opsub = qdpsub(opeo);
  QDP_ColorVector **Aout = r;
  // Ax[i][opsub] = A x[i][insub]
  for(int i=0; i<nx; i++) {
    QOP_asqtad_dslash_qdp(NULL, gl_fla, 0, Ax[i], x[i], opeo, ineo);
  }
  QLA_Real x2[nx+1], Ax2[nx+1];
  QDP_r_veq_norm2_V(x2, x, insub, nx);
  QDP_r_veq_norm2_V(Ax2, Ax, opsub, nx);
  for(int i=0; i<ni; i++) {
    int ko = oi[i];
    int kx = 0;
    if(nx>0) kx = xi[i];
    QLA_Complex z, Az, s, sx;
    QLA_c_eq_r(sx, 0);
    // A-orthogonalize out[ko] against x[kx]
    double d = 0;
    if(nx>0) {
      QDP_c_eq_V_dot_V(&z, x[kx], out[ko], insub);
      QDP_c_eq_V_dot_V(&Az, Ax[kx], Aout[ko], opsub);
      d = Ax2[kx] + m2[ko]*x2[kx];
      //QOP_printf0("d: %g\n", d);
    }
    VERB(MED, "REFINE: d: %g\n", d);
    if(d!=0) {
      double di = 1/d;
      double dim2 = di*m2[ko];
      QLA_c_eq_r_times_c(s, dim2, z);
      QLA_c_peq_r_times_c(s, di, Az);
      QDP_V_meq_c_times_V(out[ko], &s, x[kx], insub);
      QDP_V_meq_c_times_V(Aout[ko], &s, Ax[kx], opsub);
      // get scale factor for x[kx] to minimize the residual
      QDP_c_eq_V_dot_V(&z, x[kx], in[ko], insub);
      QLA_c_eq_r_times_c(sx, di, z);
      //QOP_printf0("sx: %g  %g\n", QLA_real(sx), QLA_imag(sx));
    }
    // get scale factor for out[ko] to minimize the residual
    QLA_Real o2, Ao2;
    QDP_c_eq_V_dot_V(&z, out[ko], in[ko], insub);
    QDP_r_eq_norm2_V(&o2, out[ko], insub);
    QDP_r_eq_norm2_V(&Ao2, Aout[ko], opsub);
    d = Ao2 + m2[kx]*o2;
    VERB(MED, "REFINE: d: %g  o2: %g  Ao2: %g\n", d, o2, Ao2);
    if(d!=0) {
      double di = 1/d;
      QLA_c_eq_r_times_c(s, di, z);
      QLA_c_eq_c(scale[ko], s);
      // out[ko] <- out[ko] + (sx/s)*x[kx]
      if(nx>0) {
	QLA_c_eq_c_div_c(z, sx, s);
	QDP_V_peq_c_times_V(out[ko], &z, x[kx], insub);
	QDP_V_peq_c_times_V(Aout[ko], &z, Ax[kx], opsub);
      }
    } else {
      QLA_c_eq_c(scale[ko], sx);
      if(nx>0) {
	QDP_V_eq_V(out[ko], x[kx], insub);
	QDP_V_eq_V(Aout[ko], Ax[kx], opsub);
      }
    }
  }
  // get residuals (r <- in - D2 scale out)
  for(int i=0; i<no; i++) {
    QOP_asqtad_dslash_qdp(NULL, gl_fla, 0, tv, Aout[i], ineo, opeo);
    QDP_V_meq_r_times_V(tv, &m2[i], out[i], insub);
    QDP_V_eq_c_times_V_plus_V(r[i], &scale[i], tv, in[i], insub);
  }
  QDP_r_veq_norm2_V(r2, r, insub, no);
}

void
QOP_asqtad_solve_multi_qdp(QOP_info_t *info,
			   QOP_FermionLinksAsqtad *fla,
			   QOP_invert_arg_t *inv_arg,
			   QOP_resid_arg_t *res_arg[],
			   REAL masses[],
			   QDP_ColorVector *out[],
			   QDP_ColorVector *in[],
			   int nsolve)
{
#define NC QDP_get_nc(in[0])
  QOP_evenodd_t ineo = inv_arg->evenodd;
  if(ineo==QOP_EVENODD) {
    QOP_asqtad_solve_multiA(info, fla, inv_arg,	res_arg, masses,
			    out, in, nsolve);
    return;
  }
  QDP_ColorVector *tv, *r[nsolve];
  QLA_Complex scale[nsolve];
  double dtime = 0;
  QLA_Real in2[nsolve];
  QLA_Real r2[nsolve];
  QLA_Real r2stop[nsolve];
  QLA_Real mass2[nsolve];
  QLA_Real m2s[nsolve];
  int oi[nsolve], xi[nsolve];
  int nm = 0;
  int iter = 0, nrestart = -1;
  int max_iter = inv_arg->max_iter;
  int max_restarts = inv_arg->max_restarts;
  QDP_Subset insub = qdpsub(ineo);
  QOP_evenodd_t opeo = oppsub(ineo);
#if QOP_Precision == 'D'
  double mixedrsq = inv_arg->mixed_rsq;
#endif

  gl_fla = fla;
  gl_mass = 0;
  gl_eo = ineo;
  //gl_tmp2 = QOP_asqtad_dslash_get_tmp(fla, ineo, 1);
  QDP_ColorVector *cgp = QOP_asqtad_dslash_get_tmp(fla, opeo, 2);

  QDP_r_veq_norm2_V(in2, in, insub, nsolve);

  tv = QDP_create_V();
  for(int i=0; i<nsolve; i++) {
    r2stop[i] = res_arg[i]->rsqmin * in2[i];
    mass2[i] = masses[i]*masses[i];
    oi[i] = i;
    int j = nm;
    while(j>0 && mass2[i]<m2s[j-1]) j--;
    if(j==nm || m2s[j]!=mass2[i]) {
      for(int k=nm; k>j; k--) m2s[k] = m2s[k-1];
      m2s[j] = mass2[i];
      nm++;
    }
    r[i] = QDP_create_V();
    QOP_asqtad_dslash_qdp(NULL, fla, 0, r[i], out[i], opeo, ineo);
  }
  for(int i=0; i<nsolve; i++) {
    int j = 0;
    while(j<nm-1 && m2s[j]!=mass2[i]) j++;
    xi[i] = j;
  }
  QOP_resid_arg_t *xresarg[nm];
  QDP_ColorVector *x[nm];
  double m20 = m2s[0];
  gl_mass = sqrt(m20);
  for(int i=0; i<nm; i++) {
    m2s[i] -= m20;
    x[i] = QDP_create_V();
  }

#if 1
  refine(out, in, r, r2, mass2, scale, nsolve, ineo, NULL, NULL, 0,
	 oi, NULL, nsolve, tv);
#else
  for(int i=0; i<nsolve; i++) {
    int oi = 1;
    refine(out+i, in+i, r+i, r2+i, mass2+i, scale+i, 1, ineo, in+i, &tv, 1,
	   &oi, &oi, 1, tv);
  }
#endif

  while(1) {
    // find least converged (or smaller mass if tie)
    int imax = 0;
    double r2max=1, r2stopmax=1, in2max=1;
    VERB(MED, "SOLVE: its: %i\n", iter);
    for(int i=0; i<nsolve; i++) {
      // find the stopping criterion that is closer to convergence for this 'i'
      double r2i = r2[i];
      double r2stopi = r2stop[i];
      double in2i = in2[i];
      if(res_arg[i]->final_rel*r2stopi<r2i*res_arg[i]->relmin) {
	r2i = res_arg[i]->final_rel;
	r2stopi = res_arg[i]->relmin;
	in2i = 1;
      }
      // update if less converged
      double a = r2i*r2stopmax;
      double b = r2max*r2stopi;
      if(a>b || (a==b && mass2[i]<mass2[imax])) {
	imax = i;
	r2max = r2i;
	r2stopmax = r2stopi;
	in2max = in2i;
      }
      VERB(MED, "SOLVE[%i]: rsq = %g rel = %g\n", i,
	   r2[i]/in2[i], res_arg[i]->final_rel);
    }
    VERB(LOW, "SOLVE: its: %i worst: %i mass: %g rsq: %g/%g rel: %g/%g\n",
	 iter, imax, masses[imax], r2[imax]/in2[imax], res_arg[imax]->rsqmin,
	 res_arg[imax]->final_rel, res_arg[imax]->relmin);
    // if least converged has converged, we're done
    if( r2max <= r2stopmax ||
	nrestart >= max_restarts ||
	iter >= max_iter ) break;

    // solve
    for(int i=0; i<nm; i++) {
      xresarg[i] = res_arg[imax];
      QDP_V_eq_zero(x[i], QDP_all);
    }
    double rsqminold = res_arg[imax]->rsqmin;
    //res_arg[imax]->rsqmin = r2stop[imax]/r2[imax];
    res_arg[imax]->rsqmin = r2stopmax/r2max;
    if(res_arg[imax]->rsqmin>=1) {
      res_arg[imax]->rsqmin = in2max*r2stopmax/r2max;
    }
    inv_arg->max_iter = max_iter - iter;
    //QOP_verbose(3);
#if QOP_Precision == 'D'
    int mixed = (r2max > mixedrsq*in2max);
    if(mixed) {
      QOP_F_FermionLinksAsqtad *flaf = QOP_FD_asqtad_create_L_from_L(fla);
      gl_flaf = flaf;
      float m2sf[nm];
      QDP_F_ColorVector *xf[nm];
      QDP_F_ColorVector *rf = QDP_F_create_V();
      QDP_FD_V_eq_V(rf, r[imax], insub);
      for(int i=0; i<nm; i++) {
	m2sf[i] = m2s[i];
	xf[i] = QDP_F_create_V();
	//QDP_FD_V_eq_V(xf[i], x[i], insub);
      }
      QDP_F_ColorVector *cgpf = QOPFC(asqtad_dslash_get_tmp)(flaf, opeo, 2);
      if(mixedrsq>rsqminold) {
	res_arg[imax]->rsqmin = mixedrsq*in2max/r2max;
      }
      dtime -= QOP_time();
      QOPFC(invert_cgms_V)(QOP_F_asqtad_invert_d2_norm2, inv_arg, xresarg,
			   m2sf, nm, xf, rf, cgpf, insub);
      dtime += QOP_time();
      for(int i=0; i<nm; i++) {
	QDP_DF_V_eq_V(x[i], xf[i], insub);
	QDP_F_destroy_V(xf[i]);
      }
      QDP_F_destroy_V(rf);
      QOP_F_asqtad_destroy_L(flaf);
    } else
#endif
      {
	dtime -= QOP_time();
	QOP_invert_cgms_V(QOP_asqtad_invert_d2_norm2, inv_arg, xresarg,
			  m2s, nm, x, r[imax], cgp, insub);
	dtime += QOP_time();
      }
    iter += res_arg[imax]->final_iter;
    nrestart++;
    res_arg[imax]->rsqmin = rsqminold;

#if 0
    for(int i=0; i<nm; i++) {
      QOP_asqtad_dslash_qdp(NULL, fla, 0, x[i], x[i], opeo, ineo);
      QOP_asqtad_dslash_qdp(NULL, fla, 0, tv, x[i], ineo, opeo);
      QDP_V_meq_r_times_V(tv, &mass2[i], x[i], insub);
      QDP_V_eq_V_plus_V(r[i], tv, r[imax], insub);
    }
    QDP_r_veq_norm2_V(r2, r, insub, nm);
    for(int i=0; i<nm; i++) {
      QOP_printf0(" r2/in2[%i]: %g\n", i, r2[i]/in2[i]);
    }
#endif

    refine(out, in, r, r2, mass2, scale, nsolve, ineo, x, x, nm,
	   oi, xi, nsolve, tv);
    for(int i=0; i<nsolve; i++) {
      if(res_arg[i]->relmin>0) {
	res_arg[i]->final_rel = QOP_relnorm2_V(&tv, &out[i], insub, 1);
      } else {
	res_arg[i]->final_rel = 1;
      }
    }
  }

  inv_arg->max_iter = max_iter;

  for(int i=0; i<nsolve; i++) {
    QLA_Complex s;
    QLA_c_eq_r_times_c(s, masses[i], scale[i]);
    QDP_V_eq_c_times_V(tv, &s, out[i], insub);
    QDP_V_eq_V(out[i], tv, insub);
    res_arg[i]->final_iter = iter;
    res_arg[i]->final_rsq = r2[i]/in2[i];
    res_arg[i]->final_restart = nrestart;
    if(QOP_common.verbosity>=QOP_VERB_MED) {
      QLA_Real o2;
      QDP_r_eq_norm2_V(&o2, out[i], insub);
      QOP_printf0("SOLVE[%i]: in2: %g  out2: %g\n", i, in2[i], o2);
    }
  }

  double nflop = 0.5 * (20*QLA_Nc+2*8*QLA_Nc*QLA_Nc*fla->nlinks);
  double nflopm = 0.5 * 10*QLA_Nc; /* per extra mass */
  info->final_sec = dtime;
  info->final_flop = (nflop+nflopm*(nm-1))*iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  for(int i=0; i<nm; i++) {
    QDP_destroy_V(x[i]);
  }
  for(int i=0; i<nsolve; i++) {
    QDP_destroy_V(r[i]);
  }
  QDP_destroy_V(tv);
#undef NC
}

void
QOP_asqtad_invert_qdp(QOP_info_t *info,
		      QOP_FermionLinksAsqtad *fla,
		      QOP_invert_arg_t *inv_arg,
		      QOP_resid_arg_t *res_arg,
		      REAL mass,
		      QDP_ColorVector *out,
		      QDP_ColorVector *in)
{
  if(inv_arg->evenodd!=QOP_EVENODD) {
    QOP_asqtad_solve_multi_qdp(info, fla, inv_arg, &res_arg,
			       &mass, &out, &in, 1);
    return;
  }
#define NC QDP_get_nc(in)
  /* cg has 5 * 12 = 60 flops/site/it */
  /* MdagM -> 2*(66+72*15)+12 = 2304 flops/site */
  /* MdagM without naik -> 2*(66+72*7)+12 = 1152 flops/site */
  double dtime = 0;
  double nflop = 0.5 * (20*QLA_Nc+2*8*QLA_Nc*QLA_Nc*fla->nlinks);
  double rsqminold, relminold;
  QLA_Real rsq, rsqstop, relnorm2, insq, pinsq, rsqfac;
  QDP_ColorVector *qdpin, *qdpout;
  QDP_ColorVector *cgp, *cgr;
  QDP_Subset insub, cgsub;
  QOP_evenodd_t ineo, cgeo;
  int iter = 0;
  int max_iter_old = inv_arg->max_iter;
  int max_restarts_old = inv_arg->max_restarts;
  int nrestart = -1, max_restarts = inv_arg->max_restarts;
  if(max_restarts<=0) max_restarts = 5;

  ASQTAD_INVERT_BEGIN;

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

  // calculate input norm2 on even and odd subsets
  QLA_Real insqe=0, insqo=0;
  if(ineo!=QOP_ODD) { // even sites
    QDP_r_eq_norm2_V(&insqe, in, QDP_even);
  }
  if(ineo!=QOP_EVEN) { // odd sites
    QDP_r_eq_norm2_V(&insqo, in, QDP_odd);
  }
  insq = insqe + insqo;

  // run preconditioned solve on larger input subset
  if(insqe>=insqo) cgeo = QOP_EVEN;
  else cgeo = QOP_ODD;
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  // get temporary vectors
  cgp = QOP_asqtad_dslash_get_tmp(fla, oppsub(cgeo), 1);
  cgr = QOP_asqtad_dslash_get_tmp(fla, oppsub(cgeo), 2);
  gl_tmp = cgr;
  gl_tmp2 = QOP_asqtad_dslash_get_tmp(fla, cgeo, 1);

  // copy and project input vector
  QDP_V_eq_V(gl_tmp2, in, insub);
  if(ineo!=QOP_EVENODD) QDP_V_eq_zero(gl_tmp2, qdpsub(oppsub(ineo)));
  project(fla, mass, qdpin, gl_tmp2, cgeo);
  if(ineo==cgeo) {
    QDP_V_eq_zero(qdpin, qdpsub(oppsub(cgeo)));
  } else {
    QDP_V_eq_V(qdpin, in, qdpsub(oppsub(cgeo)));
  }
  // get projected input norm2
  QDP_r_eq_norm2_V(&pinsq, qdpin, cgsub);
  rsqfac = mass*mass*insq/pinsq;  // conversion factor inner rsq/outer rsq

  // copy output vector
  //QOP_asqtad_diaginv_qdp(NULL, fla, mass, cgp, qdpin, cgeo);
  //QDP_V_eq_V(qdpin, cgp, cgsub);
  //if(ineo!=QOP_EVENODD && ineo!=cgeo) {
  //QDP_V_eq_zero(qdpout, cgsub);
  //reconstruct(fla, mass, qdpout, out->cv, qdpout, oppsub(ineo));
  //}
  QDP_V_eq_zero(qdpout, QDP_all);
  QDP_V_eq_V(qdpout, out, cgsub);

  rsqstop = insq * res_arg->rsqmin;
  VERB(LOW, "ASQTAD_INVERT: rsqstop = %g\n", rsqstop);
  rsq = 0;
  relnorm2 = 1.;
  rsqminold = res_arg->rsqmin;
  relminold = res_arg->relmin;
  res_arg->rsqmin = 0.99 * rsqminold * rsqfac;
  res_arg->relmin = 0.99 * relminold;
  inv_arg->max_restarts = 0;
  do {
    dtime -= QOP_time();
    if(QOP_asqtad.cgtype==1) {
      QOP_invert_eigcg_V(QOP_asqtad_invert_d2, inv_arg, res_arg,
			 qdpout, qdpin, cgp, cgsub, &fla->eigcg);
    } else {
      QOP_invert_cg_V(QOP_asqtad_invert_d2, inv_arg, res_arg,
		      qdpout, qdpin, cgp, cgsub);
    }
    dtime += QOP_time();
    iter += res_arg->final_iter;

    reconstruct(fla, mass, qdpout, qdpout, qdpin, oppsub(cgeo));
    //QOP_asqtad_dslash_qdp(NULL, fla, mass, cgr, qdpout, oppsub(ineo), QOP_EVENODD);
    //QDP_r_eq_norm2_V(&rsq, cgr, qdpsub(oppsub(ineo)));
    //printf("0 ?= %g\n", rsq);

    // get final residual and relmin
    QOP_asqtad_dslash_qdp(NULL, fla, mass, cgr, qdpout, ineo, QOP_EVENODD);
    QDP_V_meq_V(cgr, in, insub);
    QDP_r_eq_norm2_V(&rsq, cgr, insub);
    res_arg->rsqmin *= 0.99;
    if(res_arg->relmin>0) {
      relnorm2 = QOP_relnorm2_V(&cgr, &qdpout, insub, 1);
      res_arg->relmin = 0.99*res_arg->final_rel/relnorm2 * relminold;
    }
    //printf("%i %i rsq = %g\tprec rsq = %g\trsqstop = %g\n", nrestart,
    //res_arg->final_iter, rsq, res_arg->final_rsq, rsqstop);
    VERB(LOW, "ASQTAD_INVERT: iter %i rsq = %g rel = %g\n", iter, rsq,
	 relnorm2);
    inv_arg->max_iter = max_iter_old - iter;
  } while( (++nrestart < max_restarts) &&
	   (inv_arg->max_iter>0) &&
	   (rsq > rsqstop) &&
	   (relnorm2 > res_arg->relmin) );

  QDP_V_eq_V(out, qdpout, insub);

  QDP_destroy_V(qdpin);
  QDP_destroy_V(qdpout);

  inv_arg->max_iter = max_iter_old;
  inv_arg->max_restarts = max_restarts_old;
  res_arg->rsqmin = rsqminold;
  res_arg->relmin = relminold;
  res_arg->final_iter = iter;
  res_arg->final_rsq = rsq/insq;
  res_arg->final_rel = relnorm2;
  res_arg->final_restart = nrestart;

  info->final_sec = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;

  //dumpvec(out, insub);

  ASQTAD_INVERT_END;
#undef NC
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
  QDP_ColorVector *in[nsrc], **out[nsrc];
  for(int i=0; i<nsrc; i++) {
    in[i] = in_pt[i]->cv;
    out[i] = (QDP_ColorVector **) malloc(nmass[i]*sizeof(QDP_ColorVector *));
    for(int j=0; j<nmass[i]; j++) {
      out[i][j] = out_pt[i][j]->cv;
    }
  }
  QOP_asqtad_invert_multi_qdp(info, fla, inv_arg, res_arg, masses, nmass,
			      out, in, nsrc);
  for(int i=0; i<nsrc; i++) {
    free(out[i]);
  }
}

void
QOP_asqtad_invert_multi_qdp(QOP_info_t *info,
			    QOP_FermionLinksAsqtad *fla,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t **res_arg[],
			    QOP_Real *masses[], int nmass[],
			    QDP_ColorVector **out_pt[],
			    QDP_ColorVector *in_pt[],
			    int nsrc)
{
  if(inv_arg->evenodd!=QOP_EVENODD) {
    int nsolve = 0;
    for(int i=0; i<nsrc; i++) nsolve += nmass[i];
    QOP_resid_arg_t *ra[nsolve];
    QOP_Real ms[nsolve];
    QDP_ColorVector *out[nsolve], *in[nsolve];
    int k = 0;
    for(int i=0; i<nsrc; i++) {
      for(int j=0; j<nmass[i]; j++) {
	ra[k] = res_arg[i][j];
	ms[k] = masses[i][j];
	out[k] = out_pt[i][j];
	in[k] = in_pt[i];
	k++;
      }
    }
    QOP_asqtad_solve_multi_qdp(info, fla, inv_arg, ra, ms, out, in, nsolve);
    return;
  }
#define NC QDP_get_nc(in_pt[0])
  /* cg has 5 * 12 = 60 flops/site/it */
  /* MdagM -> 2*(66+72*15)+12 = 2304 flops/site */
  /* MdagM without naik -> 2*(66+72*7)+12 = 1152 flops/site */
  double dtime;
  double nflop = 0.5 * (20*QLA_Nc+2*8*QLA_Nc*QLA_Nc*fla->nlinks);
  double nflopm = 0.5 * 10*QLA_Nc; /* per extra mass */
  QDP_ColorVector *cgp, *cgr;
  QDP_Subset cgsub;
  QOP_evenodd_t ineo, cgeo;
  int i, j;

  ASQTAD_INVERT_BEGIN;

  ineo = inv_arg->evenodd;

  //QOP_printf0("masses[0][0] = %g\n", masses[0][0]);
  //dumpvec(in_pt[0], insub);
  //QOP_asqtad_dslash_qdp(NULL,fla,masses[0][0],out_pt[0][0],in_pt[0],QOP_EVENODD,ineo);
  //dumpvec(out_pt[0][0], QDP_all);

  gl_fla = fla;

  cgeo = ineo;
  if(ineo==QOP_EVENODD) {
    cgeo = QOP_EVEN;
    printf("warning: QOP_asqtad_invert_multi on ALL subset not supported yet!\n");
  }
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  cgp = QOP_asqtad_dslash_get_tmp(fla, oppsub(cgeo), 1);
  cgr = QOP_asqtad_dslash_get_tmp(fla, oppsub(cgeo), 2);
  gl_tmp = cgr;
  gl_tmp2 = QOP_asqtad_dslash_get_tmp(fla, cgeo, 1);

  info->final_flop = 0;
  dtime = -QDP_time();

#if 0  // fake version
  for(i=0; i<nsrc; i++) {
    for(j=0; j<nmass[i]; j++) {
      gl_mass = masses[i][j];

      QOP_invert_cg_V(QOP_asqtad_invert_d2_norm2, inv_arg, res_arg[i][j],
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
      incv[i] = in_pt[i];
      outcv[i] = out_pt[i][0];
      ra[i] = res_arg[i][0];
    }
    gl_mass = masses[imin][0];
    for(i=0; i<2; i++) shifts[i] -= gl_mass*gl_mass;

    st = 1/(shifts[1]-shifts[0]);
    QDP_V_eq_V_minus_V(x0, incv[1], incv[0], cgsub);
    QDP_V_eq_r_times_V(x0, &st, x0, cgsub);
    QDP_V_eq_V(cgp, x0, cgsub);
    QOP_asqtad_invert_d2(src, cgp, cgsub);
    QDP_V_eq_V_minus_V(src, incv[imin], src, cgsub);

    QDP_r_eq_norm2_V(&rsq1, src, cgsub);
    for(i=0; i<2; i++) {
      QDP_r_eq_norm2_V(&rsq0, incv[i], cgsub);
      rsqsave[i] = res_arg[i][0]->rsqmin;
      res_arg[i][0]->rsqmin *= rsq0/rsq1;
    }

    QOP_invert_cgms_V(QOP_asqtad_invert_d2_norm2, inv_arg, ra, shifts, 2,
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
	cv[j] = out_pt[i][j];
      }
      gl_mass = masses[i][jmin];
      for(j=0; j<nmass[i]; j++) shifts[j] -= gl_mass*gl_mass;

      QOP_invert_cgms_V(QOP_asqtad_invert_d2_norm2, inv_arg, res_arg[i],
			shifts, nmass[i], cv, in_pt[i], cgp, cgsub);

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
      QDP_V_eq_r_times_V(out_pt[i][j], &m, out_pt[i][j], cgsub);
    }
  }

  info->final_sec = dtime;
  info->status = QOP_SUCCESS;

  //dumpvec(out_pt[0][0], insub);

  ASQTAD_INVERT_END;
#undef NC
}
