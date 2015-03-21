#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=4;
static int nsrcmin=1;
static int nsrcmax=8;
static QLA_Real kappa=0.136;
static int clover=0;
static int verb=0;
#ifdef _OPENMP
static int bsmin=(512*1024*1024), bsmax=(512*1024*1024), *bsa, bsn;
#else
static int bsmin=32, bsmax=8192, *bsa, bsn;
#endif
static int check=0;

double
bench_force(QOP_info_t *info, QOP_FermionLinksWilson *flw,
	    QDP_ColorMatrix *cm[], QDP_DiracFermion *x[],
	    QDP_DiracFermion *y[], int nsrc)
{
  double sec=0, flop=0, mf=0;
  QLA_Real kappav[nsrc], eps[nsrc], sumeps;

  for(int i=0; i<ndim; i++) {
    QDP_M_eq_zero(cm[i], QDP_all);
  }
  //QOP_Force *force = QOP_create_F_from_qdp(cm);

  sumeps = 0;
  for(int i=0; i<nsrc; i++) {
    kappav[i] = kappa;
    eps[i] = 0.1*(i+1);
    sumeps += eps[i];
  }

  for(int i=0; i<=nit; i++) {
    QMP_barrier();
    //QOP_wilson_force_prec_multi_qdp(info, flw, force, kappav, eps, x, y, nsrc);
    QOP_wilson_force_prec_multi_qdp(info, flw, cm, kappav, eps, x, y, nsrc);
    QMP_barrier();
    if(i>0) {
      sec += info->final_sec;
      flop += info->final_flop;
    }
  }

  if(check) {
    //QOP_extract_F_to_qdp(cm, force);
    QLA_Real nrm2 = 0, tr = 0;
    QLA_ColorMatrix qcm, qcm2;
    for(int i=0; i<ndim; i++) {
      QLA_Real t;
      QDP_r_eq_norm2_M(&t, cm[i], QDP_all);
      nrm2 += t;
      QDP_m_eq_sum_M(&qcm, cm[i], QDP_all);
      QLA_M_eq_M_times_M(&qcm2, &qcm, &qcm);
      QLA_R_eq_re_trace_M(&t, &qcm2);
      tr += t;
    }
    printf0("|F|^2 = %12g   scaled = %g\n", nrm2, nrm2/(sumeps*sumeps));
    printf0("ReTr(sum^2F) = %12g   scaled = %g\n", tr, tr/(sumeps*sumeps));
  }
  //QOP_destroy_F(force);

  mf = 1;
  QMP_sum_double(&mf);
  QMP_sum_double(&sec);
  QMP_sum_double(&flop);
  info->final_sec = sec/(mf*nit);
  info->final_flop = flop/(mf*nit);
  mf = info->final_flop/(1e6*info->final_sec);
  return mf;
}

void
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  QDP_ColorMatrix *u[ndim];
  QDP_DiracFermion *x[nsrcmax], *y[nsrcmax];
  QOP_GaugeField *gf;
  QOP_wilson_coeffs_t coeffs = QOP_WILSON_COEFFS_ZERO;
  QOP_FermionLinksWilson *flw;
  int best_bs, best_nsrc;

  for(int i=0; i<ndim; i++) u[i] = QDP_create_M();
  get_random_links(u, ndim, 0.2);
  if(clover) {
    //clov = QDP_create_P();
    //QDP_P_eq_zero(clov, QDP_all);
    coeffs.clov_s = 0.1;
    coeffs.clov_t = 0.1;
  } else {
    //clov = NULL;
    coeffs.clov_s = 0;
    coeffs.clov_t = 0;
  }
  coeffs.aniso = 1;

  plaq = get_plaq(u);
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);

  QDP_ColorMatrix *cm[4];
  for(int i=0; i<4; i++) {
    cm[i] = QDP_create_M();
    QDP_M_eq_zero(cm[i], QDP_all);
  }

  for(int i=0; i<nsrcmax; i++) {
    x[i] = QDP_create_D();
    y[i] = QDP_create_D();
    QDP_D_eq_gaussian_S(x[i], rs, QDP_all);
    QDP_D_eq_gaussian_S(y[i], rs, QDP_all);
  }

  QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
  qoplayout.latdim = ndim;
  qoplayout.latsize = (int *) malloc(ndim*sizeof(int));
  for(int i=0; i<ndim; i++) {
    qoplayout.latsize[i] = lattice_size[i];
  }
  qoplayout.machdim = -1;

  QOP_info_t info = QOP_INFO_ZERO;
  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_init(&qoplayout);
  QOP_verbose(verb);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  //flw = QOP_wilson_create_L_from_qdp(u, clov);
  gf = QOP_create_G_from_qdp(u);
  flw = QOP_wilson_create_L_from_G(&info, &coeffs, gf);
  if(QDP_this_node==0) { printf("begin invert\n"); fflush(stdout); }

  best_mf = 0;
  best_nsrc = 0;
  best_bs = bsa[0];
  for(int nsrc=nsrcmin; nsrc<=nsrcmax; nsrc++) {
    for(int bsi=0; bsi<bsn; bsi++) {
      int bs = bsa[bsi];
      QDP_set_block_size(bs);
      mf = bench_force(&info, flw, cm, x, y, nsrc);
      printf0("FF: nsrc%3i bs%5i sec%7.4f mflops = %g\n", nsrc, bs, info.final_sec, mf);
      if(mf>best_mf) {
	best_mf = mf;
	best_nsrc = nsrc;
	best_bs = bs;
      }
    }
  }

  QDP_set_block_size(best_bs);
  QDP_profcontrol(1);
  mf = bench_force(&info, flw, cm, x, y, best_nsrc);
  QDP_profcontrol(0);
  printf0("prof: FF: nsrc%3i bs%5i sec%7.4f mflops = %g\n", best_nsrc, best_bs, info.final_sec, mf);
  printf0("best: FF: nsrc%3i bs%5i mflops = %g\n", best_nsrc, best_bs, best_mf);

  if(QDP_this_node==0) { printf("begin unload links\n"); fflush(stdout); }
  //QOP_wilson_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  QOP_finalize();
}

void
usage(char *s)
{
  printf("%s [n#] [s#] [S#] [x# [# ...]]\n",s);
  printf("\n");
  printf("b\tmin QDP blocksize\n");
  printf("B\tmax QDP blocksize\n");
  printf("k\tkappa\n");
  printf("n\tnumber of iterations\n");
  printf("s\tseed\n");
  printf("S\tnsrc\n");
  printf("v\tverbosity\n");
  printf("w\tclover term (0=off, 1=on)\n");
  printf("x\tlattice sizes (Lx, [Ly], ..)\n");
  printf("\n");
  exit(1);
}

int
main(int argc, char *argv[])
{
  int i, j;

  QDP_initialize(&argc, &argv);
  QDP_profcontrol(0);

  seed = time(NULL);
  j = 0;
  for(i=1; i<argc; i++) {
    switch(argv[i][0]) {
    case 'b' : bsmin=atoi(&argv[i][1]); break;
    case 'B' : bsmax=atoi(&argv[i][1]); break;
    case 'k' : kappa=atof(&argv[i][1]); break;
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 's' : seed=atoi(&argv[i][1]); break;
    case 'S' : nsrcmax=atoi(&argv[i][1]); break;
    case 'v' : verb=atoi(&argv[i][1]); break;
    case 'w' : clover=atoi(&argv[i][1]); break;
    case 'x' : j=i; while((i+1<argc)&&(isdigit(argv[i+1][0]))) ++i; break;
    default : usage(argv[0]);
    }
  }

  lattice_size = (int *) malloc(ndim*sizeof(int));
  if(j==0) {
    for(i=0; i<ndim; ++i) lattice_size[i] = 8;
  } else {
    if(!isdigit(argv[j][1])) usage(argv[0]);
    lattice_size[0] = atoi(&argv[j][1]);
    for(i=1; i<ndim; ++i) {
      if((++j<argc)&&(isdigit(argv[j][0]))) {
        lattice_size[i] = atoi(&argv[j][0]);
      } else {
        lattice_size[i] = lattice_size[i-1];
      }
    }
  }
  QDP_set_latsize(ndim, lattice_size);
  QDP_create_layout();

  for(i=0,j=bsmin; j<=bsmax; i++,j*=2) bsn = i+1;
  bsa = (int *) malloc(bsn*sizeof(*bsa));
  for(i=0,j=bsmin; j<=bsmax; i++,j*=2) bsa[i] = j;

  if(QDP_this_node==0) {
    print_layout();
    printf("kappa = %g\n", kappa);
    printf("seed = %i\n", seed);
    printf("clover = %i\n", clover);
  }

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}
