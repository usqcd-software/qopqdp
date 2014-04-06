#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static QDP_Lattice *lat=NULL;
static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=5;

//static const int sta[] = {0, 1};
//static const int sta[] = {1};
//static const int stn = sizeof(sta)/sizeof(int);
//static const int nsa[] = {2, 4, 8, 16};
//static const int nsn = sizeof(nsa)/sizeof(int);
//static const int nma[] = {2, 4, 8, 16};
//static const int nma[] = {0};
//static const int nmn = sizeof(nma)/sizeof(int);
static const int bsa[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
static const int bsn = sizeof(bsa)/sizeof(int);
static QDP_ColorMatrix **u, **cm;
static QOP_GaugeField *gauge;
static double flops, secs;
const double beta = 6.0;
double act=0, act2=0;
double plq=0, plq2=0;
int nact=0;

double
bench_action(QOP_gauge_coeffs_t *coeffs, QOP_Force *out)
{
  double sec=0, flop=0, mf=0;
  QLA_Real acts, actt;
  QOP_info_t info = QOP_INFO_ZERO;

  for(int i=0; i<=nit; i++) {
    QOP_symanzik_1loop_gauge_action(&info, gauge, &acts, &actt, coeffs);

    if(i>0) {
      sec += info.final_sec;
      flop += info.final_flop;
      mf += info.final_flop/(1e6*info.final_sec);
    }
  }

#if 1
  printf0("action s: %g  t: %g  tot: %g\n", acts, actt, acts+actt);

  coeffs->plaquette /= 4;
  coeffs->rectangle /= 6;
  coeffs->parallelogram /= 6;
  coeffs->adjoint_plaquette /= 8;
  QLA_Real eps=1;
  QOP_verbose(QOP_VERB_DEBUG);
  QOP_symanzik_1loop_gauge_force(&info, gauge, out, coeffs, eps);
  QOP_verbose(QOP_VERB_OFF);
  coeffs->plaquette *= 4;
  coeffs->rectangle *= 6;
  coeffs->parallelogram *= 6;
  coeffs->adjoint_plaquette *= 8;
#endif

  secs = sec/nit;
  flops = flop/nit;
  return mf/nit;
}

void
check_staple(QOP_gauge_coeffs_t *coeffs, QOP_Force *out)
{
  QOP_info_t info = QOP_INFO_ZERO;
  for(int i=0; i<ndim; i++) QDP_M_eq_zero(cm[i], QDP_all);
  QOP_symanzik_1loop_gauge_deriv_qdp(&info, u, cm, coeffs, -QLA_Nc, 0);

  int nsub=2;
  QDP_Subset *subs = QDP_even_and_odd_L(lat);
  int imp = (coeffs->rectangle!=0) || (coeffs->parallelogram!=0);
  if(imp) {
    nsub = 32;
    QDP_Subset *QOP_get_sub32(QDP_Lattice *);
    subs = QOP_get_sub32(lat);
  }
  QDP_ColorMatrix *st = QDP_create_M_L(lat);
  for(int i=0; i<ndim; i++) {
    QDP_M_eq_zero(st, QDP_all_L(lat));
    for(int subi=0; subi<nsub; subi++) {
      QOP_symanzik_1loop_gauge_staple_qdp(&info, u, st, i, coeffs,
					  subs, subi);
    }
    QLA_Real nrm2;
    QDP_M_meq_M(st, cm[i], QDP_all);
    QDP_r_eq_norm2_M(&nrm2, st, QDP_all);
    printf0("%i norm2: %g\n", i, nrm2);
  }
  QDP_destroy_M(st);
}

double
bench_heatbath(QOP_gauge_coeffs_t *coeffs)
{
  double sec=0, flop=0, mf=0;
  QOP_info_t info = QOP_INFO_ZERO;

  for(int i=0; i<=nit; i++) {
    QOP_symanzik_1loop_gauge_heatbath_qdp(&info, u, beta, coeffs, rs, 1, 1, 1);

    if(i>0) {
      sec += info.final_sec;
      flop += info.final_flop;
      mf += info.final_flop/(1e6*info.final_sec);
    }
    QLA_Real plaq = get_plaq(u);
    plq += plaq;
    plq2 += plaq*plaq;
    QLA_Real acts, actt;
    QOP_symanzik_1loop_gauge_action_qdp(&info, u, &acts, &actt, coeffs);
    double at = acts + actt;
    act += at;
    act2 += at*at;
    nact++;
  }

  secs = sec/nit;
  flops = flop/nit;
  return mf/nit;
}

void
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  int i, bs, bsi, best_bs;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<ndim; i++) u[i] = QDP_create_M();
  get_random_links(u, ndim, 0.3);

  plaq = get_plaq(u);
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);

  QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
  qoplayout.latdim = ndim;
  qoplayout.latsize = (int *) malloc(ndim*sizeof(int));
  for(i=0; i<ndim; i++) {
    qoplayout.latsize[i] = lattice_size[i];
  }
  qoplayout.machdim = -1;

  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_init(&qoplayout);

  gauge = QOP_create_G_from_qdp(u);

  QOP_Force *force;

  QDP_ColorMatrix *cm0[ndim];
  cm = cm0;
  for(i=0; i<ndim; i++) {
    cm[i] = QDP_create_M();
    QDP_M_eq_zero(cm[i], QDP_all);
  }

  QOP_gauge_coeffs_t gcoeffs = QOP_GAUGE_COEFFS_ZERO;
  gcoeffs.plaquette  = 0.4;
  gcoeffs.rectangle  = 0.3;
  gcoeffs.parallelogram   = 0.2;
  gcoeffs.adjoint_plaquette = 0;

  force = QOP_create_F_from_qdp(cm);
  mf = bench_action(&gcoeffs, force);
  QOP_destroy_F(force);
  printf0("action: sec%7.4f mflops = %g\n", secs, mf);

  force = QOP_create_F_from_qdp(cm);
  check_staple(&gcoeffs, force);
  QOP_destroy_F(force);
  gcoeffs.rectangle = 0;
  gcoeffs.parallelogram = 0;
  force = QOP_create_F_from_qdp(cm);
  check_staple(&gcoeffs, force);
  QOP_destroy_F(force);

  if(QDP_this_node==0) { printf("begin force\n"); fflush(stdout); }

  gcoeffs.plaquette  = 1;
  gcoeffs.rectangle  = 0.3;
  gcoeffs.parallelogram   = 0.2;
  best_mf = 0;
  best_bs = bsa[0];
  for(bsi=0; bsi<bsn; bsi++) {
    bs = bsa[bsi];
    QDP_set_block_size(bs);
    force = QOP_create_F_from_qdp(cm);
    mf = bench_heatbath(&gcoeffs);
    QOP_destroy_F(force);
    printf0("GF: bs%5i sec%7.4f mflops = %g\n", bs, secs, mf);
    if(mf>best_mf) {
      best_mf = mf;
      best_bs = bs;
    }
    //plaq = get_plaq(u);
    //if(QDP_this_node==0) printf("plaquette = %g\n", plaq);
  }

  if(QDP_this_node==0) {
    double pa = plq/nact;
    double p2 = plq2/nact;
    double pe = sqrt((p2-pa*pa)/(nact-1));
    printf("average plq: %g\tvar: %g\n", pa, pe);
    double aa = beta*act/nact;
    double a2 = beta*beta*act2/nact;
    double av = sqrt(a2-aa*aa);
    printf("average act: %g\tvar: %g\n", aa, av);
  }

  QDP_set_block_size(best_bs);
  QDP_profcontrol(1);
  force = QOP_create_F_from_qdp(cm);
  //mf = bench_staple(&gcoeffs, force);
  QDP_profcontrol(0);
  printf0("prof: GF: bs%5i sec%7.4f mflops = %g\n", best_bs, secs, mf);

  printf0("best: GF: bs%5i mflops = %g\n", best_bs, best_mf);

  if(QDP_this_node==0) { printf("begin unload links\n"); fflush(stdout); }
  //QOP_asqtad_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  QOP_finalize();
}

void
usage(char *s)
{
  printf("%s [n#] [s#] [x# [# ...]]\n",s);
  printf("\n");
  printf("n\tnumber of iterations\n");
  printf("s\tseed\n");
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
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 's' : seed=atoi(&argv[i][1]); break;
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
  lat = QDP_get_default_lattice();

  if(QDP_this_node==0) {
    printf("size = %i", lattice_size[0]);
    for(i=1; i<ndim; i++) {
      printf(" %i", lattice_size[i]);
    }
    printf("\n");
    printf("seed = %i\n", seed);
  }

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}
