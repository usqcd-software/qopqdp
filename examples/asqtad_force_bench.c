#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=5;
static QLA_Real mass=-1;

static const int sta[] = {0, 1};
//static const int sta[] = {1};
static const int stn = sizeof(sta)/sizeof(int);
static const int nsa[] = {2, 4, 8, 16};
static const int nsn = sizeof(nsa)/sizeof(int);
static const int nma[] = {2, 4, 8, 16};
//static const int nma[] = {0};
static const int nmn = sizeof(nma)/sizeof(int);
static const int bsa[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
static const int bsn = sizeof(bsa)/sizeof(int);
static QOP_GaugeField *gauge;
static double flops, secs;

double
bench_force(QOP_asqtad_coeffs_t *asqcoef,
	    QOP_Force *out, QDP_ColorVector *in)
{
  double sec=0, flop=0, mf=0;
  int i;
  QOP_ColorVector *qopin[2];
  QLA_Real eps[2];
  QOP_info_t info;

  eps[0] = 0.1;
  eps[1] = 0.2;
  qopin[0] = QOP_create_V_from_qdp(in);
  qopin[1] = QOP_create_V_from_qdp(in);
  for(i=0; i<=nit; i++) {

    QOP_asqtad_force_multi(&info, gauge, out, asqcoef, eps, qopin, 2);

    if(i>0) {
      sec += info.final_sec;
      flop += info.final_flop;
      mf += info.final_flop/(1e6*info.final_sec);
    }
  }
  //printf("test1\n");
  QOP_destroy_V(qopin[0]);
  QOP_destroy_V(qopin[1]);
  //printf("test2\n");

  secs = sec/nit;
  flops = flop/nit;
  return mf/nit;
}

void
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  QDP_ColorMatrix **u;
  int i, bs, bsi, best_bs;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<ndim; i++) u[i] = QDP_create_M();
  get_random_links(u, ndim, 0.3);

  plaq = get_plaq(u);
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);

  QDP_ColorVector *in;
  in = QDP_create_V();
  QDP_V_eq_gaussian_S(in, rs, QDP_all);

  QOP_layout_t qoplayout;
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
  QDP_ColorMatrix *cm[4];
  for(i=0; i<4; i++) {
    cm[i] = QDP_create_M();
    QDP_M_eq_zero(cm[i], QDP_all);
  }
  force = QOP_create_F_from_qdp(cm);

  QOP_asqtad_coeffs_t asqcoef;
  asqcoef.one_link = 1;
  asqcoef.three_staple = 0.1;
  asqcoef.five_staple = 0.03;
  asqcoef.seven_staple = 0.01;
  asqcoef.lepage = 0.1;
  asqcoef.naik = 0.1;

  if(QDP_this_node==0) { printf("begin force\n"); fflush(stdout); }

  best_mf = 0;
  best_bs = bsa[0];
  for(bsi=0; bsi<bsn; bsi++) {
    bs = bsa[bsi];
    QDP_set_block_size(bs);
    mf = bench_force(&asqcoef, force, in);
    printf0("FF: bs%5i sec%7.4f mflops = %g\n", bs, secs, mf);
    if(mf>best_mf) {
      best_mf = mf;
      best_bs = bs;
    }
  }

  QDP_set_block_size(best_bs);
  QDP_profcontrol(1);
  mf = bench_force(&asqcoef, force, in);
  QDP_profcontrol(0);
  printf0("prof: FF: bs%5i sec%7.4f mflops = %g\n", best_bs, secs, mf);

  printf0("best: FF: bs%5i mflops = %g\n", best_bs, best_mf);

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
    case 'm' : mass=atof(&argv[i][1]); break;
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

  if(mass<0) {
    mass = 0.05 * pow(QDP_volume(),0.1);
  }

  if(QDP_this_node==0) {
    printf("size = %i", lattice_size[0]);
    for(i=1; i<ndim; i++) {
      printf(" %i", lattice_size[i]);
    }
    printf("\n");
    printf("mass = %g\n", mass);
    printf("seed = %i\n", seed);
  }

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}

