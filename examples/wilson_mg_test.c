#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int test_restart=0;

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=4;
static QLA_Real kappa=0.136;
static QLA_Real kappanv=0.136;
static double rsqmin=1e-8;
static int style=-1;
static int cgtype=-1;
static int clover=0;
static int verb=0;

//static const int sta[] = {0, 1, 2, 3};
//static const int sta[] = {1};
//static const int stn = sizeof(sta)/sizeof(int);
//static const int nsa[] = {1, 2, 4, 8};
//static const int nsn = sizeof(nsa)/sizeof(int);
//static const int nma[] = {1, 2, 4, 8};
//static const int nma[] = {0};
//static const int nmn = sizeof(nma)/sizeof(int);
//static const int bsa[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
//static const int bsn = sizeof(bsa)/sizeof(int);
static int bsmin=32, bsmax=8192, *bsa, bsn;
QOP_FermionLinksWilson *flw;
QOP_WilsonMg *wilmg;

double
bench(QOP_info_t *info, QOP_invert_arg_t *inv_arg,
      QOP_resid_arg_t *res_arg, QDP_DiracFermion *out,
      QDP_DiracFermion *in, int usemg)
{
  static QLA_Real r2s=-1, r2;
  double sec=0, flop=0, mf=0;
  int iter=0;

  for(int i=0; i<=nit; i++) {
    info->final_flop = 0;
    QDP_D_eq_zero(out, QDP_all);
    QMP_barrier();
    if(usemg) {
      QOP_wilsonMgSolve(info, wilmg, flw, inv_arg, res_arg, kappa, out, in);
    } else {
      QOP_wilson_invert_qdp(info, flw, inv_arg, res_arg, kappa, out, in);
    }
    QMP_barrier();
    if(i>0) {
      iter += res_arg->final_iter;
      sec += info->final_sec;
      flop += info->final_flop;
    }
  }
  QDP_r_eq_norm2_D(&r2, out, QDP_even);
  if(r2s<0) r2s = r2;
  if(fabs(1-r2/r2s)>1e-3) {
    printf0("first norm = %g  this norm = %g\n", r2s, r2);
  }
  mf = 1;
  QMP_sum_double(&mf);
  QMP_sum_double(&sec);
  QMP_sum_double(&flop);
  res_arg->final_iter = iter/nit;
  info->final_sec = sec/(mf*nit);
  info->final_flop = flop/(mf*nit);
  mf = info->final_flop/(1e6*info->final_sec);
  return mf;
}

void
mg_setup(void)
{
#define BLOCK0 {2,2,2,2}
#define MAXITER 2000
#define MGMAXITER 50
#define FLOATINNER 1
#define SETUPRES 0.4
#define SETUPMAXIT 100
#define SETUPCHANGEFAC 0.5
#define NVECS0 8

  QOP_wilsonMgSet(wilmg, -1, "nlevels", 1);
  QOP_wilsonMgSet(wilmg, -1, "verbose", verb);
  QOP_wilsonMgSet(wilmg, -1, "profile", 1);
  QOP_wilsonMgSet(wilmg, -2, "verbose", verb-1);
  QOP_wilsonMgSet(wilmg, -1, "kappa", kappa);
  QOP_wilsonMgSet(wilmg, -1, "kappanv", kappanv);
  QOP_wilsonMgSet(wilmg, -1, "itmax", MGMAXITER);
  QOP_wilsonMgSet(wilmg, -1, "ngcr", 8);
  QOP_wilsonMgSet(wilmg, -1, "nc", QLA_Nc);

  double dlat[4], vol=1, vol0=1;
  int lat[4], lat0[4], bl0[4]=BLOCK0;
  int nvecs0;
  QDP_latsize(lat);
  for(int i=0; i<4; i++) {
    if(i==3) {
      double fac = lat[0]/lat0[0];
      lat0[i] = lat[i];
      while(fac>=2 && lat0[i]>1) {
        int t = 2;
        if(lat0[i]%3==0) t = 3;
        lat0[i] /= t;
        fac /= t;
      }
    } else {
      if(lat[i]%3==0) lat0[i] = lat[i]/3;
      else lat0[i] = lat[i]/2;
    }
    lat0[i] = lat[i]/bl0[i];
    vol *= lat[i];
    vol0 *= lat0[i];
  }
  printf0("vol = %i\n", (int)vol);
  printf0("vol0 = %i\n", (int)vol0);
  //nvecs0 = floor(0.99+sqrt(8.*vol/(3.*vol0)))-3;
  nvecs0 = NVECS0;

  for(int i=0; i<4; i++) dlat[i] = lat0[i];
  QOP_wilsonMgSetArray(wilmg, 0, "lattice", dlat, 4);
  QOP_wilsonMgSet(wilmg, 0, "nvecs", nvecs0);
  //QOP_wilsonMgSetArray(wilmg, 0, "lattice", (double[]){12,12,12,32}, 4);
  //QOP_wilsonMgSetArray(wilmg, 0, "lattice", (double[]){8,8,8,8}, 4);
  //QOP_wilsonMgSet(wilmg, 0, "nvecs", 28);
  printf0("nvecs0 = %i\n", nvecs0);

  int npre = 0;
  int npost = 4;

  QOP_wilsonMgSet(wilmg, 0, "setup_res", SETUPRES);
  QOP_wilsonMgSet(wilmg, 0, "setup_maxit", SETUPMAXIT);
  QOP_wilsonMgSet(wilmg, 0, "setup_change_fac", SETUPCHANGEFAC);
  //QOP_wilsonMgSet(wilmg, 0, "setup_nvecs", nvecs0 + SETUPNVECSEXTRA);

  QOP_wilsonMgSet(wilmg, 0, "npre", npre);
  QOP_wilsonMgSet(wilmg, 0, "npost", npost);
  QOP_wilsonMgSet(wilmg, 0, "scale", 1);
  QOP_wilsonMgSet(wilmg, 0, "cres", 0.1);
  QOP_wilsonMgSet(wilmg, 0, "itmax", MGMAXITER);
  QOP_wilsonMgSet(wilmg, 0, "ngcr", 8);  // can't later be set to more than this

  printf0("\nsetup Wilson MG\n");
  double t0 = QDP_time();
  //#if QOP_Precision == 'F'
  QOP_wilsonMgSetLinks(wilmg, flw);
  //#else
  //QOP_F3_FermionLinksWilson *fflw = QOP_FD3_wilson_create_L_from_L(flw);
  //QOP_wilsonMgSetLinks(wilmg, fflw);
  //#endif
  QOP_wilsonMgSetup(wilmg);
  t0 = QDP_time() - t0;
  printf0("setup time: %g seconds\n", t0);
}

void
start(void)
{
  double mf;
  QLA_Real plaq;
  QDP_ColorMatrix **u;
  //QDP_DiracPropagator *clov;
  QDP_DiracFermion *out, *in;
  QOP_GaugeField *gf;
  QOP_wilson_coeffs_t coeffs = QOP_WILSON_COEFFS_ZERO;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
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

  out = QDP_create_D();
  in = QDP_create_D();
  QDP_D_eq_gaussian_S(in, rs, QDP_all);

  QOP_layout_t qoplayout = QOP_LAYOUT_ZERO;
  qoplayout.latdim = ndim;
  qoplayout.latsize = (int *) malloc(ndim*sizeof(int));
  for(int i=0; i<ndim; i++) {
    qoplayout.latsize[i] = lattice_size[i];
  }
  qoplayout.machdim = -1;

  QOP_info_t info = QOP_INFO_ZERO;
  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  res_arg.rsqmin = rsqmin;
  //res_arg.relmin = rsqmin;
  inv_arg.max_iter = 600;
  inv_arg.restart = 200;
  inv_arg.max_restarts = 5;
  inv_arg.evenodd = QOP_EVENODD;

  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_init(&qoplayout);
  QOP_verbose(verb);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  //flw = QOP_wilson_create_L_from_qdp(u, clov);
  gf = QOP_create_G_from_qdp(u);
  flw = QOP_wilson_create_L_from_G(&info, &coeffs, gf);
  if(QDP_this_node==0) { printf("begin invert\n"); fflush(stdout); }
  //QDP_set_block_size(bs);

  printf0("\ncreate Wilson MG\n");
  wilmg = QOP_wilsonMgNew();
  mg_setup();

  QDP_D_eq_gaussian_S(in, rs, QDP_all);

  mf = bench(&info, &inv_arg, &res_arg, out, in, 0);
  printf0("CG: iter%5i sec%7.4f mflops = %g\n",
	  res_arg.final_iter, info.final_sec, mf);
  mf = bench(&info, &inv_arg, &res_arg, out, in, 1);
  printf0("MG: iter%5i sec%7.4f mflops = %g\n",
	  res_arg.final_iter, info.final_sec, mf);

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
  printf("c\tcgtype (0=CGNE, 1=BiCGStab)\n");
  printf("k\tkappa\n");
  printf("n\tnumber of iterations\n");
  printf("r\ttest restart\n");
  printf("s\tseed\n");
  printf("S\tstyle\n");
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
    case 'c' : cgtype=atoi(&argv[i][1]); break;
    case 'k' : kappa=atof(&argv[i][1]); break;
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 'r' : test_restart=atoi(&argv[i][1]); break;
    case 's' : seed=atoi(&argv[i][1]); break;
    case 'S' : style=atoi(&argv[i][1]); break;
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

  kappanv = kappa;
  if(QDP_this_node==0) {
    print_layout();
    printf("kappa = %g\n", kappa);
    printf("kappanv = %g\n", kappanv);
    printf("seed = %i\n", seed);
    printf("cgtype = %i\n", cgtype);
    printf("clover = %i\n", clover);
  }

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}
