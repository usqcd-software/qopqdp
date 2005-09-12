#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=5;
static double kappa=0;

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

double
bench_inv(QOP_invert_arg *inv_arg, QDP_DiracFermion *out, QDP_DiracFermion *in)
{
  double sec=0, flop=0, mf=0;
  int i, iter=0;

  for(i=0; i<=nit; i++) {
    int invit;
    QDP_D_eq_zero(out, QDP_all);
    invit = QOP_F_wilson_inv_qdp(inv_arg, out, in);
    if(i>0) {
      iter += inv_arg->final_iter;
      sec += inv_arg->final_sec;
      flop += inv_arg->final_flop;
      mf += inv_arg->final_flop/(1e6*inv_arg->final_sec);
    }
  }
  inv_arg->final_iter = iter/nit;
  inv_arg->final_sec = sec/nit;
  inv_arg->final_flop = flop/nit;
  return mf/nit;
}

void
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  QDP_ColorMatrix **u;
  int i, st, ns, nm, bs, sti, nsi, nmi, bsi,
    best_st, best_ns, best_nm, best_bs;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<ndim; i++) u[i] = QDP_create_M();
  get_random_links(u, ndim, 0.2);

  plaq = get_plaq(u);
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);

  QDP_ColorMatrix *fatlinks[4], *longlinks[4];
  QDP_DiracFermion *out, *in;
  out = QDP_create_D();
  in = QDP_create_D();
  for(i=0; i<4; i++) {
    fatlinks[i] = QDP_create_M();
    QDP_M_eq_M(fatlinks[i], u[i], QDP_all);
    longlinks[i] = QDP_create_M();
    QDP_M_eq_M(longlinks[i], u[i], QDP_all);
  }
  QDP_D_eq_gaussian_S(in, rs, QDP_all);

  QOP_layout qoplayout;
  qoplayout.ndims = ndim;
  for(i=0; i<ndim; i++) {
    qoplayout.sites[i] = lattice_size[i];
    qoplayout.bc[i] = 1;
  }

  QOP_invert_arg inv_arg;
  inv_arg.mass = kappa;
  inv_arg.rsqmin = 1e-4;
  inv_arg.max_iter = 600;
  inv_arg.restart = 200;
  inv_arg.evenodd = QOP_EVEN;

  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_wilson_invert_init(&qoplayout);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  QOP_wilson_invert_load_links_qdp(u);
  if(QDP_this_node==0) { printf("begin invert\n"); fflush(stdout); }

  best_mf = 0;
  best_st = sta[0];
  best_ns = nsa[0];
  best_nm = nma[0];
  best_bs = bsa[0];
  for(sti=0; sti<stn; sti++) {
    st = sta[sti];
    if(QOP_wilson_set_opt("st", st)==QOP_FAIL) continue;
    for(nsi=0; nsi<nsn; nsi++) {
      ns = nsa[nsi];
      if(QOP_wilson_set_opt("ns", ns)==QOP_FAIL) continue;
      for(nmi=0; nmi<nmn; nmi++) {
	nm = nma[nmi];
	if(nm==0) nm = ns;
	if(QOP_wilson_set_opt("nm", nm)==QOP_FAIL) continue;
	for(bsi=0; bsi<bsn; bsi++) {
	  bs = bsa[bsi];
	  QDP_set_block_size(bs);
	  mf = bench_inv(&inv_arg, out, in);
	  printf0("CONGRAD: st%2i ns%2i nm%2i bs%5i iter%5i sec%7.4f mflops = %g\n", st,
		  ns, nm, bs, inv_arg.final_iter, inv_arg.final_sec, mf);
	  if(mf>best_mf) {
	    best_mf = mf;
	    best_st = st;
	    best_ns = ns;
	    best_nm = nm;
	    best_bs = bs;
	  }
	}
      }
    }
  }

  printf0("best:\n");
  printf0("CONGRAD: st%2i ns%2i nm%2i bs%5i mflops = %g\n", best_st, best_ns,
	  best_nm, best_bs, best_mf);

  if(QDP_this_node==0) { printf("begin unload links\n"); fflush(stdout); }
  QOP_F_wilson_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  QOP_wilson_invert_finalize();
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

  seed = time(NULL);
  j = 0;
  for(i=1; i<argc; i++) {
    switch(argv[i][0]) {
    case 'm' : kappa=atof(&argv[i][1]); break;
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

  if(QDP_this_node==0) {
    printf("size = %i", lattice_size[0]);
    for(i=1; i<ndim; i++) {
      printf(" %i", lattice_size[i]);
    }
    printf("\n");
    printf("seed = %i\n", seed);
  }

  QDP_set_latsize(ndim, lattice_size);
  QDP_create_layout();

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}
