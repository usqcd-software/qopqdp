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
static double rsqmin=1e-8;
static int style=-1;
static int cgtype=-1;
static int clover=0;
static int verb=0;

static const int sta[] = {0, 1, 2, 3};
//static const int sta[] = {1};
static const int stn = sizeof(sta)/sizeof(int);
static const int nsa[] = {1, 2, 4, 8};
static const int nsn = sizeof(nsa)/sizeof(int);
static const int nma[] = {1, 2, 4, 8};
//static const int nma[] = {0};
static const int nmn = sizeof(nma)/sizeof(int);
//static const int bsa[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
//static const int bsn = sizeof(bsa)/sizeof(int);
#ifdef _OPENMP
static int bsmin=(512*1024*1024), bsmax=(512*1024*1024), *bsa, bsn;
#else
static int bsmin=32, bsmax=8192, *bsa, bsn;
#endif
QOP_FermionLinksWilson *flw;

double
bench_inv(QOP_info_t *info, QOP_invert_arg_t *inv_arg,
	  QOP_resid_arg_t *res_arg, QDP_DiracFermion *out,
	  QDP_DiracFermion *in)
{
  static QLA_Real r2s=-1, r2;
  double sec=0, flop=0, mf=0;
  int i, iter=0;
  QOP_DiracFermion *qopout=NULL, *qopin=NULL;

  //QDP_D_eq_gaussian_S(in, rs, QDP_all);
  QDP_D_eq_zero(out, QDP_all);
  for(i=0; i<=nit; i++) {
    if(i==0 || !test_restart) {
      qopout = QOP_create_D_from_qdp(out);
      qopin = QOP_create_D_from_qdp(in);
    }
    QMP_barrier();
    QOP_wilson_invert(info, flw, inv_arg, res_arg, kappa, qopout, qopin);
    QMP_barrier();
    if(i>0) {
      iter += res_arg->final_iter;
      sec += info->final_sec;
      flop += info->final_flop;
      //mf += info->final_flop/(1e6*info->final_sec);
    }
    if(i==nit || !test_restart) {
      QOP_destroy_D(qopout);
      QOP_destroy_D(qopin);
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
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  QDP_ColorMatrix **u;
  //QDP_DiracPropagator *clov;
  QDP_DiracFermion *out, *in;
  QOP_GaugeField *gf;
  QOP_wilson_coeffs_t coeffs = QOP_WILSON_COEFFS_ZERO;
  int i, st, ns, nm, bs, sti, nsi, nmi, bsi,
    best_st, best_ns, best_nm, best_bs;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<ndim; i++) u[i] = QDP_create_M();
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
  for(i=0; i<ndim; i++) {
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

  if(cgtype>=0) {
    QOP_opt_t optcg;
    optcg.tag = "cg";
    optcg.value = cgtype;
    QOP_wilson_invert_set_opts(&optcg, 1);
  }
  QDP_D_eq_gaussian_S(in, rs, QDP_all);

  best_mf = 0;
  best_st = sta[0];
  best_ns = nsa[0];
  best_nm = nma[0];
  best_bs = bsa[0];
  QOP_opt_t optst;
  optst.tag = "st";
  QOP_opt_t optns;
  optns.tag = "ns";
  QOP_opt_t optnm;
  optnm.tag = "nm";
  for(sti=0; sti<stn; sti++) {
    if((style>=0)&&(sti!=style)) continue;
    st = sta[sti];
    optst.value = st;
    if(QOP_wilson_invert_set_opts(&optst, 1)==QOP_FAIL) continue;
    for(nsi=0; nsi<nsn; nsi++) {
      ns = nsa[nsi];
      optns.value = ns;
      if(QOP_wilson_invert_set_opts(&optns, 1)==QOP_FAIL) continue;
      for(nmi=0; nmi<nmn; nmi++) {
	nm = nma[nmi];
	if(nm==0) nm = ns;
	optnm.value = nm;
	if(QOP_wilson_invert_set_opts(&optnm, 1)==QOP_FAIL) continue;
	for(bsi=0; bsi<bsn; bsi++) {
	  bs = bsa[bsi];
	  QDP_set_block_size(bs);
	  mf = bench_inv(&info, &inv_arg, &res_arg, out, in);
	  printf0("CONGRAD: st%2i ns%2i nm%2i bs%5i iter%5i sec%7.4f mflops = %g\n", st,
		  ns, nm, bs, res_arg.final_iter, info.final_sec, mf);
	  fflush(stdout);
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

  optst.value = best_st;
  optns.value = best_ns;
  optnm.value = best_nm;
  QOP_wilson_invert_set_opts(&optst, 1);
  QOP_wilson_invert_set_opts(&optns, 1);
  QOP_wilson_invert_set_opts(&optnm, 1);
  QDP_set_block_size(best_bs);
  QDP_profcontrol(1);
  mf = bench_inv(&info, &inv_arg, &res_arg, out, in);
  QDP_profcontrol(0);
  printf0("prof: CONGRAD: st%2i ns%2i nm%2i bs%5i iter%5i sec%7.4f mflops = %g\n",
          best_st, best_ns, best_nm, best_bs,
          res_arg.final_iter, info.final_sec, mf);

  printf0("best: CONGRAD: st%2i ns%2i nm%2i bs%5i mflops = %g\n",
          best_st, best_ns, best_nm, best_bs, best_mf);

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

  if(QDP_this_node==0) {
    print_layout();
    printf("kappa = %g\n", kappa);
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
