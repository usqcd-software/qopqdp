#include <test_common.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=5;
static QLA_Real mass=-1;
static int style=-1;
static int nmass=1;
static int cgtype=0;
static int verb=0;

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
QOP_FermionLinksAsqtad *fla;

double
bench_inv(QOP_info_t *info, QOP_invert_arg_t *inv_arg,
	  QOP_resid_arg_t *res_arg, QDP_ColorVector *out, QDP_ColorVector *in)
{
  double sec=0, flop=0, mf=0;
  int i, j, iter=0;
  QOP_ColorVector *qopout[nmass], **pqo, *qopin;
  QLA_Real masses[nmass], *pm;
  QOP_resid_arg_t *ra[nmass], **pra;

  pqo = qopout;
  pm = masses;
  pra = ra;
  for(i=0; i<=nit; i++) {
    QDP_V_eq_zero(out, QDP_all);
    qopin = QOP_create_V_from_qdp(in);
    for(j=0; j<nmass; j++) {
      qopout[j] = QOP_create_V_from_qdp(out);
      masses[j] = mass*(j+1);
      ra[j] = res_arg;
    }
    QMP_barrier();
    if(nmass == 1) {
      QOP_asqtad_invert(info, fla, inv_arg, res_arg, mass, qopout[0], qopin);
    } else {
      QOP_asqtad_invert_multi(info, fla, inv_arg, &pra, &pm, &nmass, &pqo, &qopin, 1);
    }
    QMP_barrier();
    if(i>0) {
      iter += res_arg->final_iter;
      sec += info->final_sec;
      flop += info->final_flop;
      //mf += info->final_flop/(1e6*info->final_sec);
    }
    for(j=0; j<nmass; j++) {
      QOP_destroy_V(qopout[j]);
    }
    QOP_destroy_V(qopin);
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
  int i, st, ns, nm, bs, sti, nsi, nmi, bsi,
    best_st, best_ns, best_nm, best_bs;

  u = (QDP_ColorMatrix **) malloc(ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<ndim; i++) u[i] = QDP_create_M();
  get_random_links(u, ndim, 0.3);

  plaq = get_plaq(u);
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);

  QDP_ColorMatrix *fatlinks[4], *longlinks[4];
  QDP_ColorVector *out, *in;
  out = QDP_create_V();
  in = QDP_create_V();
  for(i=0; i<4; i++) {
    fatlinks[i] = QDP_create_M();
    QDP_M_eq_M(fatlinks[i], u[i], QDP_all);
    longlinks[i] = QDP_create_M();
    QDP_M_eq_M(longlinks[i], u[i], QDP_all);
  }
  QDP_V_eq_gaussian_S(in, rs, QDP_all);

  QOP_layout_t qoplayout;
  qoplayout.latdim = ndim;
  qoplayout.latsize = (int *) malloc(ndim*sizeof(int));
  for(i=0; i<ndim; i++) {
    qoplayout.latsize[i] = lattice_size[i];
  }
  qoplayout.machdim = -1;

  QOP_GaugeField *gf;
  QOP_asqtad_coeffs_t coeffs;
  coeffs.one_link = 1;
  coeffs.three_staple = 0.1;
  coeffs.five_staple = 0.1;
  coeffs.seven_staple = 0.1;
  coeffs.lepage = 0.1;
  coeffs.naik = 0.1;

  QOP_info_t info;
  QOP_invert_arg_t inv_arg = QOP_INVERT_ARG_DEFAULT;
  QOP_resid_arg_t res_arg = QOP_RESID_ARG_DEFAULT;
  res_arg.rsqmin = 1e-4;
  inv_arg.max_iter = 600;
  inv_arg.restart = 200;
  inv_arg.max_restarts = 5;
  inv_arg.evenodd = QOP_EVEN;

  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_init(&qoplayout);
  QOP_verbose(verb);
  QOP_opt_t optcg;
  optcg.tag = "cg";
  optcg.value = cgtype;
  QOP_asqtad_invert_set_opts(&optcg, 1);
  if(QDP_this_node==0) { printf("convert gauge field\n"); fflush(stdout); }
  gf = QOP_convert_G_from_qdp(u);
  if(QDP_this_node==0) { printf("rephase gauge field\n"); fflush(stdout); }
  QOP_Complex phase[4] = {{1,0},{1,0},{1,0},{-1,0}};
  int signmask[4] = {0,1,3,7};
  QOP_bc_t bc;
  QOP_staggered_sign_t ss;
  bc.phase = phase;
  ss.signmask = signmask;
  QOP_rephase_G(gf, &bc, &ss);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  //fla = QOP_asqtad_create_L_from_qdp(fatlinks, longlinks);
  QDP_profcontrol(1);
  fla = QOP_asqtad_create_L_from_G(&info, &coeffs, gf);
  QDP_profcontrol(0);
  if(QDP_this_node==0) { printf("load links: secs = %g\t mflops = %g\n", info.final_sec, info.final_flop/(1e6*info.final_sec)); }
  if(QDP_this_node==0) { printf("begin invert\n"); fflush(stdout); }
  //QOP_asqtad_destroy_L(fla);

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
    if(QOP_asqtad_invert_set_opts(&optst, 1)==QOP_FAIL) continue;
    for(nsi=0; nsi<nsn; nsi++) {
      ns = nsa[nsi];
      optns.value = ns;
      if(QOP_asqtad_invert_set_opts(&optns, 1)==QOP_FAIL) continue;
      for(nmi=0; nmi<nmn; nmi++) {
	nm = nma[nmi];
	if(nm==0) nm = ns;
	optnm.value = nm;
	if(QOP_asqtad_invert_set_opts(&optnm, 1)==QOP_FAIL) continue;
	for(bsi=0; bsi<bsn; bsi++) {
	  bs = bsa[bsi];
	  QDP_set_block_size(bs);
	  //fla = QOP_asqtad_create_L_from_G(&info, &coeffs, gf);
	  mf = bench_inv(&info, &inv_arg, &res_arg, out, in);
	  //QOP_asqtad_destroy_L(fla);
	  printf0("CONGRAD: st%2i ns%2i nm%2i bs%5i iter%5i sec%7.4f mflops = %g\n", st,
		  ns, nm, bs, res_arg.final_iter, info.final_sec, mf);
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
  //fla = QOP_asqtad_create_L_from_G(&info, &coeffs, gf);

  optst.value = best_st;
  optns.value = best_ns;
  optnm.value = best_nm;
  QOP_asqtad_invert_set_opts(&optst, 1);
  QOP_asqtad_invert_set_opts(&optns, 1);
  QOP_asqtad_invert_set_opts(&optnm, 1);
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
  //QOP_asqtad_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  //QOP_asqtad_invert_finalize();
}

void
usage(char *s)
{
  printf("%s [c#] [m#] [M#] [n#] [s#] [S#] [x# [# ...]] [v#]\n",s);
  printf("\n");
  printf("c\tcgtype\n");
  printf("m\tlightest mass\n");
  printf("M\tnumber of masses\n");
  printf("n\tnumber of iterations\n");
  printf("s\tseed\n");
  printf("S\tstyle\n");
  printf("x\tlattice sizes (Lx, [Ly], ..)\n");
  printf("v\tverbose\n");
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
    case 'c' : cgtype=atoi(&argv[i][1]); break;
    case 'm' : mass=atof(&argv[i][1]); break;
    case 'M' : nmass=atof(&argv[i][1]); break;
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 's' : seed=atoi(&argv[i][1]); break;
    case 'S' : style=atoi(&argv[i][1]); break;
    case 'x' : j=i; while((i+1<argc)&&(isdigit(argv[i+1][0]))) ++i; break;
    case 'v' : verb=atoi(&argv[i][1]); break;
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
    mass = 0.2 * pow(QDP_volume(),0.1);
  }

  if(QDP_this_node==0) {
    print_layout();
    printf("mass = %g\n", mass);
    printf("seed = %i\n", seed);
  }

  rs = QDP_create_S();
  seed_rand(rs, seed);

  start();

  QDP_finalize();
  return 0;
}
