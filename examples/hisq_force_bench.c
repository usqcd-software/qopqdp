#include <test_common.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

static int ndim=4;
static int *lattice_size;
static int seed;
static int nit=5;
static int nsrcmin=1;
static int nsrcmax=8;
static int check=0;

//static const int fsma[] = {MAX_NSRC};
//static const int fsmn = sizeof(fsma)/sizeof(int);
//static const int bsa[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
//static const int bsn = sizeof(bsa)/sizeof(int);
static int *fsma, fsmn=-1;
static int bsmin=32, bsmax=8192, *bsa, bsn;
static QOP_GaugeField *gauge;
static double flops, secs;

double
bench_force(QOP_FermionLinksHisq *flh, QOP_hisq_coeffs_t *coeffs,
	    QDP_ColorMatrix *cm[], QDP_ColorVector *in, int nsrc)
{
  double sec=0, flop=0, mf=0;
  QOP_ColorVector *qopin[nsrc];
  QLA_Real eps[nsrc], sumeps;
  QOP_info_t info;

  for(int i=0; i<ndim; i++) {
    QDP_M_eq_zero(cm[i], QDP_all);
  }
  QOP_Force *force = QOP_create_F_from_qdp(cm);

  sumeps = 0;
  for(int i=0; i<nsrc; i++) {
    eps[i] = 0.1*(i+1);
    sumeps += eps[i];
    qopin[i] = QOP_create_V_from_qdp(in);
  }

  for(int i=0; i<=nit; i++) {
    QOP_hisq_force_multi(&info, flh, force, coeffs, eps, qopin, &nsrc);
    if(i>0) {
      sec += info.final_sec;
      flop += info.final_flop;
      mf += info.final_flop/(1e6*info.final_sec);
    }
  }

  if(check) {
    QOP_extract_F_to_qdp(cm, force);
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
  QOP_destroy_F(force);
  for(int i=0; i<nsrc; i++) {
    QOP_destroy_V(qopin[i]);
  }

  secs = sec/nit;
  flops = flop/nit;
  printf("flops: %.0f\n", flops);
  return mf/nit;
}

void
start(void)
{
  double mf, best_mf;
  QLA_Real plaq;
  QDP_ColorMatrix **u;
  int i, fsm, fsmi, nsrc, bs, bsi, best_fsm, best_nsrc, best_bs;

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

  QDP_ColorMatrix *cm[4];
  for(i=0; i<4; i++) {
    cm[i] = QDP_create_M();
    QDP_M_eq_zero(cm[i], QDP_all);
  }

  QOP_hisq_coeffs_t coeffs;
  coeffs.n_naiks = 1;
  coeffs.eps_naik[0] = 0;
  //coeffs.eps_naik[1] = 0.1;
  coeffs.ugroup = QOP_UNITARIZE_U3;
  //coeffs.ugroup = QOP_UNITARIZE_SU3;
  coeffs.umethod = QOP_UNITARIZE_RATIONAL;
  //coeffs.umethod = QOP_UNITARIZE_ANALYTIC;
  coeffs.fat7_one_link = 1;
  coeffs.fat7_three_staple = 0.1;
  coeffs.fat7_five_staple = 0.1;
  coeffs.fat7_seven_staple = 0.1;
  coeffs.asqtad_one_link = 1;
  coeffs.asqtad_three_staple = 0.1;
  coeffs.asqtad_five_staple = 0.1;
  coeffs.asqtad_seven_staple = 0.1;
  coeffs.asqtad_lepage = 0.1;
  coeffs.asqtad_naik = 0.1;
  coeffs.difference_one_link = 1;
  coeffs.difference_naik = 1;
#if 0
  //coeffs.asqtad_one_link = 0;
  //coeffs.asqtad_three_staple = 0;
  //coeffs.asqtad_five_staple = 0;
  //coeffs.asqtad_seven_staple = 0;
  coeffs.asqtad_lepage = 0;
  //coeffs.asqtad_naik = 0;
#endif

  QOP_info_t info;
  QOP_FermionLinksHisq *flh;
  //QOP_verbose(verb);
  //if(QDP_this_node==0) { printf("convert gauge field\n"); fflush(stdout); }
  //gf = QOP_convert_G_from_qdp(u);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  //fla = QOP_asqtad_create_L_from_qdp(fatlinks, longlinks);
  QDP_profcontrol(1);
  //{ QOP_opt_t ol = {.tag="reunit_allow_svd_only",.value=0}; QOP_hisq_links_set_opts(&ol,1); }
  //{ QOP_opt_t ol = {.tag="reunit_svd_only",.value=1}; QOP_hisq_links_set_opts(&ol,1); }
  flh = QOP_hisq_create_L_from_G(&info, &coeffs, gauge);
  //fla = QOP_get_asqtad_links_from_hisq(flh)[0];
  QDP_profcontrol(0);
  if(QDP_this_node==0) { printf("load links: secs = %g\t mflops = %g\n", info.final_sec, info.final_flop/(1e6*info.final_sec)); }
  if(QDP_this_node==0) { printf("begin force\n"); fflush(stdout); }

  best_mf = -1;
  best_fsm = fsma[0];
  best_nsrc = 0;
  best_bs = bsa[0];
  QOP_opt_t optfsm;
  optfsm.tag = "fnmat_src_min";
  for(fsmi=0; fsmi<fsmn; fsmi++) {
    fsm = fsma[fsmi];
    optfsm.value = fsm;
    if(QOP_hisq_force_set_opts(&optfsm, 1)==QOP_FAIL) continue;
    for(nsrc=nsrcmin; nsrc<=nsrcmax; nsrc++) {
      for(bsi=0; bsi<bsn; bsi++) {
	bs = bsa[bsi];
	QDP_set_block_size(bs);
	mf = bench_force(flh, &coeffs, cm, in, nsrc);
	printf0("FF: fsm%3i nsrc%3i bs%5i sec%7.4f mflops = %g\n", fsm, nsrc, bs, secs, mf);
	if(mf>best_mf) {
	  best_mf = mf;
	  best_fsm = fsm;
	  best_nsrc = nsrc;
	  best_bs = bs;
	}
      }
    }
  }

  optfsm.value = best_fsm;
  QOP_hisq_force_set_opts(&optfsm, 1);
  QDP_set_block_size(best_bs);
  QDP_profcontrol(1);
  mf = bench_force(flh, &coeffs, cm, in, best_nsrc);
  QDP_profcontrol(0);
  printf0("prof: FF: fsm%3i nsrc%3i bs%5i sec%7.4f mflops = %g\n", best_fsm, best_nsrc, best_bs, secs, mf);
  printf0("best: FF: fsm%3i nsrc%3i bs%5i mflops = %g\n", best_fsm, best_nsrc, best_bs, best_mf);

  if(QDP_this_node==0) { printf("begin unload links\n"); fflush(stdout); }
  //QOP_hisq_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  QOP_finalize();
}

void
usage(char *s)
{
  printf("%s [n#] [s#] [S#] [x# [# ...]]\n",s);
  printf("\n");
  printf("n\tnumber of iterations\n");
  printf("s\tseed\n");
  printf("S\tnsrc\n");
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
    case 'c' : check=atoi(&argv[i][1]); break;
    case 'f' : fsmn=atoi(&argv[i][1]); break;
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 's' : seed=atoi(&argv[i][1]); break;
    case 'S' : nsrcmax=atoi(&argv[i][1]); break;
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

  if(nsrcmax<0) nsrcmin = nsrcmax = -nsrcmax;

  if(fsmn<0) fsmn = nsrcmax;
  fsma = (int *) malloc(sizeof(*fsma));
  fsma[0] = fsmn;
  fsmn = 1;

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

