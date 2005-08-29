#include <ctype.h>
#include <time.h>
#include <math.h>
#include <qdp.h>
#include <qop.h>

typedef struct param_struct {
  QLA_Real beta;
  int ndim;
  int nc;
  int *lattice_size;
  int seed;
} param_t;

static param_t params;

typedef struct field_struct {
  QDP_RandomState *rs;
} field_t;

static field_t fields;

static int nit;
static char *ifn, *ofn;

void
lex_int(QLA_Int *li, int coords[])
{
  int i,t;

  t = coords[0];
  for(i=1; i<params.ndim; i++) {
    t = t*QDP_coord_size(i) + coords[i];
  }
  *li = t;
}

void
seed_rand(QDP_RandomState *rs, int seed)
{
  QDP_Int *li;

  li = QDP_create_I();

  QDP_I_eq_func(li, lex_int, QDP_all);
  QDP_S_eq_seed_i_I(rs, seed, li, QDP_all);

  QDP_destroy_I(li);
}

QLA_Complex
det(QLA_ColorMatrix *m)
{
  QLA_ColorMatrix tm;
  QLA_Complex z1;
  int i, j, c;

  for(j = 0; j < QDP_Nc; j++)
    for(i = 0; i < QDP_Nc; i++)
      QLA_elem_M(tm,j,i) = QLA_elem_M(*m,j,i);

  for(j = 0; j < QDP_Nc; j++) {
    for(i = 0; i <= j; i++) {
      QLA_Complex t2;
      t2 = QLA_elem_M(tm,j,i);
      for(c = 0; c < i; c++)
        QLA_c_meq_c_times_c(t2, QLA_elem_M(tm,c,i), QLA_elem_M(tm,j,c));

      QLA_elem_M(tm,j,i) = t2;
    }

    for(i = (j+1); i < QDP_Nc; i++) {
      QLA_Complex t2;
      t2 = QLA_elem_M(tm,j,i);
      for(c = 0; c < j; c++)
        QLA_c_meq_c_times_c(t2, QLA_elem_M(tm,c,i), QLA_elem_M(tm,j,c));

      QLA_c_eq_c_div_c(QLA_elem_M(tm,j,i), t2, QLA_elem_M(tm,j,j));
    }
  }

  /* The determinant */
  z1 = QLA_elem_M(tm,0,0);
  for(c = 1; c < QDP_Nc; c++) {
    QLA_Complex z;
    QLA_c_eq_c_times_c(z, z1, QLA_elem_M(tm,c,c));
    z1 = z;
  }

  return z1;
}

void
normalize(QLA_ColorMatrix *m, int r)
{
  QLA_Real n, t;
  int c;

  n = 0;
  for(c=0; c<QDP_Nc; c++) {
    QLA_R_eq_norm2_C(&t, &QLA_elem_M(*m, r, c));
    n += t;
  }
  n = 1/sqrt(n);
  for(c=0; c<QDP_Nc; c++) {
    QLA_C_eq_r_times_C(&QLA_elem_M(*m,r,c), &n, &QLA_elem_M(*m,r,c));
  }
}

void
orthogonalize(QLA_ColorMatrix *m, int r1, int r2)
{
  QLA_Complex z, t;
  int c;

  QLA_C_eq_zero(&z);
  for(c=0; c<QDP_Nc; c++) {
    QLA_C_eq_C_dot_C(&t, &QLA_elem_M(*m,r1,c), &QLA_elem_M(*m,r2,c));
    QLA_C_peq_C(&z, &t);
  }
  for(c=0; c<QDP_Nc; c++) {
    QLA_C_meq_C_times_C(&QLA_elem_M(*m,r2,c), &z, &QLA_elem_M(*m,r1,c));
  }
}

void
make_unitary_func(QLA_ColorMatrix *m, int coords[])
{
  QLA_Complex z1, z2;
  QLA_Real r;
  int i, j, c;

  for(i=0; i<QDP_Nc; i++) {
    for(j=0; j<i; j++) {
      orthogonalize(m, j, i);
    }
    normalize(m, i);
  }

  z1 = det(m);

  r = QLA_norm_c(z1);
  QLA_c_eq_ca(z2, z1);
  QLA_c_eq_c_div_r(z1, z2, r);
  for(c = 0; c < QDP_Nc; ++c) {
    QLA_Complex z;
    QLA_c_eq_c_times_c(z, QLA_elem_M(*m,c,QDP_Nc-1), z1);
    QLA_elem_M(*m,c,QDP_Nc-1) = z;
  }
  //z1 = det(m);
  //printf("%g\t%g\n", QLA_real(z1), QLA_imag(z1));
}

void
make_unitary(QDP_ColorMatrix **m, int n)
{
  int i;

  for(i=0; i<n; i++) {
    QDP_M_eq_func(m[i], make_unitary_func, QDP_all);
  }
}

QLA_Real
get_plaq(QDP_ColorMatrix *link[])
{
  int mu, nu;
  QLA_Real plaq;
  QDP_ColorMatrix *temp1, *temp2, *temp3, *temp4;

#ifdef LOCAL_SUM
  QDP_Real *treal1, *treal2;
  treal1 = QDP_create_R();
  treal2 = QDP_create_R();
  QDP_R_eq_zero(treal2, QDP_all);
#else
  QLA_Real tplaq;
  plaq = 0;
#endif

  temp1 = QDP_create_M();
  temp2 = QDP_create_M();
  temp3 = QDP_create_M();
  temp4 = QDP_create_M();

  for(mu=0; mu<params.ndim-1; ++mu) {
    for(nu=mu+1; nu<params.ndim; ++nu) {

      QDP_M_eq_sM(temp1, link[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_sM(temp2, link[mu], QDP_neighbor[nu], QDP_forward, QDP_all);

      QDP_M_eq_Ma_times_M(temp3, link[nu], link[mu], QDP_all);

      QDP_M_eq_M_times_M(temp4, temp3, temp1, QDP_all);
      QDP_discard_M(temp1);

#ifdef LOCAL_SUM
      QDP_R_eq_re_M_dot_M(treal1, temp2, temp4, QDP_all);
      QDP_discard_M(temp2);
      QDP_R_peq_R(treal2, treal1, QDP_all);
#else
      QDP_r_eq_re_M_dot_M(&tplaq, temp2, temp4, QDP_all);
      QDP_discard_M(temp2);
      plaq += tplaq;
#endif

    }
  }

#ifdef LOCAL_SUM
  QDP_r_eq_sum_R(&plaq, treal2, QDP_all);
  QDP_destroy_R(treal1);
  QDP_destroy_R(treal2);
#endif

  QDP_destroy_M(temp1);
  QDP_destroy_M(temp2);
  QDP_destroy_M(temp3);
  QDP_destroy_M(temp4);

  return plaq/(0.5*params.ndim*(params.ndim-1)*QDP_volume());
}

void
start(void)
{
  int i;
  QDP_ColorMatrix **u;
  QDP_RandomState *rs;
  QLA_Real plaq;

  rs = QDP_create_S();
  seed_rand(rs, params.seed);
  fields.rs = rs;

  u = (QDP_ColorMatrix **) malloc(params.ndim*sizeof(QDP_ColorMatrix *));
  for(i=0; i<params.ndim; i++) {
    u[i] = QDP_create_M();
    if(ifn==NULL) {
      if(1) {
	QDP_M_eq_gaussian_S(u[i], rs, QDP_all);
      } else {
	QLA_Complex z;
	QLA_real(z) = 1;
	QLA_imag(z) = 0;
	QDP_M_eq_c(u[i], &z, QDP_all);
      }
    }
  }
  if(ifn!=NULL) {
    QDP_String *md;
    QDP_Reader *qr;
    md = QDP_string_create();
    qr = QDP_open_read(md, ifn);
    //QDP_FN_vread_M(params.nc, qr, md, u, params.ndim);
    QDP_F3_vread_M(qr, md, u, params.ndim);
    QDP_close_read(qr);
    QDP_string_destroy(md);
  }
  make_unitary(u, params.ndim);

  //plaq = get_plaq(u)/(0.5*params.ndim*(params.ndim-1));
  plaq = get_plaq(u);
  //printf("%-8i%-14g%-14g\n", 0, plaq, QLA_real(get_ploop(u)));
  if(QDP_this_node==0) printf("plaquette = %g\n", plaq);
  //plaq = get_plaq2(u)/(0.5*params.ndim*(params.ndim-1));
  //plaq = get_plaq2(u);
  //printf("%-8i%-14g\n", 0, plaq);

#if 0
  for(i=1; i<=nit; i++) {
    update(u);
    make_unitary(u, params.ndim);
    //plaq = get_plaq(u)/(0.5*params.ndim*(params.ndim-1));
    plaq = get_plaq(u);
    if(QDP_this_node==0) {
      printf("%-8i%-14g%-14g\n", i, plaq, QLA_real(get_ploop(u)));
    }
  }
#endif

  QOP_layout qoplayout;
  qoplayout.ndims = params.ndim;
  for(i=0; i<params.ndim; i++) {
    qoplayout.sites[i] = params.lattice_size[i];
    qoplayout.bc[i] = 1;
  }

  QDP_ColorMatrix *fatlinks[4], *longlinks[4];
  QDP_ColorVector *out, *in;
  out = QDP_create_V();
  in = QDP_create_V();
  for(i=0; i<4; i++) {
    fatlinks[i] = u[i];
    longlinks[i] = u[i];
  }
  QDP_V_eq_gaussian_S(in, rs, QDP_all);
  QDP_V_eq_zero(out, QDP_all);

  QOP_invert_arg inv_arg;
  inv_arg.mass = 0.1;
  inv_arg.rsqmin = 1e-3;
  inv_arg.max_iter = 500;
  inv_arg.restart = 100;
  inv_arg.evenodd = QOP_EVEN;

  if(QDP_this_node==0) { printf("begin init\n"); fflush(stdout); }
  QOP_asqtad_invert_init(&qoplayout);
  if(QDP_this_node==0) { printf("begin load links\n"); fflush(stdout); }
  QOP_asqtad_invert_load_links_qdp(fatlinks, longlinks);
  if(QDP_this_node==0) { printf("begin invert\n"); fflush(stdout); }
  nit = QOP_F_asqtad_inv_qdp(&inv_arg, out, in);
  if(QDP_this_node==0) { printf("begin unload links\n"); fflush(stdout); }
  QOP_F_asqtad_invert_unload_links();
  if(QDP_this_node==0) { printf("begin finalize\n"); fflush(stdout); }
  QOP_asqtad_invert_finalize();

  if(QDP_this_node==0) printf("number of iterations = %i\n", nit);

  if(ofn!=NULL) {
    QDP_String *md;
    QDP_Writer *qw;
    md = QDP_string_create();
    QDP_string_set(md, "test");
    qw = QDP_open_write(md, ofn, QDP_SINGLEFILE);
    //QDP_FN_vwrite_M(params.nc, qw, md, u, params.ndim);
    QDP_F3_vwrite_M(qw, md, u, params.ndim);
    QDP_close_write(qw);
    QDP_string_destroy(md);
  }

}

void
usage(char *s)
{
  printf("%s [b#] [c#] [d#] [i<fn>] [n#] [o<fn>] [s#] [x# [# ...]]\n",s);
  printf("\n");
  printf("b\tbeta\n");
  printf("c\tnumber of colors\n");
  printf("d\tnumber of dimensions\n");
  printf("i\tinput file\n");
  printf("n\tnumber of iterations\n");
  printf("o\toutput file\n");
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

  params.ndim = 4;
  params.nc = 3;
  params.beta = 6.0;
  nit = 10;
  params.seed = time(NULL);
  ifn = NULL;
  ofn = NULL;
  j = 0;
  for(i=1; i<argc; i++) {
    switch(argv[i][0]) {
    case 'b' : params.beta=atof(&argv[i][1]); break;
    case 'c' : params.nc=atoi(&argv[i][1]); break;
    case 'd' : params.ndim=atoi(&argv[i][1]); break;
    case 'i' : ifn=&argv[i][1]; break;
    case 'n' : nit=atoi(&argv[i][1]); break;
    case 'o' : ofn=&argv[i][1]; break;
    case 's' : params.seed=atoi(&argv[i][1]); break;
    case 'x' : j=i; while((i+1<argc)&&(isdigit(argv[i+1][0]))) ++i; break;
    default : usage(argv[0]);
    }
  }

  params.lattice_size = (int *) malloc(params.ndim*sizeof(int));
  if(j==0) {
    for(i=0; i<params.ndim; ++i) params.lattice_size[i] = 8;
  } else {
    if(!isdigit(argv[j][1])) usage(argv[0]);
    params.lattice_size[0] = atoi(&argv[j][1]);
    for(i=1; i<params.ndim; ++i) {
      if((++j<argc)&&(isdigit(argv[j][0]))) {
        params.lattice_size[i] = atoi(&argv[j][0]);
      } else {
        params.lattice_size[i] = params.lattice_size[i-1];
      }
    }
  }

  if(QDP_this_node==0) {
    printf("ndims = %i\n", params.ndim);
    printf("size = %i", params.lattice_size[0]);
    for(i=1; i<params.ndim; i++) {
      printf(" %i", params.lattice_size[i]);
    }
    printf("\n");
    printf("Nc = %i\n", params.nc);
    printf("seed = %i\n", params.seed);
  }

  QDP_set_latsize(params.ndim, params.lattice_size);
  QDP_create_layout();

  start();

  return 0;
}

