#include <test_common.h>
#include <math.h>

QDP_RandomState *rs;

static void
lex_int(QLA_Int *li, int coords[])
{
  int i,t;

  t = coords[0];
  for(i=1; i<QDP_ndim(); i++) {
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

static QLA_Complex
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

static void
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

static void
orthogonalize(QLA_ColorMatrix *m, int r1, int r2)
{
  QLA_Complex z, t;
  int c;

  QLA_C_eq_zero(&z);
  for(c=0; c<QDP_Nc; c++) {
    QLA_C_eq_C_dot_C(&t, &QLA_elem_M(*m,r1,c), &QLA_elem_M(*m,r2,c));
    QLA_c_peq_c(z, t);
  }
  for(c=0; c<QDP_Nc; c++) {
    QLA_C_meq_C_times_C(&QLA_elem_M(*m,r2,c), &z, &QLA_elem_M(*m,r1,c));
  }
}

static void
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

  for(mu=0; mu<QDP_ndim()-1; ++mu) {
    for(nu=mu+1; nu<QDP_ndim(); ++nu) {

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

  return plaq/(0.5*QDP_ndim()*(QDP_ndim()-1)*QDP_volume());
}

void
get_random_links(QDP_ColorMatrix **u, int n, QLA_Real r)
{
  int i;
  QDP_ColorMatrix *cm = QDP_create_M();

  for(i=0; i<n; i++) {
    QLA_Complex z;
    QLA_c_eq_r(z, 1);
    QDP_M_eq_c(u[i], &z, QDP_all);
    QDP_M_eq_gaussian_S(cm, rs, QDP_all);
    QDP_M_peq_r_times_M(u[i], &r, cm, QDP_all);
  }
  QDP_destroy_M(cm);
  make_unitary(u, n);
}
