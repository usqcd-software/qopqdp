/* The Fermilab relative residual */
//#include <math.h>

static double g_rel;
static QLA_Real *g_r;

#define NC nc
static void
relnorm_func(NCPROT vector(*x), int i)
{
  vector(*r) = (vector(*)) g_r;
  QLA_Real xnrm2;
  r_eq_norm2_v(&xnrm2, x);
  if(xnrm2 == 0.) {
    g_rel += 1.;
  } else {
    QLA_Real rnrm2;
    r_eq_norm2_v(&rnrm2, &r[i]);
    g_rel += rnrm2/xnrm2;
  }
}
#undef NC

QLA_Real
QOPPCV(relnorm2)(Vector2 **rsd, Vector2 **out, QDP_Subset subset, int nv)
{
  g_rel = 0;
  for(int i=0; i<nv; i++) {
    g_r = expose_V(rsd[i]);
    V_eq_funci(out[i], relnorm_func, subset);
    reset_V(rsd[i]);
  }
  QMP_sum_double(&g_rel);
  QLA_Real rel = g_rel/QDP_subset_len(subset);
  return rel;
}
