/* The Fermilab relative residual */
#include <math.h>

static vector *g_r;
static double g_rel;

static void
relnorm_func(vector *x, int i)
{
  QLA_Real xnrm2;
  r_eq_norm2_v(&xnrm2, x);
  if(xnrm2 == 0.) {
    g_rel += 1.;
  } else {
    QLA_Real rnrm2;
    r_eq_norm2_v(&rnrm2, &g_r[i]);
    g_rel += rnrm2/xnrm2;
  }
}

QLA_Real
QOPPCV(relnorm2)(Vector2 **rsd, Vector2 **out, QDP_Subset subset, int nv)
{
  QLA_Real rel;

  g_rel = 0;
  for(int i=0; i<nv; i++) {
    g_r = expose_V(rsd[i]);
    V_eq_funci(out[i], relnorm_func, subset);
    reset_V(rsd[i]);
  }
  QMP_sum_double(&g_rel);
  rel = g_rel/QDP_subset_len(subset);
  return rel;
}

#if 0
void
QOPPCV(relnorm2)(QLA_Real *rel, Vector *rsd, Vector *out, 
		 QDP_Subset subset)
{
  QDP_Real *relm = QDP_create_R();
  QLA_Real *relx;
  QDP_Int *sub = QDP_create_I();
  vector *rsdx, *outx;
  QLA_Real rsdm, outm;
  QLA_Int one = 1;
  QLA_Int *subx;
  int i;

  /* Make subset mask */
  QDP_I_eq_zero(sub, QDP_all);
  QDP_I_eq_i(sub, &one, subset);
  subx = QDP_expose_I(sub);

  /* Compute relative norm over subset */
  rsdx = expose_V(rsd);
  outx = expose_V(out);
  relx = QDP_expose_R(relm);
  for(i = 0; i < QDP_sites_on_node; i++)
    if(subx[i]){
      r_eq_norm2_v(&rsdm, rsdx+i);
      r_eq_norm2_v(&outm, outx+i);
      relx[i] = (outm==0.) ? 1.0 : (rsdm/outm);
    }
  QDP_reset_R(relm);
  QDP_r_eq_sum_R(rel, relm, subset);
  *rel = *rel/QDP_subset_len(subset);

  QDP_destroy_I(sub);
  QDP_destroy_R(relm);
  reset_V(rsd);
  reset_V(out);
}

#endif
