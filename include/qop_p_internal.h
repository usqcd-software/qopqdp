#ifndef _QOP_P_INTERNAL_H
#define _QOP_P_INTERNAL_H

#ifdef HAVE_NC3
#ifdef HAVE_NCN
#include <qop_mg_internal.h>

// V-cycle

typedef struct QOP_IP_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_P_Gcr *gcr;
  QOP_P_Cgls *cgls;
  QDP_N_ColorVector **r;
  QDP_N_ColorVector **p;
  QDP_N_ColorVector **Ap;
} QOP_IP_MgVcycleArgs;

void QOP_IP_mgVcycle(QDP_N_ColorVector **out, QDP_N_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN
#endif // HAVE_NC3

#if QOP_Precision == _QOP_Precision
#  include <qop_p_internal_generic.h>
#endif

#endif /* _QOP_P_INTERNAL_H */
