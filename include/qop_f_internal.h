// DO NOT EDIT
// generated from qop_p_internal.h
#ifndef _QOP_F_INTERNAL_H
#define _QOP_F_INTERNAL_H

#ifdef HAVE_NC3
#ifdef HAVE_NCN
#include <qop_mg_internal.h>

// V-cycle

typedef struct QOP_F_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_F_Gcr *gcr;
  QOP_F_Cgls *cgls;
  QDP_FN_ColorVector **r;
  QDP_FN_ColorVector **p;
  QDP_FN_ColorVector **Ap;
} QOP_F_MgVcycleArgs;

void QOP_F_mgVcycle(QDP_FN_ColorVector **out, QDP_FN_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN
#endif // HAVE_NC3

#if QOP_Precision == 'F'
#  include <qop_f_internal_generic.h>
#endif

#endif /* _QOP_F_INTERNAL_H */
