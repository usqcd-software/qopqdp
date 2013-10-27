// DO NOT EDIT
// generated from qop_p_internal.h
#ifndef _QOP_D_INTERNAL_H
#define _QOP_D_INTERNAL_H

#ifdef HAVE_NC3
#ifdef HAVE_NCN
#include <qop_mg_internal.h>

// V-cycle

typedef struct QOP_D_MgVcycleArgs {
  double s;
  double tpre;
  double tcoarse;
  double tpost;
  double delta;
  void (*op)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *opargs;
  void (*cop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *copargs;
  void (*sop)(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);
  void *sopargs;
  int nv;
  int npre;
  int npost;
  int indent;
  int verbose;
  int count;
  QDP_Subset sub;
  QDP_Subset sub2;
  QOP_D_Gcr *gcr;
  QOP_D_Cgls *cgls;
  QDP_DN_ColorVector **r;
  QDP_DN_ColorVector **p;
  QDP_DN_ColorVector **Ap;
} QOP_D_MgVcycleArgs;

void QOP_D_mgVcycle(QDP_DN_ColorVector **out, QDP_DN_ColorVector **in, int sign, void *args);

#endif // HAVE_NCN
#endif // HAVE_NC3

#if QOP_Precision == 'D'
#  include <qop_d_internal_generic.h>
#endif

#endif /* _QOP_D_INTERNAL_H */
