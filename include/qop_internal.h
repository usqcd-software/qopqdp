#ifndef _QOPQDP_INTERNAL_H
#define _QOPQDP_INTERNAL_H

#include <qop.h>

#define printf0 if(QDP_this_node==0) printf

#if QOP_Precision == 2
#define PREC(x) QOP_D_##x
#define REAL double
#else
#define PREC(x) QOP_F_##x
#define REAL float
#endif

#ifdef __cplusplus
extern "C" {
#endif

double dclock(void);

#ifdef __cplusplus
}
#endif

#endif /* _QOPQDP_INTERNAL_H */
