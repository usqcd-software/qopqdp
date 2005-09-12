#ifndef _QOPQDP_INTERNAL_H
#define _QOPQDP_INTERNAL_H

#include <qop.h>

#define printf0 if(QDP_this_node==0) printf

#ifdef __cplusplus
extern "C" {
#endif

double dclock(void);

#ifdef __cplusplus
}
#endif

#endif /* _QOPQDP_INTERNAL_H */
