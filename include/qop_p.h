#ifndef _QOP_P_H
#define _QOP_P_H

#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif

QOP_P_test

  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == _QOP_Precision
#  include <qop_p_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_P_H */
