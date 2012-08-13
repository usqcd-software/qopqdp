#define QDP_Precision 'D'
#define QDP_Nc 2
#define QOP_Precision 'D'
#define QOP_Nc 2
#include <qop_internal.h>
#include <generic_vD.h>

#ifdef HAVE_BLAS
#include "invert_cg_blas_p.c"
#else
#include "invert_cg_p.c"
#endif
