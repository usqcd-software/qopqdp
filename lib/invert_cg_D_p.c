#include <qop_internal.h>
#include <generic_D.h>

#ifdef HAVE_BLAS
#include "invert_cg_blas_p.c"
#else
#include "invert_cg_p.c"
#endif
