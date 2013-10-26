#ifndef _TEST_COMMON_H
#define _TEST_COMMON_H

#include <qop_qdp.h>
#include <qmp.h>

#if QOP_Colors == 'N'

#define NCPROT int nc,
#define NCPROT1 int nc,
#define NCARG QLA_Nc,
#define NCARG1 QLA_Nc,
#ifndef QLA_ColorMatrix
#if QOP_Precision == 'F'
#define QLA_ColorMatrix(x) QLA_FN_ColorMatrix(QLA_Nc,(x))
#else
#define QLA_ColorMatrix(x) QLA_DN_ColorMatrix(QLA_Nc,(x))
#endif
#define QLA_F_ColorMatrix(x) QLA_FN_ColorMatrix(QLA_Nc,(x))
#define QLA_D_ColorMatrix(x) QLA_DN_ColorMatrix(QLA_Nc,(x))
#endif

#else

#define NCPROT
#define NCPROT1
#define NCARG
#define NCARG1
#ifndef QLA_ColorMatrix
#define QLA_ColorMatrix(x) QLA_ColorMatrix x
#define QLA_F_ColorMatrix(x) QLA_F_ColorMatrix x
#define QLA_D_ColorMatrix(x) QLA_D_ColorMatrix x
#endif

#if QOP_Colors == 1
#undef NCPROT1
#define NCPROT1 int nc,
#undef NCARG1
#define NCARG1 QLA_Nc,
#endif

#endif

#define printf0 if(QDP_this_node==0) printf

extern QDP_RandomState *rs;

void print_layout(void);
void seed_rand(QDP_RandomState *rs, int seed);
void make_unitary(QDP_ColorMatrix **m, int n);
QLA_Real get_plaq(QDP_ColorMatrix *link[]);
void get_random_links(QDP_ColorMatrix **u, int n, QLA_Real r);
void get_latsize(int *ndim, int **ls, char *fn);
void load_lattice(QDP_ColorMatrix *gauge[], char *fn);
void point_source_V(QDP_ColorVector *v, int *x, int c);

#endif /* _TEST_COMMON_H */
