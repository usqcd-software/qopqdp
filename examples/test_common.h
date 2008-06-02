#ifndef _TEST_COMMON_H
#define _TEST_COMMON_H

#include <qop_qdp.h>
#include <qmp.h>

#define printf0 if(QDP_this_node==0) printf

extern QDP_RandomState *rs;

void print_layout(void);
void seed_rand(QDP_RandomState *rs, int seed);
void make_unitary(QDP_ColorMatrix **m, int n);
QLA_Real get_plaq(QDP_ColorMatrix *link[]);
void get_random_links(QDP_ColorMatrix **u, int n, QLA_Real r);
void get_latsize(int *ndim, int **ls, char *fn);
void load_lattice(QDP_ColorMatrix *gauge[], char *fn);

#endif /* _TEST_COMMON_H */
