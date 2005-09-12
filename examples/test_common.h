#ifndef _TEST_COMMON_H
#define _TEST_COMMON_H

#include <qdp.h>
#include <qop.h>

#define printf0 if(QDP_this_node==0) printf

extern QDP_RandomState *rs;

void seed_rand(QDP_RandomState *rs, int seed);
void make_unitary(QDP_ColorMatrix **m, int n);
QLA_Real get_plaq(QDP_ColorMatrix *link[]);
void get_random_links(QDP_ColorMatrix **u, int n, QLA_Real r);

#endif /* _TEST_COMMON_H */
