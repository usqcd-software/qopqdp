/* Some data manipulations for FermionLinksWilson */

/* We need definitions for both precisions */

#include <qop.h>
#include <qop_qdp.h>
#include <qmp.h>

/* We have to redefine these here, since the present design of the
   qop_internal.h header does not permit defining both single and
   double precision types in the same source file. */

#define QOP_malloc(var, type, num)					\
  (var) = (type *) malloc(num*sizeof(type));				\
  if(!(var)) {								\
    QMP_error("Error: QOP ran out of memory in function %s\n", __func__); \
    exit(1);								\
  }

struct QOP_D3_GaugeField_struct {
  QDP_D3_ColorMatrix **links;
  QLA_D3_ColorMatrix **raw;
};


struct QOP_F3_GaugeField_struct {
  QDP_F3_ColorMatrix **links;
  QLA_F3_ColorMatrix **raw;
};


typedef struct {
  QDP_D3_DiracFermion **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_D3_eigcg_t_D;

struct QOP_D3_FermionLinksWilson_struct {
  double clovinvkappa;
  int dblstored;
  QDP_D3_ColorMatrix **links;
  QDP_D3_ColorMatrix **bcklinks;
  QDP_D3_ColorMatrix **dbllinks;
  QOP_D3_GaugeField *qopgf;
  QDP_D3_DiracPropagator *qdpclov;
  double *clov, *clovinv;
  double **rawlinks, *rawclov;
  QOP_D3_eigcg_t_D eigcg;
};

typedef struct {
  QDP_F3_DiracFermion **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_F3_eigcg_t_D;

struct QOP_F3_FermionLinksWilson_struct {
  float clovinvkappa;
  int dblstored;
  QDP_F3_ColorMatrix **links;
  QDP_F3_ColorMatrix **bcklinks;
  QDP_F3_ColorMatrix **dbllinks;
  QOP_F3_GaugeField *qopgf;
  QDP_F3_DiracPropagator *qdpclov;
  float *clov, *clovinv;
  float **rawlinks, *rawclov;
  QOP_F3_eigcg_t_D eigcg;
};

typedef struct {
  int inited;
  int verbosity;
  int proflevel;
  int we_inited_qdp;
  int ndim;
  QDP_Shift neighbor3[4];
  QDP_ShiftDir shiftfwd[8], shiftbck[8];
} QOP_common_t;
extern QOP_common_t QOP_common;

#define QOP_printf0 if(QDP_this_node==0) printf

#define CLOV_REALS (2*6*6) // 2 packed 6x6 Hermitian matrices

/* Create a single-precision copy of a double-precision gauge field */

QOP_F3_GaugeField *
QOP_FD3_create_G_from_G(QOP_D3_GaugeField *qopgf_double){
  QOP_F3_GaugeField *qopgf_single;
  int i,x;

  QOP_malloc(qopgf_single, QOP_F3_GaugeField, 1);

  /* Create and copy the gauge links */
  if(qopgf_double->raw == NULL){
    qopgf_single->raw = NULL;
    QOP_malloc(qopgf_single->links, QDP_F3_ColorMatrix *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      qopgf_single->links[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(qopgf_single->links[i], qopgf_double->links[i], QDP_all);
    }
  } else {
    qopgf_single->links = NULL;
    QOP_printf0("QOP_FD3_create_G_from_G: Warning: raw member is not supported\n");
    QOP_malloc(qopgf_single->raw, QLA_F3_ColorMatrix *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(qopgf_single->raw[i], QLA_F3_ColorMatrix, QDP_sites_on_node);
      for(x=0; x<QDP_sites_on_node; x++)
	QLA_FD3_M_eq_M(qopgf_single->raw[i]+x, qopgf_double->raw[i]+x);
    }
  }
  return qopgf_single;
}

/* Create a single-precision copy of double-precision Wilson fermion links */

QOP_F3_FermionLinksWilson *
QOP_FD3_wilson_create_L_from_L(QOP_D3_FermionLinksWilson *flw_double){
  QOP_F3_FermionLinksWilson *flw_single;
  int i;

  /* Create the parent struct */
  QOP_malloc(flw_single, QOP_F3_FermionLinksWilson, 1);

  /* Copy scalar values */
  flw_single->dblstored = flw_double->dblstored;
  flw_single->clovinvkappa = flw_double->clovinvkappa;

  /* Create and copy the gauge field. */
  /* Here we keep two copies of the pointers to the QDP color matrix
     fields, one in the QOP gauge field member and one in the links
     member. */
  if(flw_double->qopgf != NULL){
    flw_single->links = NULL;
    flw_single->qopgf = QOP_FD3_create_G_from_G(flw_double->qopgf);
    QOP_malloc(flw_single->links, QDP_F3_ColorMatrix *, 4);
    for(i=0; i<4; i++){
      flw_single->links[i] = flw_single->qopgf->links[i];
    }
  } else {
    QOP_printf0("QOP_FD3_wilson_create_L_from_L: Error: missing the gauge field\n");
  }

  /* Create and copy backward links */
  QOP_malloc(flw_single->bcklinks, QDP_F3_ColorMatrix *, 4);
  for(i=0; i<4; i++){
    if(flw_double->dblstored != 0 && flw_double->bcklinks[i] != NULL){
      flw_single->bcklinks[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(flw_single->bcklinks[i], 
		     flw_double->bcklinks[i], QDP_all);
    } else {
      flw_double->bcklinks[i] = NULL;
    }
  }
  
  /* Create and copy double links */
  QOP_malloc(flw_single->dbllinks, QDP_F3_ColorMatrix *, 8);
  for(i=0; i<4; i++) {
    if(flw_double->dblstored != 0){
      flw_single->dbllinks[2*i] = flw_single->links[i];
      flw_single->dbllinks[2*i+1] = flw_single->bcklinks[i];
    } else {
      flw_single->dbllinks[2*i] = NULL;
      flw_single->dbllinks[2*i+1] = NULL;
    }
  }

  /* Create and copy clover term */
  if(flw_double->clov != NULL){
    int size = QDP_sites_on_node*CLOV_REALS;
    int k;
    QOP_malloc(flw_single->clov, float, size);
    for(k=0; k<size; k++){
      flw_single->clov[k] = flw_double->clov[k];
    }
  } else {
    flw_single->clov = NULL;
  }

  /* Create and copy clovinv term */
  if(flw_double->clovinv != NULL){
    int size = QDP_sites_on_node*CLOV_REALS;
    int k;
    QOP_malloc(flw_single->clovinv, float, size);
    for(k=0; k<size; k++){
      flw_single->clovinv[k] = flw_double->clovinv[k];
    }
  } else {
    flw_single->clovinv = NULL;
  }

  /* Do we really need to copy these as well? */
  flw_single->rawlinks = NULL;
  flw_single->rawclov = NULL;
  flw_single->qdpclov = NULL;
  flw_single->eigcg.u = NULL;

  return flw_single;
}

