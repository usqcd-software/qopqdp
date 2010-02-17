/* Created in emulation of wilson_dslash */

#include <string.h>
#include <qop.h>
#include <qop_qdp.h>
#include <qmp.h>

#define QOP_malloc(var, type, num)					\
  (var) = (type *) malloc(num*sizeof(type));				\
  if(!(var)) {								\
    QMP_error("Error: QOP ran out of memory in function %s\n", __func__); \
    exit(1);								\
  }

/* We have to redefine these here, since the present design of the
   qop_internal.h header does not permit defining both single and
   double precision types in the same source file. */

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

struct QOP_D3_FermionLinksDW_struct {
  int dblstored;
  QDP_D3_ColorMatrix **links;
  QDP_D3_ColorMatrix **bcklinks;
  QDP_D3_ColorMatrix **dbllinks;
  QOP_D3_GaugeField *qopgf;
  //double **rawlinks;
  //QOP_D3_eigcg_t_D eigcg;
};

typedef struct {
  QDP_F3_DiracFermion **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_F3_eigcg_t_D;

struct QOP_F3_FermionLinksDW_struct {
  int dblstored;
  QDP_F3_ColorMatrix **links;
  QDP_F3_ColorMatrix **bcklinks;
  QDP_F3_ColorMatrix **dbllinks;
  QOP_F3_GaugeField *qopgf;
  //float **rawlinks;
  //QOP_F3_eigcg_t_D eigcg;
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

// Create a single-precision copy of a double-precision gauge field
QOP_F3_GaugeField *
QOP_FD3_create_G_from_G(QOP_D3_GaugeField *qopgf_double){
  QOP_F3_GaugeField *qopgf_single;
  int i,x;

  QOP_malloc(qopgf_single, QOP_F3_GaugeField, 1);

  /* Create and copy the gauge links */
  if (qopgf_double->raw == NULL) {
    qopgf_single->raw = NULL;
    QOP_malloc(qopgf_single->links, QDP_F3_ColorMatrix *, QOP_common.ndim);
    for (i=0; i<QOP_common.ndim; i++) {
      qopgf_single->links[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(qopgf_single->links[i], qopgf_double->links[i],
                     QDP_all);
    }
  } else {
    qopgf_single->links = NULL;
    QOP_printf0("QOP_FD3_create_G_from_G: Warning: raw member is not supported\n");
    QOP_malloc(qopgf_single->raw, QLA_F3_ColorMatrix *, QOP_common.ndim);
    for (i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(qopgf_single->raw[i], QLA_F3_ColorMatrix,
                 QDP_sites_on_node);
      for (x=0; x<QDP_sites_on_node; x++)
	      QLA_FD3_M_eq_M(qopgf_single->raw[i]+x, qopgf_double->raw[i]+x);
    }
  }
  return qopgf_single;
}

/* Create a single-precision copy of double-precision Wilson fermion links */

QOP_F3_FermionLinksDW *
QOP_FD3_dw_create_L_from_L(QOP_D3_FermionLinksDW *fldw_double){
  QOP_F3_FermionLinksDW *fldw_single;
  int i;

  /* Create the parent struct */
  QOP_malloc(fldw_single, QOP_F3_FermionLinksDW, 1);

  /* Copy scalar values */
  fldw_single->dblstored = fldw_double->dblstored;

  /* Create and copy the gauge field. */
  /* Here we keep two copies of the pointers to the QDP color matrix
     fields, one in the QOP gauge field member and one in the links
     member. */
  if (fldw_double->qopgf != NULL){
    fldw_single->links = NULL;
    fldw_single->qopgf = QOP_FD3_create_G_from_G(fldw_double->qopgf);
    QOP_malloc(fldw_single->links, QDP_F3_ColorMatrix *, 4);
    for (i=0; i<4; i++){
      fldw_single->links[i] = fldw_single->qopgf->links[i];
    }
  } else {
    QOP_printf0("QOP_FD3_dw_create_L_from_L: Error: missing the gauge field\n");
  }

  /* Create and copy backward links */
  QOP_malloc(fldw_single->bcklinks, QDP_F3_ColorMatrix *, 4);
  for (i=0; i<4; i++){
    if (fldw_double->dblstored != 0 && fldw_double->bcklinks[i] != NULL){
      fldw_single->bcklinks[i] = QDP_F3_create_M();
      QDP_FD3_M_eq_M(fldw_single->bcklinks[i], 
		     fldw_double->bcklinks[i], QDP_all);
    } else {
      fldw_double->bcklinks[i] = NULL;
    }
  }
  
  /* Create and copy double links */
  QOP_malloc(fldw_single->dbllinks, QDP_F3_ColorMatrix *, 8);
  for (i=0; i<4; i++) {
    if(fldw_double->dblstored != 0){
      fldw_single->dbllinks[2*i] = fldw_single->links[i];
      fldw_single->dbllinks[2*i+1] = fldw_single->bcklinks[i];
    } else {
      fldw_single->dbllinks[2*i] = NULL;
      fldw_single->dbllinks[2*i+1] = NULL;
    }
  }

  return fldw_single;
}

