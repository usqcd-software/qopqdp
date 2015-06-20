/* Some data manipulations for FermionLinksWilson */

/* We need definitions for both precisions */

#include <qop.h>
#include <qop_qdp.h>
#include <qop_internal.h>

#define CLOV_REALS (2*6*6) // 2 packed 6x6 Hermitian matrices
#define CLOV_SIZE (CLOV_REALS*sizeof(REAL)) 

/* Create a single-precision copy of a double-precision gauge field */

QOP_F3_GaugeField *
QOP_FD3_create_G_from_G(QOP_D3_GaugeField *qopgf_double){
  QOP_F3_GaugeField *qopgf_single;
  QDP_Lattice *lat = QDP_D_get_lattice_M(qopgf_double->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  int sites_on_node = QDP_sites_on_node_L(lat);
  int i,x;

  QOP_malloc(qopgf_single, QOP_F3_GaugeField, 1);

  /* Create and copy the gauge links */
  if(qopgf_double->raw == NULL){
    qopgf_single->raw = NULL;
    QOP_malloc(qopgf_single->links, QDP_F3_ColorMatrix *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      qopgf_single->links[i] = QDP_F3_create_M_L(lat);
      QDP_FD3_M_eq_M(qopgf_single->links[i], qopgf_double->links[i], all);
    }
  } else {
    qopgf_single->links = NULL;
    QOP_printf0("QOP_FD3_create_G_from_G: Warning: raw member is not supported\n");
    QOP_malloc(qopgf_single->raw, QLA_F_Real *, QOP_common.ndim);
    for(i=0; i<QOP_common.ndim; i++) {
      QOP_malloc(qopgf_single->raw[i], QLA_F_Real, 18*sites_on_node);
      for(x=0; x<18*sites_on_node; x++)
	  QLA_FD_R_eq_R(qopgf_single->raw[i]+x, qopgf_double->raw[i]+x);
    }
  }
  return qopgf_single;
}

#if 0

/* Create a single-precision copy of double-precision Wilson fermion links */

QOP_F3_FermionLinksWilson *
QOP_FD3_wilson_create_L_from_L(QOP_D3_FermionLinksWilson *flw_double){

  QOP_F3_FermionLinksWilson *flw_single;
  QDP_Lattice *lat = QDP_D_get_lattice_M(qopgf_double->links[0]);
  QDP_Subset all = QDP_all_L(lat);
  int sites_on_node = QDP_sites_on_node_L(lat);
  int i;

  /* Create the parent struct */
  flw_single = QOP_wilson_initialize_gauge_L();

  /* Copy scalar values */
  flw_single->dblstored = flw_double->dblstored;
  flw_single->clovinvkappa = flw_double->clovinvkappa;

  /* Create and copy the modified gauge links */

  if(flw_double->links != NULL) {
    QOP_malloc(flw_single->links, QDP_F3_ColorMatrix *, 4);
    for(i=0; i<4; i++){
      flw_single->links[i] = QDP_F3_create_M_L(lat);
      QDP_FD3_M_eq_M(flw_single->links[i], flw_double->links[i], all);
    }
  } else {
    QOP_printf0("QOP_FD3_wilson_create_L_from_L: Error: missing gauge links\n");
  }

  /* Create and copy the modified backward gauge links */
  QOP_malloc(flw_single->bcklinks, QDP_F3_ColorMatrix *, 4);
  for(i=0; i<4; i++){
    if(flw_double->dblstored != 0 && flw_double->bcklinks[i] != NULL){
      flw_single->bcklinks[i] = QDP_F3_create_M_L(lat);
      QDP_FD3_M_eq_M(flw_single->bcklinks[i], 
		     flw_double->bcklinks[i], all);
    } else {
      flw_double->bcklinks[i] = NULL;
    }
  }
  
  /* Create and copy the double links */
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

  /* Create and copy the clover term */
  if(flw_double->clov != NULL){
    int size = sites_on_node*CLOV_REALS;
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
    int size = sites_on_node*CLOV_REALS;
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

#endif
