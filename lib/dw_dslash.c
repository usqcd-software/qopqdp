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
  QOP_D3_FermionLinksWilson *flw;
};

typedef struct {
  QDP_F3_DiracFermion **u;
  QLA_Real *l;
  int numax, nu, nev, m, nv;
  int nn, addvecs;
} QOP_F3_eigcg_t_D;

struct QOP_F3_FermionLinksDW_struct {
  QOP_F3_FermionLinksWilson *flw;
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

// Create a single-precision copy of double-precision Wilson links
QOP_F3_FermionLinksDW *
QOP_FD3_dw_create_L_from_L( QOP_D3_FermionLinksDW *fldw_double ) {

  // Create the parent struct
  QOP_F3_FermionLinksDW *fldw_single;
  QOP_malloc(fldw_single, QOP_F3_FermionLinksDW, 1);

  // Invoke the Wilson conversion tool
  fldw_single->flw = QOP_FD3_wilson_create_L_from_L(fldw_double->flw);

  return fldw_single;
}

