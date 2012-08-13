#ifndef _QOP_POC_H
#define _QOP_POC_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_PC_GaugeField *
QOP_POC_create_G_from_G(QOP_OC_GaugeField *qopgf_double);

QOP_PC_FermionLinksWilson *
QOP_POC_wilson_create_L_from_L(QOP_OC_FermionLinksWilson *flw_double);

QOP_PC_FermionLinksAsqtad *
QOP_POC_asqtad_create_L_from_L(QOP_OC_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == _QOP_Colors
#  include <qop_dfc_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_POC_H */
