// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_FDN_H
#define _QOP_FDN_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_FN_GaugeField *
QOP_FDN_create_G_from_qdp(QDP_DN_ColorMatrix *links[]);

QOP_FN_GaugeField *
QOP_FDN_create_G_from_G(QOP_DN_GaugeField *qopgf_double);

QOP_FN_FermionLinksWilson *
QOP_FDN_wilson_create_L_from_L(QOP_DN_FermionLinksWilson *flw_double);

QOP_FN_FermionLinksAsqtad *
QOP_FDN_asqtad_create_L_from_L(QOP_DN_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 'N'
#  include <qop_dfn_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_FDN_H */
// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_DFN_H
#define _QOP_DFN_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_DN_GaugeField *
QOP_DFN_create_G_from_qdp(QDP_FN_ColorMatrix *links[]);

QOP_DN_GaugeField *
QOP_DFN_create_G_from_G(QOP_FN_GaugeField *qopgf_double);

QOP_DN_FermionLinksWilson *
QOP_DFN_wilson_create_L_from_L(QOP_FN_FermionLinksWilson *flw_double);

QOP_DN_FermionLinksAsqtad *
QOP_DFN_asqtad_create_L_from_L(QOP_FN_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 'N'
#  include <qop_dfn_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_DFN_H */
