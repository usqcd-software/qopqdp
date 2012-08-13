// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_FD3_H
#define _QOP_FD3_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_F3_GaugeField *
QOP_FD3_create_G_from_G(QOP_D3_GaugeField *qopgf_double);

QOP_F3_FermionLinksWilson *
QOP_FD3_wilson_create_L_from_L(QOP_D3_FermionLinksWilson *flw_double);

QOP_F3_FermionLinksAsqtad *
QOP_FD3_asqtad_create_L_from_L(QOP_D3_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 3
#  include <qop_df3_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_FD3_H */
// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_DF3_H
#define _QOP_DF3_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_D3_GaugeField *
QOP_DF3_create_G_from_G(QOP_F3_GaugeField *qopgf_double);

QOP_D3_FermionLinksWilson *
QOP_DF3_wilson_create_L_from_L(QOP_F3_FermionLinksWilson *flw_double);

QOP_D3_FermionLinksAsqtad *
QOP_DF3_asqtad_create_L_from_L(QOP_F3_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 3
#  include <qop_df3_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_DF3_H */
