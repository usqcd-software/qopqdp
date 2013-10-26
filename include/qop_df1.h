// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_FD1_H
#define _QOP_FD1_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_F1_GaugeField *
QOP_FD1_create_G_from_G(QOP_D1_GaugeField *qopgf_double);

QOP_F1_FermionLinksWilson *
QOP_FD1_wilson_create_L_from_L(QOP_D1_FermionLinksWilson *flw_double);

QOP_F1_FermionLinksAsqtad *
QOP_FD1_asqtad_create_L_from_L(QOP_D1_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 1
#  include <qop_df1_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_FD1_H */
// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_DF1_H
#define _QOP_DF1_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_D1_GaugeField *
QOP_DF1_create_G_from_G(QOP_F1_GaugeField *qopgf_double);

QOP_D1_FermionLinksWilson *
QOP_DF1_wilson_create_L_from_L(QOP_F1_FermionLinksWilson *flw_double);

QOP_D1_FermionLinksAsqtad *
QOP_DF1_asqtad_create_L_from_L(QOP_F1_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 1
#  include <qop_df1_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_DF1_H */
