// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_FD2_H
#define _QOP_FD2_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_F2_GaugeField *
QOP_FD2_create_G_from_G(QOP_D2_GaugeField *qopgf_double);

QOP_F2_FermionLinksWilson *
QOP_FD2_wilson_create_L_from_L(QOP_D2_FermionLinksWilson *flw_double);

QOP_F2_FermionLinksAsqtad *
QOP_FD2_asqtad_create_L_from_L(QOP_D2_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 2
#  include <qop_df2_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_FD2_H */
// DO NOT EDIT
// generated from qop_poc.h
#ifndef _QOP_DF2_H
#define _QOP_DF2_H

#ifdef __cplusplus
extern "C" {
#endif

  /* create QOP object from one of different precision */

QOP_D2_GaugeField *
QOP_DF2_create_G_from_G(QOP_F2_GaugeField *qopgf_double);

QOP_D2_FermionLinksWilson *
QOP_DF2_wilson_create_L_from_L(QOP_F2_FermionLinksWilson *flw_double);

QOP_D2_FermionLinksAsqtad *
QOP_DF2_asqtad_create_L_from_L(QOP_F2_FermionLinksAsqtad *fla_src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Colors == 2
#  include <qop_df2_color_generic.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* _QOP_DF2_H */
