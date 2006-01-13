#ifndef _QOP_QDP_H
#define _QOP_QDP_H

#include <qop.h>
#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif


  /**********************/
  /*  General routines  */
  /**********************/

  /* single precision */

QOP_F3_ColorVector  *QOP_F3_create_V_from_qdp(QDP_F3_ColorVector *src);
QOP_F3_DiracFermion *QOP_F3_create_D_from_qdp(QDP_F3_DiracFermion *src);
QOP_F3_GaugeField   *QOP_F3_create_G_from_qdp(QDP_F3_ColorMatrix *src[]);
QOP_F3_Force        *QOP_F3_create_F_from_qdp(QDP_F3_ColorMatrix *src[]);

void QOP_F3_extract_V_to_qdp(QDP_F3_ColorVector *d, QOP_F3_ColorVector *src);
void QOP_F3_extract_D_to_qdp(QDP_F3_DiracFermion *d, QOP_F3_DiracFermion *src);
void QOP_F3_extract_G_to_qdp(QDP_F3_ColorMatrix *d[], QOP_F3_GaugeField *src);
void QOP_F3_extract_F_to_qdp(QDP_F3_ColorMatrix *d[], QOP_F3_Force *src);

QOP_F3_ColorVector  *QOP_F3_convert_V_from_qdp(QDP_F3_ColorVector *src);
QOP_F3_DiracFermion *QOP_F3_convert_D_from_qdp(QDP_F3_DiracFermion *src);
QOP_F3_GaugeField   *QOP_F3_convert_G_from_qdp(QDP_F3_ColorMatrix *src[]);
QOP_F3_Force        *QOP_F3_convert_F_from_qdp(QDP_F3_ColorMatrix *src[]);

QDP_F3_ColorVector   *QOP_F3_convert_V_to_qdp(QOP_F3_ColorVector *src);
QDP_F3_DiracFermion  *QOP_F3_convert_D_to_qdp(QOP_F3_DiracFermion *src);
QDP_F3_ColorMatrix  **QOP_F3_convert_G_to_qdp(QOP_F3_GaugeField *src);
QDP_F3_ColorMatrix  **QOP_F3_convert_F_to_qdp(QOP_F3_Force *src);

  /* double precision */

QOP_D3_ColorVector  *QOP_D3_create_V_from_qdp(QDP_D3_ColorVector *src);
QOP_D3_DiracFermion *QOP_D3_create_D_from_qdp(QDP_D3_DiracFermion *src);
QOP_D3_GaugeField   *QOP_D3_create_G_from_qdp(QDP_D3_ColorMatrix *src[]);
QOP_D3_Force        *QOP_D3_create_F_from_qdp(QDP_D3_ColorMatrix *src[]);

void QOP_D3_extract_V_to_qdp(QDP_D3_ColorVector *d, QOP_D3_ColorVector *src);
void QOP_D3_extract_D_to_qdp(QDP_D3_DiracFermion *d, QOP_D3_DiracFermion *src);
void QOP_D3_extract_G_to_qdp(QDP_D3_ColorMatrix *d[], QOP_D3_GaugeField *src);
void QOP_D3_extract_F_to_qdp(QDP_D3_ColorMatrix *d[], QOP_D3_Force *src);

QOP_D3_ColorVector  *QOP_D3_convert_V_from_qdp(QDP_D3_ColorVector *src);
QOP_D3_DiracFermion *QOP_D3_convert_D_from_qdp(QDP_D3_DiracFermion *src);
QOP_D3_GaugeField   *QOP_D3_convert_G_from_qdp(QDP_D3_ColorMatrix *src[]);
QOP_D3_Force        *QOP_D3_convert_F_from_qdp(QDP_D3_ColorMatrix *src[]);

QDP_D3_ColorVector   *QOP_D3_convert_V_to_qdp(QOP_D3_ColorVector *src);
QDP_D3_DiracFermion  *QOP_D3_convert_D_to_qdp(QOP_D3_DiracFermion *src);
QDP_D3_ColorMatrix  **QOP_D3_convert_G_to_qdp(QOP_D3_GaugeField *src);
QDP_D3_ColorMatrix  **QOP_D3_convert_F_to_qdp(QOP_D3_Force *src);


  /*********************/
  /*  Asqtad routines  */
  /*********************/

  /* single precision */

QOP_F3_FermionLinksAsqtad *
  QOP_F3_asqtad_create_L_from_qdp(QDP_F3_ColorMatrix *fatlinks[],
				  QDP_F3_ColorMatrix *longlinks[]);
void QOP_F3_asqtad_extract_L_to_qdp(QDP_F3_ColorMatrix *fatlinks[],
				    QDP_F3_ColorMatrix *longlinks[],
				    QOP_F3_FermionLinksAsqtad *src);
QOP_F3_FermionLinksAsqtad *
  QOP_F3_asqtad_convert_L_from_qdp(QDP_F3_ColorMatrix *fatlinks[],
				   QDP_F3_ColorMatrix *longlinks[]);
void QOP_F3_asqtad_convert_L_to_qdp(QDP_F3_ColorMatrix ***fatlinks,
				    QDP_F3_ColorMatrix ***longlinks,
				    QOP_F3_FermionLinksAsqtad *src);

  /* double precision */

QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_create_L_from_qdp(QDP_D3_ColorMatrix *fatlinks[],
				  QDP_D3_ColorMatrix *longlinks[]);
void QOP_D3_asqtad_extract_L_to_qdp(QDP_D3_ColorMatrix *fatlinks[],
				    QDP_D3_ColorMatrix *longlinks[],
				    QOP_D3_FermionLinksAsqtad *src);
QOP_D3_FermionLinksAsqtad *
  QOP_D3_asqtad_convert_L_from_qdp(QDP_D3_ColorMatrix *fatlinks[],
				   QDP_D3_ColorMatrix *longlinks[]);
void QOP_D3_asqtad_convert_L_to_qdp(QDP_D3_ColorMatrix ***fatlinks,
				    QDP_D3_ColorMatrix ***longlinks,
				    QOP_D3_FermionLinksAsqtad *src);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* single precision */

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_qdp(QDP_F3_ColorMatrix *links[]);
void QOP_F3_wilson_extract_L_to_qdp(QDP_F3_ColorMatrix *links[],
				    QOP_F3_FermionLinksWilson *src);
QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_qdp(QDP_F3_ColorMatrix *links[]);
QDP_F3_ColorMatrix **
  QOP_F3_wilson_convert_L_to_qdp(QOP_F3_FermionLinksWilson *src);

  /* double precision */

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_qdp(QDP_D3_ColorMatrix *links[]);
void QOP_D3_wilson_extract_L_to_qdp(QDP_D3_ColorMatrix *links[],
				    QOP_D3_FermionLinksWilson *src);
QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_qdp(QDP_D3_ColorMatrix *links[]);
QDP_D3_ColorMatrix **
  QOP_D3_wilson_convert_L_to_qdp(QOP_D3_FermionLinksWilson *src);


  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QOP_Precision == 1

#define QOP_create_V_from_qdp QOP_F3_create_V_from_qdp
#define QOP_create_D_from_qdp QOP_F3_create_D_from_qdp
#define QOP_create_G_from_qdp QOP_F3_create_G_from_qdp
#define QOP_create_F_from_qdp QOP_F3_create_F_from_qdp

#define QOP_extract_V_to_qdp QOP_F3_extract_V_to_qdp
#define QOP_extract_D_to_qdp QOP_F3_extract_D_to_qdp
#define QOP_extract_G_to_qdp QOP_F3_extract_G_to_qdp
#define QOP_extract_F_to_qdp QOP_F3_extract_F_to_qdp

#define QOP_convert_V_from_qdp QOP_F3_convert_V_from_qdp
#define QOP_convert_D_from_qdp QOP_F3_convert_D_from_qdp
#define QOP_convert_G_from_qdp QOP_F3_convert_G_from_qdp
#define QOP_convert_F_from_qdp QOP_F3_convert_F_from_qdp

#define QOP_convert_V_to_qdp QOP_F3_convert_V_to_qdp
#define QOP_convert_D_to_qdp QOP_F3_convert_D_to_qdp
#define QOP_convert_G_to_qdp QOP_F3_convert_G_to_qdp
#define QOP_convert_F_to_qdp QOP_F3_convert_F_to_qdp

#define QOP_asqtad_create_L_from_qdp  QOP_F3_asqtad_create_L_from_qdp
#define QOP_asqtad_extract_L_to_qdp   QOP_F3_asqtad_extract_L_to_qdp
#define QOP_asqtad_convert_L_from_qdp QOP_F3_asqtad_convert_L_from_qdp
#define QOP_asqtad_convert_L_to_qdp   QOP_F3_asqtad_convert_L_to_qdp

#define QOP_wilson_create_L_from_qdp  QOP_F3_wilson_create_L_from_qdp
#define QOP_wilson_extract_L_to_qdp   QOP_F3_wilson_extract_L_to_qdp
#define QOP_wilson_convert_L_from_qdp QOP_F3_wilson_convert_L_from_qdp
#define QOP_wilson_convert_L_to_qdp   QOP_F3_wilson_convert_L_to_qdp

#else

#define QOP_create_V_from_qdp QOP_D3_create_V_from_qdp
#define QOP_create_D_from_qdp QOP_D3_create_D_from_qdp
#define QOP_create_G_from_qdp QOP_D3_create_G_from_qdp
#define QOP_create_F_from_qdp QOP_D3_create_F_from_qdp

#define QOP_extract_V_to_qdp QOP_D3_extract_V_to_qdp
#define QOP_extract_D_to_qdp QOP_D3_extract_D_to_qdp
#define QOP_extract_G_to_qdp QOP_D3_extract_G_to_qdp
#define QOP_extract_F_to_qdp QOP_D3_extract_F_to_qdp

#define QOP_convert_V_from_qdp QOP_D3_convert_V_from_qdp
#define QOP_convert_D_from_qdp QOP_D3_convert_D_from_qdp
#define QOP_convert_G_from_qdp QOP_D3_convert_G_from_qdp
#define QOP_convert_F_from_qdp QOP_D3_convert_F_from_qdp

#define QOP_convert_V_to_qdp QOP_D3_convert_V_to_qdp
#define QOP_convert_D_to_qdp QOP_D3_convert_D_to_qdp
#define QOP_convert_G_to_qdp QOP_D3_convert_G_to_qdp
#define QOP_convert_F_to_qdp QOP_D3_convert_F_to_qdp

#define QOP_asqtad_create_L_from_qdp  QOP_D3_asqtad_create_L_from_qdp
#define QOP_asqtad_extract_L_to_qdp   QOP_D3_asqtad_extract_L_to_qdp
#define QOP_asqtad_convert_L_from_qdp QOP_D3_asqtad_convert_L_from_qdp
#define QOP_asqtad_convert_L_to_qdp   QOP_D3_asqtad_convert_L_to_qdp

#define QOP_wilson_create_L_from_qdp  QOP_D3_wilson_create_L_from_qdp
#define QOP_wilson_extract_L_to_qdp   QOP_D3_wilson_extract_L_to_qdp
#define QOP_wilson_convert_L_from_qdp QOP_D3_wilson_convert_L_from_qdp
#define QOP_wilson_convert_L_to_qdp   QOP_D3_wilson_convert_L_to_qdp

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_QDP_H */
