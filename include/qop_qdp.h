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

void QOP_F3_asqtad_load_L_from_qdp(QOP_F3_FermionLinksAsqtad *asqtad,
				   QDP_F3_ColorMatrix *fatlinks[],
				   QDP_F3_ColorMatrix *longlinks[]);

void QOP_F3_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_F3_FermionLinksAsqtad *asqtad,
			      float mass,
			      QDP_F3_ColorVector *out,
			      QDP_F3_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_F3_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_F3_FermionLinksAsqtad *asqtad,
			       float mass,
			       QDP_F3_ColorVector *out,
			       QDP_F3_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_F3_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_F3_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      float mass,
			      QDP_F3_ColorVector *out,
			      QDP_F3_ColorVector *in);

void QOP_F3_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_F3_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    float *masses[],
				    int nmass[],
				    QDP_F3_ColorVector **out[],
				    QDP_F3_ColorVector *in[],
				    int nsrc);

void QOP_F3_asqtad_get_eigcg(QOP_F3_FermionLinksAsqtad *asqtad,
			     QDP_F3_ColorVector **evecs,
			     QLA_F_Real *evals, int *nv);

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

void QOP_D3_asqtad_load_L_from_qdp(QOP_D3_FermionLinksAsqtad *asqtad,
				   QDP_D3_ColorMatrix *fatlinks[],
				   QDP_D3_ColorMatrix *longlinks[]);

void QOP_D3_asqtad_dslash_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      double mass,
			      QDP_D3_ColorVector *out,
			      QDP_D3_ColorVector *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_D3_asqtad_diaginv_qdp(QOP_info_t *info,
			       QOP_D3_FermionLinksAsqtad *asqtad,
			       double mass,
			       QDP_D3_ColorVector *out,
			       QDP_D3_ColorVector *in,
			       QOP_evenodd_t eo);

void QOP_D3_asqtad_invert_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksAsqtad *asqtad,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      double mass,
			      QDP_D3_ColorVector *out,
			      QDP_D3_ColorVector *in);

void QOP_D3_asqtad_invert_multi_qdp(QOP_info_t *info,
				    QOP_D3_FermionLinksAsqtad *asqtad,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    double *masses[],
				    int nmass[],
				    QDP_D3_ColorVector **out[],
				    QDP_D3_ColorVector *in[],
				    int nsrc);

void QOP_D3_asqtad_get_eigcg(QOP_D3_FermionLinksAsqtad *asqtad,
			     QDP_D3_ColorVector **evecs,
			     QLA_D_Real *evals, int *nv);


  /*******************/
  /*  HISQ routines  */
  /*******************/

  /* single precision */

void QOP_F3_hisq_dslash_qdp(QOP_info_t *info,
			    QOP_F3_FermionLinksHisq *hisq,
			    float mass,
			    QDP_F3_ColorVector *out,
			    QDP_F3_ColorVector *in,
			    QOP_evenodd_t eo_out,
			    QOP_evenodd_t eo_in);

void QOP_F3_hisq_diaginv_qdp(QOP_info_t *info,
			     QOP_F3_FermionLinksHisq *hisq,
			     float mass,
			     QDP_F3_ColorVector *out,
			     QDP_F3_ColorVector *in,
			     QOP_evenodd_t eo);

void QOP_F3_hisq_invert_qdp(QOP_info_t *info,
			    QOP_F3_FermionLinksHisq *hisq,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t *res_arg,
			    float mass,
			    QDP_F3_ColorVector *out,
			    QDP_F3_ColorVector *in);

void QOP_F3_hisq_invert_multi_qdp(QOP_info_t *info,
				  QOP_F3_FermionLinksHisq *hisq,
				  QOP_invert_arg_t *inv_arg,
				  QOP_resid_arg_t **res_arg[],
				  float *masses[],
				  int nmass[],
				  QDP_F3_ColorVector **out[],
				  QDP_F3_ColorVector *in[],
				  int nsrc);

  /* double precision */

void QOP_D3_hisq_dslash_qdp(QOP_info_t *info,
			    QOP_D3_FermionLinksHisq *hisq,
			    double mass,
			    QDP_D3_ColorVector *out,
			    QDP_D3_ColorVector *in,
			    QOP_evenodd_t eo_out,
			    QOP_evenodd_t eo_in);

void QOP_D3_hisq_diaginv_qdp(QOP_info_t *info,
			     QOP_D3_FermionLinksHisq *hisq,
			     double mass,
			     QDP_D3_ColorVector *out,
			     QDP_D3_ColorVector *in,
			     QOP_evenodd_t eo);

void QOP_D3_hisq_invert_qdp(QOP_info_t *info,
			    QOP_D3_FermionLinksHisq *hisq,
			    QOP_invert_arg_t *inv_arg,
			    QOP_resid_arg_t *res_arg,
			    double mass,
			    QDP_D3_ColorVector *out,
			    QDP_D3_ColorVector *in);

void QOP_D3_hisq_invert_multi_qdp(QOP_info_t *info,
				  QOP_D3_FermionLinksHisq *hisq,
				  QOP_invert_arg_t *inv_arg,
				  QOP_resid_arg_t **res_arg[],
				  double *masses[],
				  int nmass[],
				  QDP_D3_ColorVector **out[],
				  QDP_D3_ColorVector *in[],
				  int nsrc);


  /*********************/
  /*  Wilson routines  */
  /*********************/

  /* single precision */

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_create_L_from_qdp(QDP_F3_ColorMatrix *links[],
				  QDP_F3_DiracPropagator *clov);

void QOP_F3_wilson_extract_L_to_qdp(QDP_F3_ColorMatrix *links[],
				    QDP_F3_DiracPropagator *clov,
				    QOP_F3_FermionLinksWilson *src);

QOP_F3_FermionLinksWilson *
  QOP_F3_wilson_convert_L_from_qdp(QDP_F3_ColorMatrix *links[],
				   QDP_F3_DiracPropagator *clov);

void QOP_F3_wilson_convert_L_to_qdp(QDP_F3_ColorMatrix ***links,
				    QDP_F3_DiracPropagator **clov,
				    QOP_F3_FermionLinksWilson *src);

void QOP_F3_wilson_load_L_from_qdp(QOP_F3_FermionLinksWilson *wilson,
				   QDP_F3_ColorMatrix *links[],
				   QDP_F3_DiracPropagator *clov);

void QOP_F3_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_F3_FermionLinksWilson *flw,
			      float kappa,
			      int sign,
			      QDP_F3_DiracFermion *out,
			      QDP_F3_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_F3_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_F3_FermionLinksWilson *flw,
			       float kappa,
			       QDP_F3_DiracFermion *out,
			       QDP_F3_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_F3_wilson_invert_qdp(QOP_info_t *info,
			      QOP_F3_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      float kappa,
			      QDP_F3_DiracFermion *out_pt,
			      QDP_F3_DiracFermion *in_pt);

void QOP_F3_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_F3_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    float *kappas[],
				    int nkappa[],
				    QDP_F3_DiracFermion **out_pt[],
				    QDP_F3_DiracFermion *in_pt[],
				    int nsrc);

  /* double precision */

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_create_L_from_qdp(QDP_D3_ColorMatrix *links[],
				  QDP_D3_DiracPropagator *clov);

void QOP_D3_wilson_extract_L_to_qdp(QDP_D3_ColorMatrix *links[],
				    QDP_D3_DiracPropagator *clov,
				    QOP_D3_FermionLinksWilson *src);

QOP_D3_FermionLinksWilson *
  QOP_D3_wilson_convert_L_from_qdp(QDP_D3_ColorMatrix *links[],
				   QDP_D3_DiracPropagator *clov);

void QOP_D3_wilson_convert_L_to_qdp(QDP_D3_ColorMatrix ***links,
				    QDP_D3_DiracPropagator **clov,
				    QOP_D3_FermionLinksWilson *src);

void QOP_D3_wilson_load_L_from_qdp(QOP_D3_FermionLinksWilson *wilson,
				   QDP_D3_ColorMatrix *links[],
				   QDP_D3_DiracPropagator *clov);

void QOP_D3_wilson_dslash_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksWilson *flw,
			      double kappa,
			      int sign,
			      QDP_D3_DiracFermion *out,
			      QDP_D3_DiracFermion *in,
			      QOP_evenodd_t eo_out,
			      QOP_evenodd_t eo_in);

void QOP_D3_wilson_diaginv_qdp(QOP_info_t *info,
			       QOP_D3_FermionLinksWilson *flw,
			       double kappa,
			       QDP_D3_DiracFermion *out,
			       QDP_D3_DiracFermion *in,
			       QOP_evenodd_t eo);

void QOP_D3_wilson_invert_qdp(QOP_info_t *info,
			      QOP_D3_FermionLinksWilson *links,
			      QOP_invert_arg_t *inv_arg,
			      QOP_resid_arg_t *res_arg,
			      double kappa,
			      QDP_D3_DiracFermion *out_pt,
			      QDP_D3_DiracFermion *in_pt);

void QOP_D3_wilson_invert_multi_qdp(QOP_info_t *info,
				    QOP_D3_FermionLinksWilson *links,
				    QOP_invert_arg_t *inv_arg,
				    QOP_resid_arg_t **res_arg[],
				    double *kappas[],
				    int nkappa[],
				    QDP_D3_DiracFermion **out_pt[],
				    QDP_D3_DiracFermion *in_pt[],
				    int nsrc);


  /**************************/
  /*  Domain Wall routines  */
  /**************************/

// Note that in our convention, free-field M5 = 1

  /* single precision */

QOP_F3_FermionLinksDW *
  QOP_F3_dw_create_L_from_qdp(QDP_F3_ColorMatrix *links[]);

void QOP_F3_dw_extract_L_to_qdp(QDP_F3_ColorMatrix *links[],
				QOP_F3_FermionLinksDW *src);

QOP_F3_FermionLinksDW *
  QOP_F3_dw_convert_L_from_qdp(QDP_F3_ColorMatrix *links[]);

void QOP_F3_dw_convert_L_to_qdp(QDP_F3_ColorMatrix ***links,
				QOP_F3_FermionLinksDW *src);

void QOP_F3_dw_load_L_from_qdp(QOP_F3_FermionLinksDW *dw,
			       QDP_F3_ColorMatrix *links[]);
			       
void QOP_F3_dw_dslash_qdp(QOP_info_t *info,
			  QOP_F3_FermionLinksDW *links,
			  float M5, 
			  float m,
			  int sign,
			  QDP_F3_DiracFermion *out_pt[],
			  QDP_F3_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);
			  
void QOP_F3_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_F3_FermionLinksDW *links,
			   float M5,
			   float m,
			   QDP_F3_DiracFermion *out_pt[],
			   QDP_F3_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);
			       
void QOP_F3_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_FermionLinksDW *fldw,
			   float M5,
			   float m,
			   QDP_DiracFermion **out,
			   QDP_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_F3_dw_invert_qdp(QOP_info_t *info,
			  QOP_F3_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  float M5,
			  float m,
			  QDP_F3_DiracFermion *out[],
			  QDP_F3_DiracFermion *in[],
			  int Ls);

void QOP_F3_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_F3_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				float *M5[],
				float *m[],
				int nmass[],
				QDP_F3_DiracFermion ***out[],
				QDP_F3_DiracFermion **in[],
				int nsrc,
				int Ls);


  /* double precision */

QOP_D3_FermionLinksDW *
  QOP_D3_dw_create_L_from_qdp(QDP_D3_ColorMatrix *links[]);

void QOP_D3_dw_extract_L_to_qdp(QDP_D3_ColorMatrix *links[],
				QOP_D3_FermionLinksDW *src);

QOP_D3_FermionLinksDW *
  QOP_D3_dw_convert_L_from_qdp(QDP_D3_ColorMatrix *links[]);

void QOP_D3_dw_convert_L_to_qdp(QDP_D3_ColorMatrix ***links,
				QOP_D3_FermionLinksDW *src);

void QOP_D3_dw_load_L_from_qdp(QOP_D3_FermionLinksDW *dw,
			       QDP_D3_ColorMatrix *links[]);

void QOP_D3_dw_dslash_qdp(QOP_info_t *info,
			  QOP_D3_FermionLinksDW *links,
			  double M5,
			  double m,
			  int sign,
			  QDP_D3_DiracFermion *out_pt[],
			  QDP_D3_DiracFermion *in_pt[],
			  int Ls,
			  QOP_evenodd_t eo_out,
			  QOP_evenodd_t eo_in);
			  
void QOP_D3_dw_dslash2_qdp(QOP_info_t *info,
			   QOP_D3_FermionLinksDW *links,
			   double M5,
			   double m,
			   QDP_D3_DiracFermion *out_pt[],
			   QDP_D3_DiracFermion *in_pt[],
			   int Ls,
			   QOP_evenodd_t eo_out,
			   QOP_evenodd_t eo_in);
			  
void QOP_D3_dw_diaginv_qdp(QOP_info_t *info,
			   QOP_FermionLinksDW *fldw,
			   double M5,
			   double m,
			   QDP_DiracFermion **out,
			   QDP_DiracFermion **in,
			   int ls,
			   QOP_evenodd_t eo);

void QOP_D3_dw_invert_qdp(QOP_info_t *info,
			  QOP_D3_FermionLinksDW *links,
			  QOP_invert_arg_t *inv_arg,
			  QOP_resid_arg_t *res_arg,
			  double M5,
			  double m,
			  QDP_D3_DiracFermion *out_pt[],
			  QDP_D3_DiracFermion *in_pt[],
			  int Ls);

void QOP_D3_dw_invert_multi_qdp(QOP_info_t *info,
				QOP_D3_FermionLinksDW *links,
				QOP_invert_arg_t *inv_arg,
				QOP_resid_arg_t **res_arg[],
				double *M5[],
				double *m[],
				int nmass[],
				QDP_D3_DiracFermion ***out_pt[],
				QDP_D3_DiracFermion **in_pt[],
				int nsrc,
				int Ls);


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
#define QOP_asqtad_load_L_from_qdp    QOP_F3_asqtad_load_L_from_qdp
#define QOP_asqtad_dslash_qdp         QOP_F3_asqtad_dslash_qdp
#define QOP_asqtad_diaginv_qdp        QOP_F3_asqtad_diaginv_qdp
#define QOP_asqtad_invert_qdp         QOP_F3_asqtad_invert_qdp
#define QOP_asqtad_invert_multi_qdp   QOP_F3_asqtad_invert_multi_qdp
#define QOP_asqtad_get_eigcg          QOP_F3_asqtad_get_eigcg

#define QOP_hisq_dslash_qdp         QOP_F3_hisq_dslash_qdp
#define QOP_hisq_diaginv_qdp        QOP_F3_hisq_diaginv_qdp
#define QOP_hisq_invert_qdp         QOP_F3_hisq_invert_qdp
#define QOP_hisq_invert_multi_qdp   QOP_F3_hisq_invert_multi_qdp

#define QOP_wilson_create_L_from_qdp  QOP_F3_wilson_create_L_from_qdp
#define QOP_wilson_extract_L_to_qdp   QOP_F3_wilson_extract_L_to_qdp
#define QOP_wilson_convert_L_from_qdp QOP_F3_wilson_convert_L_from_qdp
#define QOP_wilson_convert_L_to_qdp   QOP_F3_wilson_convert_L_to_qdp
#define QOP_wilson_load_L_from_qdp    QOP_F3_wilson_load_L_from_qdp
#define QOP_wilson_dslash_qdp         QOP_F3_wilson_dslash_qdp
#define QOP_wilson_diaginv_qdp        QOP_F3_wilson_diaginv_qdp
#define QOP_wilson_invert_qdp         QOP_F3_wilson_invert_qdp
#define QOP_wilson_invert_multi_qdp   QOP_F3_wilson_invert_multi_qdp

#define QOP_dw_create_L_from_qdp  QOP_F3_dw_create_L_from_qdp
#define QOP_dw_extract_L_to_qdp   QOP_F3_dw_extract_L_to_qdp
#define QOP_dw_convert_L_from_qdp QOP_F3_dw_convert_L_from_qdp
#define QOP_dw_convert_L_to_qdp   QOP_F3_dw_convert_L_to_qdp
#define QOP_dw_load_L_from_qdp    QOP_F3_dw_load_L_from_qdp
#define QOP_dw_dslash_qdp         QOP_F3_dw_dslash_qdp
#define QOP_dw_dslash2_qdp         QOP_F3_dw_dslash2_qdp
#define QOP_dw_diaginv_qdp        QOP_F3_dw_diaginv_qdp
#define QOP_dw_invert_qdp         QOP_F3_dw_invert_qdp
#define QOP_dw_invert_multi_qdp   QOP_F3_dw_invert_multi_qdp

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
#define QOP_asqtad_load_L_from_qdp    QOP_D3_asqtad_load_L_from_qdp
#define QOP_asqtad_dslash_qdp         QOP_D3_asqtad_dslash_qdp
#define QOP_asqtad_diaginv_qdp        QOP_D3_asqtad_diaginv_qdp
#define QOP_asqtad_invert_qdp         QOP_D3_asqtad_invert_qdp
#define QOP_asqtad_invert_multi_qdp   QOP_D3_asqtad_invert_multi_qdp
#define QOP_asqtad_get_eigcg          QOP_D3_asqtad_get_eigcg

#define QOP_hisq_dslash_qdp         QOP_D3_hisq_dslash_qdp
#define QOP_hisq_diaginv_qdp        QOP_D3_hisq_diaginv_qdp
#define QOP_hisq_invert_qdp         QOP_D3_hisq_invert_qdp
#define QOP_hisq_invert_multi_qdp   QOP_D3_hisq_invert_multi_qdp

#define QOP_wilson_create_L_from_qdp  QOP_D3_wilson_create_L_from_qdp
#define QOP_wilson_extract_L_to_qdp   QOP_D3_wilson_extract_L_to_qdp
#define QOP_wilson_convert_L_from_qdp QOP_D3_wilson_convert_L_from_qdp
#define QOP_wilson_convert_L_to_qdp   QOP_D3_wilson_convert_L_to_qdp
#define QOP_wilson_load_L_from_qdp    QOP_D3_wilson_load_L_from_qdp
#define QOP_wilson_dslash_qdp         QOP_D3_wilson_dslash_qdp
#define QOP_wilson_diaginv_qdp        QOP_D3_wilson_diaginv_qdp
#define QOP_wilson_invert_qdp         QOP_D3_wilson_invert_qdp
#define QOP_wilson_invert_multi_qdp   QOP_D3_wilson_invert_multi_qdp

#define QOP_dw_create_L_from_qdp  QOP_D3_dw_create_L_from_qdp
#define QOP_dw_extract_L_to_qdp   QOP_D3_dw_extract_L_to_qdp
#define QOP_dw_convert_L_from_qdp QOP_D3_dw_convert_L_from_qdp
#define QOP_dw_convert_L_to_qdp   QOP_D3_dw_convert_L_to_qdp
#define QOP_dw_load_L_from_qdp    QOP_D3_dw_load_L_from_qdp
#define QOP_dw_dslash_qdp         QOP_D3_dw_dslash_qdp
#define QOP_dw_dslash2_qdp        QOP_D3_dw_dslash2_qdp
#define QOP_dw_diaginv_qdp        QOP_D3_dw_diaginv_qdp
#define QOP_dw_invert_qdp         QOP_D3_dw_invert_qdp
#define QOP_dw_invert_multi_qdp   QOP_D3_dw_invert_multi_qdp

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_QDP_H */
