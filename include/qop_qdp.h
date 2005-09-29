#ifndef _QOP_QDP_H
#define _QOP_QDP_H

#include <qop.h>
#include <qdp.h>

#ifdef __cplusplus
extern "C" {
#endif


  /*  Asqtad routines  */

QOP_status_t QOP_F_asqtad_invert_load_links_qdp(QDP_F3_ColorMatrix *fatlnk[],
                                                QDP_F3_ColorMatrix *longlnk[]);
int QOP_F_asqtad_invert_qdp(QOP_invert_arg *inv_arg,
			    QDP_F3_ColorVector *out, QDP_F3_ColorVector *in);

QOP_status_t QOP_D_asqtad_invert_load_links_qdp(QDP_D3_ColorMatrix *fatlnk[],
                                                QDP_D3_ColorMatrix *longlnk[]);
int QOP_D_asqtad_invert_qdp(QOP_invert_arg *inv_arg,
			    QDP_D3_ColorVector *out, QDP_D3_ColorVector *in);


  /*  Wilson routines  */

QOP_status_t QOP_F_wilson_invert_load_links_qdp(QDP_F3_ColorMatrix *links[]);
int QOP_F_wilson_invert_qdp(QOP_invert_arg *inv_arg,
			    QDP_F3_DiracFermion *out, QDP_F3_DiracFermion *in);

QOP_status_t QOP_D_wilson_invert_load_links_qdp(QDP_D3_ColorMatrix *links[]);
int QOP_D_wilson_invert_qdp(QOP_invert_arg *inv_arg,
			    QDP_D3_DiracFermion *out, QDP_D3_DiracFermion *in);


  /* Mapping of generic names to specific precision */

#if QOP_Precision == 1

#define QOP_asqtad_invert_load_links_qdp QOP_F_asqtad_invert_load_links_qdp
#define QOP_asqtad_invert_qdp QOP_F_asqtad_invert_qdp

#define QOP_wilson_invert_load_links_qdp QOP_F_wilson_invert_load_links_qdp
#define QOP_wilson_invert_qdp QOP_F_wilson_invert_qdp

#else

#define QOP_asqtad_invert_load_links_qdp QOP_D_asqtad_invert_load_links_qdp
#define QOP_asqtad_invert_qdp QOP_D_asqtad_invert_qdp

#define QOP_wilson_invert_load_links_qdp QOP_D_wilson_invert_load_links_qdp
#define QOP_wilson_invert_qdp QOP_D_wilson_invert_qdp

#endif


#ifdef __cplusplus
}
#endif

#endif /* _QOP_QDP_H */
