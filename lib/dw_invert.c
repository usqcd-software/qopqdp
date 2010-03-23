#include <qop_internal.h>
#include <string.h>

int QOP_dw_initQ = 0;
int QOP_dw_style = 0;
int QOP_dw_nsvec = 4;
int QOP_dw_nvec = 4;
extern int QOP_wilson_style;
extern int QOP_wilson_nsvec;
extern int QOP_wilson_nvec;

QOP_status_t
QOP_dw_invert_set_opts(QOP_opt_t opts[], int nopts)
{
  int i;

  DW_INVERT_BEGIN;

  for (i=0; i<nopts; i++) {
    char *tag = opts[i].tag;
    double value = opts[i].value;

    if (!strcmp(tag,"ns")) {
      if ( value==1 || value==2 || value==4 ||
	         ((QOP_dw_style&1)==1 && value==8) ) {
	      QOP_dw_nsvec = (int) value;
      } else {
	      return QOP_FAIL;
      }
    } else if (!strcmp(tag,"nm")) {
      if ( value==1 || value==2 || value==4 ||
	         ((QOP_dw_style&1)==1 && value==8) ) {
	      QOP_dw_nvec = (int) value;
      } else {
	      return QOP_FAIL;
      }
    } else if (!strcmp(tag,"st")) {
      if ( value==0 || value==1 || value==2 || value==3 ) {
	      QOP_dw_style = (int) value;
	      if ((QOP_dw_style&1)==0) {
	        if (QOP_dw_nsvec>4) QOP_dw_nsvec = 4;
	        if (QOP_dw_nvec>4) QOP_dw_nvec = 4;
	      }
      } else {
	      return QOP_FAIL;
      }
    }
  }
  
  // Need to set the Wilson style for use in that operator
  QOP_wilson_style = QOP_dw_style;
  QOP_wilson_nsvec = QOP_dw_nsvec;
  QOP_wilson_nvec  = QOP_dw_nvec;

  DW_INVERT_END;
  return QOP_SUCCESS;
}
