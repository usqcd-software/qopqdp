#include <qop_internal.h>
#include <string.h>

int QOP_wilson_inited = 0;
int QOP_wilson_style = 0;
int QOP_wilson_nsvec = 4;
int QOP_wilson_nvec = 4;

QOP_status_t
QOP_asqtad_invert_init(void)
{
  QOP_wilson_inited = 1;
  return QOP_SUCCESS;
}

QOP_status_t
QOP_wilson_invert_finalize(void)
{
  QOP_wilson_inited = 0;
  return QOP_SUCCESS;
}

QOP_status_t
QOP_wilson_invert_set_opts(QOP_opt_t opts[], int nopts)
{
  int i;
  for(i=0; i<nopts; i++) {
    char *tag = opts[i].tag;
    double value = opts[i].value;

    if(!strcmp(tag,"ns")) {
      if( (value==1)||(value==2)||(value==4)||
	  (((QOP_wilson_style&1)==1)&&(value==8)) ) {
	QOP_wilson_nsvec = (int) value;
      } else {
	return QOP_FAIL;
      }
    } else if(!strcmp(tag,"nm")) {
      if( (value==1)||(value==2)||(value==4)||
	  (((QOP_wilson_style&1)==1)&&(value==8)) ) {
	QOP_wilson_nvec = (int) value;
      } else {
	return QOP_FAIL;
      }
    } else if(!strcmp(tag,"st")) {
      if((value==0)||(value==1)||(value==2)||(value==3)) {
	QOP_wilson_style = (int) value;
	if((QOP_wilson_style&1)==0) {
	  if(QOP_wilson_nsvec>4) QOP_wilson_nsvec = 4;
	  if(QOP_wilson_nvec>4) QOP_wilson_nvec = 4;
	}
      } else {
	return QOP_FAIL;
      }
    }
  }

  return QOP_SUCCESS;
}
