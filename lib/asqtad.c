#include <qop.h>

int QOP_style = 0;
int QOP_nsvec = 8;
int QOP_nvec = 8;

QOP_status_t
QOP_asqtad_invert_init(QOP_layout *layout)
{
  return QOP_SUCCESS;
}

QOP_status_t
QOP_asqtad_invert_finalize(void)
{
  return QOP_SUCCESS;
}

QOP_status_t
QOP_asqtad_set_opt(char *tag, double value)
{
  if(!strcmp(tag,"ns")) {
    if( (value==1)||(value==2)||(value==4)||(value==8)||
	((QOP_style==1)&&(value==16)) ) {
      QOP_nsvec = (int) value;
      return QOP_SUCCESS;
    }
  } else if(!strcmp(tag,"nm")) {
    if( (value==1)||(value==2)||(value==4)||(value==8)||
	((QOP_style==1)&&(value==16)) ) {
      QOP_nvec = (int) value;
      return QOP_SUCCESS;
    }
  } else if(!strcmp(tag,"st")) {
    if((value==0)||(value==1)) {
      QOP_style = (int) value;
      return QOP_SUCCESS;
    }
  }

  return QOP_FAIL;
}
