#include <qop_internal.h>
#include <string.h>

QOP_asqtad_t QOP_asqtad = {.inited=0};

QOP_status_t
QOP_asqtad_invert_init(void)
{
  int i;
  int disp[4]={0,0,0,0};

  for(i=0; i<4; i++) {
    QOP_asqtad.shifts[2*i] = QDP_neighbor[i];
    disp[i] = 3;
    QOP_common.neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
    QOP_asqtad.shifts[2*i+1] = QOP_common.neighbor3[i];
  }
  for(i=0; i<8; i++) {
    QOP_common.shiftfwd[i] = QDP_forward;
    QOP_common.shiftbck[i] = QDP_backward;
  }
  for(i=0; i<4; i++) {
    QOP_asqtad.shifts_dbl[4*i] = QOP_asqtad.shifts[2*i];
    QOP_asqtad.shifts_dbl[4*i+1] = QOP_asqtad.shifts[2*i+1];
    QOP_asqtad.shifts_dbl[4*i+2] = QOP_asqtad.shifts[2*i];
    QOP_asqtad.shifts_dbl[4*i+3] = QOP_asqtad.shifts[2*i+1];
    QOP_asqtad.shiftdirs_dbl[4*i] = QDP_forward;
    QOP_asqtad.shiftdirs_dbl[4*i+1] = QDP_forward;
    QOP_asqtad.shiftdirs_dbl[4*i+2] = QDP_backward;
    QOP_asqtad.shiftdirs_dbl[4*i+3] = QDP_backward;
  }

  QOP_asqtad.inited = 1;
  QOP_asqtad.style = 1;
  QOP_asqtad.nsvec = 8;
  QOP_asqtad.nvec = 8;
  return QOP_SUCCESS;
}

QOP_status_t
QOP_asqtad_invert_set_opts(QOP_opt_t opts[], int nopts)
{
  int i;
  for(i=0; i<nopts; i++) {
    char *tag = opts[i].tag;
    double value = opts[i].value;

    if(!strcmp(tag,"ns")) {
      if( (value==1)||(value==2)||(value==4)||(value==8)||
	  ((QOP_asqtad.style==1)&&(value==16)) ) {
	QOP_asqtad.nsvec = (int) value;
      } else {
	return QOP_FAIL;
      }
    } else if(!strcmp(tag,"nm")) {
      if( (value==1)||(value==2)||(value==4)||(value==8)||
	  ((QOP_asqtad.style==1)&&(value==16)) ) {
	QOP_asqtad.nvec = (int) value;
      } else {
	return QOP_FAIL;
      }
    } else if(!strcmp(tag,"st")) {
      if((value==0)||(value==1)) {
	QOP_asqtad.style = (int) value;
	if(QOP_asqtad.style==0) {
	  if(QOP_asqtad.nsvec>8) QOP_asqtad.nsvec = 8;
	  if(QOP_asqtad.nvec>8) QOP_asqtad.nvec = 8;
	}
      } else {
	return QOP_FAIL;
      }
    } else {
      return QOP_FAIL;
    }
  }

  return QOP_SUCCESS;
}
