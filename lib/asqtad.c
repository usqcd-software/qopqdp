#include <qop.h>
#include <string.h>

int QOP_asqtad_inited = 0;
int QOP_style = 0;
int QOP_nsvec = 8;
int QOP_nvec = 8;
QDP_Shift QOP_asqtad_shifts[8], QOP_neighbor3[4];
QDP_Shift QOP_asqtad_shifts_dbl[16];
QDP_ShiftDir QOP_shiftfwd[8], QOP_shiftbck[8];
QDP_ShiftDir QOP_asqtad_shiftdirs_dbl[16];

QOP_status_t
QOP_asqtad_invert_init(QOP_layout *layout)
{
  int i;
  int disp[4]={0,0,0,0};

  for(i=0; i<4; i++) {
    QOP_asqtad_shifts[2*i] = QDP_neighbor[i];
    disp[i] = 3;
    QOP_neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
    QOP_asqtad_shifts[2*i+1] = QOP_neighbor3[i];
  }
  for(i=0; i<8; i++) {
    QOP_shiftfwd[i] = QDP_forward;
    QOP_shiftbck[i] = QDP_backward;
  }
  for(i=0; i<4; i++) {
    QOP_asqtad_shifts_dbl[4*i] = QOP_asqtad_shifts[2*i];
    QOP_asqtad_shifts_dbl[4*i+1] = QOP_asqtad_shifts[2*i+1];
    QOP_asqtad_shifts_dbl[4*i+2] = QOP_asqtad_shifts[2*i];
    QOP_asqtad_shifts_dbl[4*i+3] = QOP_asqtad_shifts[2*i+1];
    QOP_asqtad_shiftdirs_dbl[4*i] = QDP_forward;
    QOP_asqtad_shiftdirs_dbl[4*i+1] = QDP_forward;
    QOP_asqtad_shiftdirs_dbl[4*i+2] = QDP_backward;
    QOP_asqtad_shiftdirs_dbl[4*i+3] = QDP_backward;
  }

  QOP_asqtad_inited = 1;
  return QOP_SUCCESS;
}

QOP_status_t
QOP_asqtad_invert_finalize(void)
{
  int i;
  for(i=0; i<4; i++) {
    QDP_destroy_shift(QOP_neighbor3[i]);
  }

  QOP_asqtad_inited = 0;
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
