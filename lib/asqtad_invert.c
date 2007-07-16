#include <string.h>
#include <qop_internal.h>

QOP_asqtad_t QOP_asqtad = {.inited=0};

#define OPTNUM (QOP_asqtad.style+2*QOP_asqtad.nsvec)

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
  QOP_asqtad.optnum = OPTNUM;
  return QOP_SUCCESS;
}

#define setvar(_var, _type, _tag, _opts, _nopts)			\
  { int i; for(i=0; i<_nopts; i++) {					\
      if(!strcmp(_opts[i].tag,_tag)) _var = (_type) _opts[i].value;	\
    } }

#define valid_style(st) ( (st>=0) && (st<=3) )
#define valid_nvec(nv,st) \
  ((nv==1)||(nv==2)||(nv==4)||(nv==8)||((st==1)&&(nv==16)))
#define fix_nvec(nv,st) { if(nv<1) nv=1; else nv=8; }

QOP_status_t
QOP_asqtad_invert_set_opts(QOP_opt_t opts[], int nopts)
{
  int st, ns, nm;
  st = QOP_asqtad.style;
  ns = QOP_asqtad.nsvec;
  nm = QOP_asqtad.nvec;

  setvar(st, int, "st", opts, nopts);
  if(!valid_style(st)) return QOP_FAIL;

  if(!valid_nvec(ns,st)) fix_nvec(ns,st);
  if(!valid_nvec(nm,st)) fix_nvec(nm,st);

  setvar(ns, int, "ns", opts, nopts);
  setvar(nm, int, "nm", opts, nopts);
  if(!valid_nvec(ns,st)) return QOP_FAIL;
  if(!valid_nvec(nm,st)) return QOP_FAIL;

  QOP_asqtad.style = st;
  QOP_asqtad.nsvec = ns;
  QOP_asqtad.nvec = nm;
  QOP_asqtad.optnum = OPTNUM;

  return QOP_SUCCESS;
}
