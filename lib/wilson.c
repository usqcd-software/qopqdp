#include <qop_internal.h>
#include <string.h>

int QOP_wilson_inited = 0;
int QOP_wilson_style = 0;
int QOP_wilson_nsvec = 8;
int QOP_wilson_nvec = 8;

QOP_status_t
QOP_wilson_invert_init(QOP_layout *layout) {
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
QOP_wilson_set_opt(char *tag, double value) {
  if((value==0)||(value==2)) return QOP_SUCCESS;
  return QOP_FAIL;
}
