#include <string.h>
#include <qop_internal.h>

QOP_asqtad_force_t QOP_asqtad_ff = {.inited=0,.style=1};

#define setvar(_var, _type, _tag, _opts, _nopts)			\
  { int i; for(i=0; i<_nopts; i++) {					\
      if(!strcmp(_opts[i].tag,_tag)) _var = (_type) _opts[i].value;	\
    } }

/* Style choices are

   0 = ASVEC  (parallel transport source vectors)
   1 = FNMAT  (parallel transport outer product - best for many sources )

   The style choice applies only to nsrc > 3.

*/

#define valid_style(st) ( (st>=0) && (st<=1) )

QOP_status_t
QOP_asqtad_force_set_opts(QOP_opt_t opts[], int nopts)
{
  int st, ns, nm;
  st = QOP_asqtad_ff.style;

  setvar(st, int, "st", opts, nopts);
  if(!valid_style(st)) return QOP_FAIL;

  QOP_asqtad_ff.style = st;

  return QOP_SUCCESS;
}
