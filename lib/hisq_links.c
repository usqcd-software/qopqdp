#include <string.h>
#include <qop_internal.h>

QOP_hisq_links_t QOP_hisq_links = {
  .inited = 0,
  .want_deps = 0,
  .want_aux=1
};

#define setvar(_var, _type, _tag, _opts, _nopts)			\
  { int i; for(i=0; i<_nopts; i++) {					\
      if(!strcmp(_opts[i].tag,_tag)) _var = (_type) _opts[i].value;	\
    } }

/* Options are these

   want_deps  true if we want the derivative wrto Naik eps
   want_aux   true if we want to keep the auxiliary links
              (We need to keep them for a fermion force
	      calculation, but if we need the links only for
	      inversions, then we can discard them to save space.)

*/

QOP_status_t
QOP_hisq_links_set_opts(QOP_opt_t opts[], int nopts)
{
  ;
  int want_deps, want_aux;

  HISQ_LINKS_BEGIN;

  want_deps = QOP_hisq_links.want_deps;
  setvar(want_deps, int, "want_deps", opts, nopts);
  QOP_hisq_links.want_deps = want_deps;

  want_aux = QOP_hisq_links.want_aux;
  setvar(want_aux, int, "want_aux", opts, nopts);
  QOP_hisq_links.want_aux = want_aux;

  HISQ_LINKS_END;

  return QOP_SUCCESS;
}
