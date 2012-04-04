#include <string.h>
#include <qop_internal.h>

QOP_hisq_links_t QOP_hisq_links = {
  .inited = 0,
  .want_deps = 0,
  .want_aux=1,
  .reunit_allow_svd = 1,
  .reunit_svd_only = 0,
  .reunit_svd_rel_error = 1e-8,
  .reunit_svd_abs_error = 1e-8,
  .svd_values_info = 1,
  .use_fat7_lepage = 0
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
  int reunit_allow_svd, reunit_svd_only;
  double reunit_svd_rel_error, reunit_svd_abs_error;
  int svd_values_info, use_fat7_lepage;

  HISQ_LINKS_BEGIN;

  want_deps = QOP_hisq_links.want_deps;
  setvar(want_deps, int, "want_deps", opts, nopts);
  QOP_hisq_links.want_deps = want_deps;

  want_aux = QOP_hisq_links.want_aux;
  setvar(want_aux, int, "want_aux", opts, nopts);
  QOP_hisq_links.want_aux = want_aux;

  reunit_allow_svd = QOP_hisq_links.reunit_allow_svd;
  setvar(reunit_allow_svd, int, "reunit_allow_svd", opts, nopts);
  QOP_hisq_links.reunit_allow_svd = reunit_allow_svd;

  reunit_svd_only = QOP_hisq_links.reunit_svd_only;
  setvar(reunit_svd_only, int, "reunit_svd_only", opts, nopts);
  QOP_hisq_links.reunit_svd_only = reunit_svd_only;

  reunit_svd_rel_error = QOP_hisq_links.reunit_svd_rel_error;
  setvar(reunit_svd_rel_error, double, "reunit_svd_rel_error", opts, nopts);
  QOP_hisq_links.reunit_svd_rel_error = reunit_svd_rel_error;

  reunit_svd_abs_error = QOP_hisq_links.reunit_svd_abs_error;
  setvar(reunit_svd_abs_error, double, "reunit_svd_abs_error", opts, nopts);
  QOP_hisq_links.reunit_svd_abs_error = reunit_svd_abs_error;

  svd_values_info = QOP_hisq_links.svd_values_info;
  setvar(svd_values_info, int, "svd_values_info", opts, nopts);
  QOP_hisq_links.svd_values_info = svd_values_info;

  use_fat7_lepage = QOP_hisq_links.use_fat7_lepage;
  setvar(use_fat7_lepage, int, "use_fat7_lepage", opts, nopts);
  QOP_hisq_links.use_fat7_lepage = use_fat7_lepage;

  HISQ_LINKS_END;

  return QOP_SUCCESS;
}
