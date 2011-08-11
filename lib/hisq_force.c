#include <string.h>
#include <qop_internal.h>

QOP_hisq_force_t QOP_hisq_ff = {
  .inited = 0,
  .fnmat_src_min = 4,
  .veclength=4
};

#define setvar(_var, _type, _tag, _opts, _nopts)			\
  { int i; for(i=0; i<_nopts; i++) {					\
      if(!strcmp(_opts[i].tag,_tag)) _var = (_type) _opts[i].value;	\
    } }

/* Options are these

   fnmat_src_min  For nsrc < value use the ASVEC algorithm.  Otherwise
                  use the FNMAT algorithm.
   veclength      The block size for ASVEC

*/

#define valid_fnmat_src_min(fsm) ( (fsm>=0) )
#define valid_veclength(vl) ( (vl>0) )

QOP_status_t
QOP_hisq_force_set_opts(QOP_opt_t opts[], int nopts)
{
  ;
//AB THIS ROUTINE IS NOT DOING ANYTHING AT THE MOMENT
#ifdef ABA
  int fsm, vl;

  HISQ_FORCE_BEGIN;

  fsm = QOP_asqtad_ff.fnmat_src_min;
  setvar(fsm, int, "fnmat_src_min", opts, nopts);
  if(!valid_fnmat_src_min(fsm)) return QOP_FAIL;
  QOP_asqtad_ff.fnmat_src_min = fsm;

  vl = QOP_asqtad_ff.veclength;
  setvar(vl, int, "veclength", opts, nopts);
  if(!valid_veclength(vl)) return QOP_FAIL;
  QOP_asqtad_ff.veclength = vl;

  HISQ_FORCE_END;
#endif
  return QOP_SUCCESS;
}
