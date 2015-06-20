/**************************************************************************
This code inverts the D-slash for the improved fermilab fermions 
It uses the bi-conjugate gradient algoritm (cgtype = 1). 
Adapted from James Osborn's "wilson_invert_p.c" code...........MBO/08/24/09
***************************************************************************/
#include <string.h>
#include <qop_internal.h>


#define IFLA
#define printf0(...)

extern int QOP_wilson_cgtype;
extern int QOP_wilson_eigcg_numax;
extern int QOP_wilson_eigcg_m;
extern int QOP_wilson_eigcg_nev;
extern QDP_DiracFermion *QOP_wilson_dslash_get_tmp(QOP_FermionLinksWilson *flw, QOP_evenodd_t eo, int n);

/* inverter stuff */

static QOP_FermionLinksWilson   *gl_flw;
static REAL                     gl_kappa;
static QOP_evenodd_t            gl_eo;
static QDP_DiracFermion         *gl_tmp, *gl_tmp2, *ctmp, *gtmp;
static QOP_wilson_ifla_coeffs_t tmc;

static void
QOP_wilson_invert_difla(QDP_DiracFermion *out, QDP_DiracFermion *in,
			   QDP_Subset subset)
{ 
  /* #ifdef IFLA  */
  QOP_wilson_ifla_dslash_qdp(NULL, gl_flw, gl_kappa, 1, &tmc, out, in, gl_eo, 
			     gl_eo);
  /* #else */
  //QOP_wilson_dslash_qdp(NULL, gl_flw, gl_kappa, 1, out, in, gl_eo, gl_eo);
  /* #endif */
}

//static int count = 1; // unused so comment out for now
	
void
QOP_wilson_ifla_invert(QOP_info_t *info,
		       QOP_FermionLinksWilson *flw,
		       QOP_invert_arg_t *inv_arg,
		       QOP_resid_arg_t *res_arg,
		       REAL kappa,
		       QOP_wilson_ifla_coeffs_t *coeffs,
		       QOP_DiracFermion *out,
		       QOP_DiracFermion *in)
{
  double dtime=0;
  double nflop;
  double rsqminold, relminold;
  QLA_Real rsq, rsqstop, relnorm2, insq;
  QDP_DiracFermion *qdpin, *qdpout;
  QDP_DiracFermion *cgp,   *cgr;
  
  QDP_Subset insub, cgsub;
  QOP_evenodd_t ineo, cgeo;

  QDP_Lattice *lat = QDP_get_lattice_M(flw->links[0]);
  
  int iter             =  0;
  int max_iter_old     = inv_arg->max_iter;    /* Max. # of iteration   */
  int max_restarts_old = inv_arg->max_restarts;/* # of restarts allowed */
  int nrestart = -1, max_restarts = inv_arg->max_restarts;
  if(max_restarts<=0) max_restarts = 5;

  //  res_arg->relmin = 1.e-14;

  WILSON_INVERT_BEGIN;
  
//#ifdef LU
//  printf("LU Decomposition : Yes  CgType = %d \n",QOP_wilson_cgtype );
//#else
//  printf("LU Decomposition : No   CgType = %d \n",QOP_wilson_cgtype );
//#endif
  
  
  /* -- IFLA COEFFICIENTS --- */
  tmc.kapifla = coeffs->kapifla;
  tmc.kappa_s = coeffs->kappa_s;
  tmc.kappa_t = coeffs->kappa_t;
  tmc.r_s     = coeffs->r_s;
  tmc.r_t     = coeffs->r_t;
  tmc.zeta    = coeffs->zeta;
  tmc.c_E     = coeffs->c_E;
  tmc.c_B     = coeffs->c_B;
  tmc.c_1     = coeffs->c_1;
  tmc.c_2     = coeffs->c_2;
  tmc.c_3     = coeffs->c_3;
  tmc.c_4     = coeffs->c_4;
  tmc.c_5     = coeffs->c_5;
  tmc.c_EE    = coeffs->c_EE;
  tmc.u0      = coeffs->u0;
  /* ------------------------ */

  
#if 0
  if(count == 0){
    printf("====================\n");
    printf("Function : QOP_wilson_ifla_invert\n");
    printf("These are the lattice action coefficients :: XX \n");
    printf("Kappa = %e\t r_s = %e\t zeta = %e\n",tmc.kapifla,tmc.r_s,tmc.zeta);
    printf("c_E   = %e\t c_B = %e\t c_EE = %e\n",tmc.c_E,tmc.c_B,tmc.c_EE);
    printf("c_1   = %e\t c_2 = %e\t c_3  = %e\n",tmc.c_1,tmc.c_2,tmc.c_3);
    printf("c_4   = %e\t c_5 = %e\t u0   = %e\n",tmc.c_4,tmc.c_5,tmc.u0);
    printf("====================\n");
  };
  count += count ;
#endif
  
  ineo   = inv_arg->evenodd; /* subset of source vectors */
  insub  = qdpsub(ineo);

  qdpin  = QDP_create_D_L(lat); 
  qdpout = QDP_create_D_L(lat);
  if(flw->clov!=NULL) ctmp = QDP_create_D_L(lat);
  

  gl_flw   = flw;   /* FermionLinksWilson */
  gl_kappa = kappa; /* Kappa Parameter    */

  /* cg has   5 * 48 = 240 flops/site/it */
  /* bicg has 9*4*24 = 864 flops/site/it */

  /* MdagM -> 2*((144+168*7)+48) = 2736 flops/site */
  /* clov -> 2*(12*(42-2)) = 960 flops/site */
  
  nflop  = 2736+240;
  if(QOP_wilson_cgtype==1) nflop += (864-240);
 
  cgeo = QOP_EVENODD;
  
  cgsub = qdpsub(cgeo);
  gl_eo = cgeo;

  cgp = QOP_wilson_dslash_get_tmp(flw, oppsub(cgeo), 1);
  cgr = QOP_wilson_dslash_get_tmp(flw, oppsub(cgeo), 2); 
 

  gl_tmp = cgr;

  QDP_D_eq_zero(qdpin, QDP_all);

  gl_tmp2 = cgr;
  if(QOP_wilson_cgtype==1) {
    QDP_D_eq_D(qdpin, in->df, insub);
  }; 
  
  QDP_D_eq_D(qdpout, out->df, insub);

 

  QDP_r_eq_norm2_D(&insq, in->df, insub);

  rsqstop = insq * res_arg->rsqmin;
  VERB(LOW, "WILSON_INVERT: rsqstop = %g\n", rsqstop);
  rsq        = 0;
  relnorm2   = 1.;
  rsqminold  = res_arg->rsqmin;
  relminold = res_arg->relmin;
  res_arg->rsqmin     *=  0.9;
  res_arg->relmin *= 0.5;
  inv_arg->max_restarts = 0;
  
  /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
  do {
    //    printf("Starting the loop\n");
    
    inv_arg->max_iter = max_iter_old - iter;

    dtime -= QOP_time();
    
    if(QOP_wilson_cgtype==1) {
      //      printf("Yes cgtype is 1\n");
      QOP_invert_bicgstab_D(QOP_wilson_invert_difla, inv_arg, res_arg,
			       qdpout, qdpin, cgp, cgr, cgsub);
    } else {
      printf0("Cgtype must be 1\n") break;
    };
    

    dtime += QOP_time();

    // get final residual
    //QDP_r_eq_norm2_D(&rsq, qdpout, QDP_all);
    //printf("nrm = %g\n", rsq);

    
    QOP_wilson_ifla_dslash_qdp(NULL, flw, kappa, 1, &tmc, cgr, qdpout, ineo, QOP_EVENODD);
    
    QDP_D_meq_D(cgr, in->df, insub);
    QDP_r_eq_norm2_D(&rsq, cgr, insub);
    
    
    if(res_arg->relmin > 0)
      relnorm2 = QOP_relnorm2_D(&cgr, &qdpout, insub, 1);
    //printf("RELNORM2 = %g\n",relnorm2);
    //printf("%i %i rsq = %g\tprec rsq = %g\trsqstop = %g\n", nrestart,
    //	   res_arg->final_iter, rsq, res_arg->final_rsq, rsqstop);
    res_arg->rsqmin = 0.9*res_arg->final_rsq*rsqstop/rsq;
    /* Do the same for the relative minimum if we are using it */
    if(res_arg->relmin > 0)
      res_arg->relmin = 0.5*res_arg->final_rel/relnorm2 * relminold;
    iter += res_arg->final_iter;
    VERB(LOW, "WILSON_INVERT: iter %i rsq = %g rel = %g\n", iter, rsq, 
	 relnorm2);
  //  } while(((rsqstop > 0 && rsq > rsqstop) &&
  //  (res_arg->relmin > 0 && relnorm2 > res_arg->relmin)) &&
  //  (nrestart++ < max_restarts));
  } while( rsq > rsqstop &&
	   relnorm2 > relminold &&
	   nrestart++ < max_restarts );

  /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

  QDP_D_eq_D(out->df, qdpout, insub);

  QDP_destroy_D(qdpin);
  QDP_destroy_D(qdpout);
  if(flw->clov!=NULL) QDP_destroy_D(ctmp);
  if(QOP_wilson_cgtype==0) QDP_destroy_D(gtmp);

  /* -------------------------------------- */
  inv_arg->max_iter      = max_iter_old;
  inv_arg->max_restarts  = max_restarts_old;
  /* -------------------------------------- */
  res_arg->rsqmin        = rsqminold;
  res_arg->relmin        = relminold;
  res_arg->final_iter    = iter;
  res_arg->final_rsq     = rsq/insq;
  res_arg->final_rel     = relnorm2;
  res_arg->final_restart = nrestart;
  /* -------------------------------------- */
  info->final_sec  = dtime;
  info->final_flop = nflop*res_arg->final_iter*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
  /* -------------------------------------- */

  WILSON_INVERT_END;
}

