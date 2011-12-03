/******* hisq_links_p.c ****************/


#include <qop_internal.h>

/* 10/1/2009 A.Bazavov, modified version from Alan Gray,
                        introduced function to allocate FermionLinksHisq */
/* 4/23/2011 C. DeTar, simplified the HISQ links structure and API */


#ifdef AB_DEBUG_QOP_HISQ

/* Deugging tools */

static void 
PrintM(QDP_ColorMatrix * field){
  QLA_ColorMatrix *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);

  mom0 = QDP_expose_M(field);
  printf("Field: %e %e %e\n", mom0[ind].e[0][0].real,
	 mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);
  printf("       %e %e %e\n", mom0[ind].e[1][0].real,
	 mom0[ind].e[1][1].real, mom0[ind].e[1][2].real);
  printf("       %e %e %e\n\n", mom0[ind].e[2][0].real,
	 mom0[ind].e[2][1].real, mom0[ind].e[2][2].real);

  QDP_reset_M(field);
}

static void 
PrintMsite(QDP_ColorMatrix * field, int *x){
  QLA_ColorMatrix *mom0;
  int ind = QDP_index(x);

  mom0 = QDP_expose_M(field);
  printf("Field: %e %e %e\n", mom0[ind].e[0][0].real,
	 mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);
  printf("       %e %e %e\n", mom0[ind].e[1][0].real,
	 mom0[ind].e[1][1].real, mom0[ind].e[1][2].real);
  printf("       %e %e %e\n\n", mom0[ind].e[2][0].real,
	 mom0[ind].e[2][1].real, mom0[ind].e[2][2].real);

  QDP_reset_M(field);
}

static void 
PrintMone(QLA_ColorMatrix m){
  printf("Field: %e %e %e\n", m.e[0][0].real,
	 m.e[0][1].real, m.e[0][2].real);
  printf("       %e %e %e\n", m.e[1][0].real,
	 m.e[1][1].real, m.e[1][2].real);
  printf("       %e %e %e\n\n", m.e[2][0].real,
	 m.e[2][1].real, m.e[2][2].real);
}

#endif


/* CD: This is now internal to QOP */
static QOPPC(FermionLinksHisq)*
allocate_hisq_fermion_links(QOP_hisq_coeffs_t *coeffs)
{

  int n_naiks = coeffs->n_naiks;
  int i;
  QOP_hisq_unitarize_group_t ugroup = coeffs->ugroup;
  QOPPC(FermionLinksHisq)* flh;

#if (QOP_Precision==1)
//ABprint  printf("Enter QOP_F3_allocate_hisq_fermion_links\n");
#else
//ABprint  printf("Enter QOP_D3_allocate_hisq_fermion_links\n");
#endif

  //AB Allocate space for the structure itself
  QOP_malloc( flh, QOPPC(FermionLinksHisq), 1 );

  if(ugroup == QOP_UNITARIZE_SU3)
    flh->WeqY = 0;
  else
    flh->WeqY = 1;

  //AB Allocate space for all intermediate links
  QOP_malloc(flh->U_links, QDPPC(ColorMatrix) *, 4); // gauge links
  QOP_malloc(flh->V_links, QDPPC(ColorMatrix) *, 4); // Fat7 smeared
  QOP_malloc(flh->Y_unitlinks, QDPPC(ColorMatrix) *, 4); // projected U(3)
  QOP_malloc(flh->W_unitlinks, QDPPC(ColorMatrix) *, 4); // projected SU(3)
  for(i=0;i<4;i++) {
    flh->U_links[i]=QDP_create_M();
    flh->V_links[i]=QDP_create_M();
    flh->Y_unitlinks[i]=QDP_create_M();
    if(flh->WeqY)
      flh->W_unitlinks[i]=flh->Y_unitlinks[i];
    else
      flh->W_unitlinks[i]=QDP_create_M();
  }

  flh->n_naiks = n_naiks;
  QOP_malloc(flh->fn, QOPPC(FermionLinksAsqtad) *, n_naiks); // Asqtad links

#if (QOP_Precision==1)
//ABprint  printf("Exit  QOP_F3_allocate_hisq_fermion_links\n");
#else
//ABprint  printf("Exit  QOP_D3_allocate_hisq_fermion_links\n");
#endif

  return flh;

}

static void
destroy_4M( QDPPC(ColorMatrix) **links ){
  int i;

  for(i = 0; i < 4; i++)
    if(links[i] != NULL)
      QDP_destroy_M(links[i]);
  
  free(links);
  
}

//create HISQ links from QOPPC(GaugeField)
//CD Here we follow closely what the simplified MILC code does

QOPPC(FermionLinksHisq) *
QOPPC(hisq_create_L_from_G)(QOP_info_t *info,
			    QOP_hisq_coeffs_t *coeffs,
			    QOPPC(GaugeField) *gauge)
{
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter QOPPC(hisq_create_L_from_G) in hisq_links_p.c\n");
#endif

  int n_naiks = coeffs->n_naiks;
  double *eps_naik = coeffs->eps_naik;
  QOP_hisq_unitarize_group_t ugroup = coeffs->ugroup;
  int want_aux = QOP_hisq_links.want_aux;
  int want_deps = QOP_hisq_links.want_deps;
  QOPPC(FermionLinksHisq) *flh;
  QOPPC(FermionLinksAsqtad) *fla;
  QOP_asqtad_coeffs_t acoeffs1, acoeffs2, acoeffs3;
  int i,inaik;
  QOPPC(GaugeField) *qopgf_tmp;
  double final_flop = 0.0;
  double dtime = -QOP_time();

  HISQ_LINKS_BEGIN;

  flh = allocate_hisq_fermion_links(coeffs);

  //AB Copy the gauge field in FermionLinksHisq structure
  if(want_aux)
    QOPPC(extract_G_to_qdp)(&(flh->U_links[0]), gauge);

#ifdef AB_DEBUG_QOP_HISQ
  printf("*** U check in hisq_links_p.c\n");
  PrintM(flh->U_links[0]);fflush(stdout);
#endif


  // fat7 stage 
  // convert hisq style coeffs to asqtad style 
  acoeffs1.one_link     = coeffs->fat7_one_link;
  acoeffs1.three_staple = coeffs->fat7_three_staple;
  acoeffs1.five_staple  = coeffs->fat7_five_staple;
  acoeffs1.seven_staple = coeffs->fat7_seven_staple;
  // by definition lepage and naik coeffs are 0 for fat7
  acoeffs1.lepage = 0.;
  acoeffs1.naik = 0.;

  acoeffs2.one_link     = coeffs->asqtad_one_link;
  acoeffs2.three_staple = coeffs->asqtad_three_staple;
  acoeffs2.five_staple  = coeffs->asqtad_five_staple;
  acoeffs2.seven_staple = coeffs->asqtad_seven_staple;
  acoeffs2.lepage       = coeffs->asqtad_lepage;
  acoeffs2.naik         = coeffs->asqtad_naik;

  acoeffs3.one_link     = coeffs->difference_one_link;
  acoeffs3.three_staple = 0;
  acoeffs3.five_staple  = 0;
  acoeffs3.seven_staple = 0;
  acoeffs3.lepage       = 0;
  acoeffs3.naik         = coeffs->difference_naik;
  
  // smear gauge with fat7 
  fla = QOP_asqtad_create_L_from_G(info, &acoeffs1, gauge);
  final_flop += info->final_flop;
  
  // convert fla back to gauge to enable further smearing 
  // since fat7 contains only fat links, can just copy fla->fatlinks to gauge
  for (i=0;i<4;i++)
    {
      QLA_Real two = 2.0; // correct for factor of 1/2 in links
      QDP_M_eq_r_times_M(flh->V_links[i], &two, fla->fatlinks[i], QDP_all);
    }
  QOP_asqtad_destroy_L(fla);
  final_flop += 198 * QDP_sites_on_node;



  //AB CAREFUL HERE: MAY NEED TO REMOVE PHASES, UNITARIZE
  //   AND THEN PUT PHASES BACK IN
  
  // projection stage 
  // project gauge fields to SU3 or U3

  if(ugroup == QOP_UNITARIZE_U3)
    for (i=0;i<4;i++){
      QOPPC(u3reunit)(info, flh->V_links[i],flh->Y_unitlinks[i]);
      final_flop += info->final_flop;
    }


  if(!want_aux){
    destroy_4M(flh->V_links);
    flh->V_links = NULL;
  }

#ifdef AB_DEBUG_QOP_HISQ
  printf("******* matrix after reunit in hisq_links_p.c\n");fflush(stdout);
  PrintM(flh->Y_unitlinks[0]);fflush(stdout);
#endif

  // SU(3) projection, if requested

  if(ugroup == QOP_UNITARIZE_SU3){
    // TO DO: WE NEED TO REMOVE THE STANDARD KS PHASES FIRST
    // THEN THIS WILL WORK.
QOP_printf0("QOP_hisq_create_L_from_G: SU(3) projection not supported for now\n");
exit(1);
    for (i=0;i<4;i++)
      QOPPC(su3reunit)(info, flh->Y_unitlinks[i],flh->W_unitlinks[i]);
  }

  if(!want_aux && !flh->WeqY){
    destroy_4M(flh->Y_unitlinks);
    flh->Y_unitlinks = NULL;
  }


  // prepare the unitarized links
  qopgf_tmp = QOP_create_G_from_qdp( &(flh->W_unitlinks[0]) );

  if(!want_aux){
    destroy_4M(flh->W_unitlinks);
    flh->W_unitlinks = NULL;
  }

  flh->fn_deps = NULL;

  if(n_naiks > 1) { // Have non-zero epsilon corrections
    // so calculate the difference with 3rd path table set first
    // and then proceed with 2nd path table (Asqtad)
    
    // convert coefficients
    //AB THIS NEEDS TO BE CHANGED SINCE WE DO NOT WANT TO DO
    //   HUGE NUMBER MULTIPLICATIONS BY 0
    
    // 3rd path table set
    flh->fn[0] = QOP_asqtad_create_L_from_G(info, &acoeffs3, qopgf_tmp);
    final_flop += info->final_flop;

    // Copy to fn_deps if we want the derivative wrto epsilon
    if(want_deps)
      flh->fn_deps = QOPPC(asqtad_create_L_from_L)(flh->fn[0]);


    for( inaik = 1; inaik < n_naiks; inaik++ ) {
      QLA_Real rr;
      rr = eps_naik[inaik];
      flh->fn[inaik] = QOPPC(asqtad_create_L_from_r_times_L)(&rr, flh->fn[0]);
      if(flh->fn[inaik]->longlinks)
	final_flop += 8*18*QDP_sites_on_node*(n_naiks-1);
      else
	final_flop += 4*18*QDP_sites_on_node*(n_naiks-1);
    }
    

    // dispose Asqtad links object
    QOP_asqtad_destroy_L(flh->fn[0]);
    
    // asqtad stage 

    // 2nd path table set 
    flh->fn[0] = QOP_asqtad_create_L_from_G(info, &acoeffs2, qopgf_tmp);
    final_flop += info->final_flop;

    for( inaik = 1; inaik < n_naiks; inaik++) {
      QOPPC(asqtad_L_peq_L)(flh->fn[inaik], flh->fn[0]);

      if(flh->fn[inaik]->longlinks)
	final_flop += 8*18*QDP_sites_on_node;
      else
	final_flop += 4*18*QDP_sites_on_node;
    }

  } else { // there are no Naik corrections, n_naiks=1 (can't be 0)

    // asqtad stage 

    flh->fn[0] = QOP_asqtad_create_L_from_G(info, &acoeffs2, qopgf_tmp);
    final_flop += info->final_flop;
    
    if(want_deps){
      flh->fn_deps = QOP_asqtad_create_L_from_G(info, &acoeffs3, qopgf_tmp);
      final_flop += info->final_flop;
    }

  }

  QOP_destroy_G(qopgf_tmp);

#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  QOPPC(hisq_create_L_from_G) in hisq_links_p.c\n");
#endif

  dtime += QOP_time();
  info->final_flop = final_flop;
  info->final_sec = dtime;

  HISQ_LINKS_END;

  return flh;
}

void
QOPPC(hisq_destroy_L)(QOPPC(FermionLinksHisq) *flh)
{
  int i;
  int n_naiks = flh->n_naiks;

  HISQ_LINKS_BEGIN;

  for(i=0; i < n_naiks; i++)
    QOP_asqtad_destroy_L(flh->fn[i]);
  free(flh->fn);

  QOP_asqtad_destroy_L(flh->fn_deps);

  if(flh->W_unitlinks != NULL)
    destroy_4M(flh->W_unitlinks);
  if(!flh->WeqY){
    if(flh->Y_unitlinks != NULL)
      destroy_4M(flh->Y_unitlinks);
  }
  else
    free(flh->Y_unitlinks);
  if(flh->V_links != NULL)
    destroy_4M(flh->V_links);
  if(flh->U_links != NULL)
    destroy_4M(flh->U_links);

  free(flh);

  HISQ_LINKS_END;

}

QOPPC(FermionLinksAsqtad) **
  QOPPC(get_asqtad_links_from_hisq)(QOPPC(FermionLinksHisq) *flh){

  return flh->fn;
}

QOPPC(FermionLinksAsqtad) *
  QOPPC(get_asqtad_deps_links_from_hisq)(QOPPC(FermionLinksHisq) *flh){

  return flh->fn_deps;
}

