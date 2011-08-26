/**************** hisq_force_fnmat2_p.c ******************************/
// HISQ force using new optimized asqtad routines

#include <qop_internal.h>

extern void 
QOPPC(get_mid)(QOP_info_t *info, QDP_ColorMatrix *mid[], QDP_Shift shifts[], int ns,
	       REAL eps[], QOP_ColorVector *x[], int nterms);

extern void
QOPPC(asqtad_deriv)(QOP_info_t *info, QDP_ColorMatrix *gauge[],
		    QDP_ColorMatrix *force[], QOP_asqtad_coeffs_t *coef,
		    QDP_ColorMatrix *mid_fat[], QDP_ColorMatrix *mid_naik[]);

void 
QOPPC(hisq_force_multi_wrapper_fnmat2)(QOP_info_t *info,  
				       QOPPC(FermionLinksHisq) *flh,
				       QOP_Force *Force, 
				       QOP_hisq_coeffs_t *hisq_coeff,
				       REAL *residues,
				       QOP_ColorVector *x[], 
				       int *n_orders_naik)
{
  if(!QOP_asqtad.inited) QOP_asqtad_invert_init();

  double dtime = QDP_time();
  int nflops = 0;
  QOP_info_t tinfo;

  QDP_ColorMatrix *force[4] =  {Force->force[0], Force->force[1], 
				Force->force[2], Force->force[3]};

  QDP_ColorMatrix *Ugf[4], *Vgf[4], *Wgf[4];
  for(int i=0; i<4; i++) {
    Ugf[i] = flh->U_links[i];
    Vgf[i] = flh->V_links[i];
    Wgf[i] = flh->W_unitlinks[i];
  }

  QDP_ColorMatrix *force_accum_0[4];
  QDP_ColorMatrix *force_accum_0_naik[4];
  QDP_ColorMatrix *force_accum_1[4];
  QDP_ColorMatrix *force_accum_1u[4];
  QDP_ColorMatrix *force_accum_2[4];
  QDP_ColorMatrix *force_final[4];
  QDP_ColorMatrix *tmat = QDP_create_M();
  for(int i=0; i<4; i++) {
     force_accum_0[i] = QDP_create_M();
     force_accum_0_naik[i] = QDP_create_M();
     force_accum_1[i] = QDP_create_M();
     force_accum_1u[i] = QDP_create_M();
     force_accum_2[i] = QDP_create_M();
     force_final[i] = QDP_create_M();
     QDP_M_eq_zero(force_accum_2[i], QDP_all);
  }

  int n_naiks = hisq_coeff->n_naiks;
  int nterms = 0;
  for(int inaik = 0; inaik < n_naiks; inaik++)
    nterms += n_orders_naik[inaik];

  // loop on different naik masses
  int n_naik_shift = 0;
  for(int inaik=0; inaik<n_naiks; inaik++) {
    int n_orders_naik_current;
    if( inaik==0 ) {
      n_orders_naik_current = nterms;
    } else {
      n_orders_naik_current = n_orders_naik[inaik];
    }

    QOPPC(get_mid)(&tinfo, force_accum_0, QDP_neighbor, 4,
		   residues+n_naik_shift, x+n_naik_shift, n_orders_naik_current);
    nflops += tinfo.final_flop;
    QOPPC(get_mid)(&tinfo, force_accum_0_naik, QOP_common.neighbor3, 4,
		   residues+n_naik_shift, x+n_naik_shift, n_orders_naik_current);
    nflops += tinfo.final_flop;

    // smearing level 0
    for(int i=0; i<4; i++)
      QDP_M_eq_zero(force_accum_1[i], QDP_all);
    if(inaik==0) {
      QOP_asqtad_coeffs_t acoef;
      acoef.one_link = hisq_coeff->asqtad_one_link;
      acoef.three_staple = hisq_coeff->asqtad_three_staple;
      acoef.five_staple = hisq_coeff->asqtad_five_staple;
      acoef.seven_staple = hisq_coeff->asqtad_seven_staple;
      acoef.lepage = hisq_coeff->asqtad_lepage;
      acoef.naik = hisq_coeff->asqtad_naik;
      QOPPC(asqtad_deriv)(&tinfo, Wgf, force_accum_1, &acoef,
			  force_accum_0, force_accum_0_naik);
      nflops += tinfo.final_flop;
    } else {
      QOP_asqtad_coeffs_t acoef;
      acoef.one_link = hisq_coeff->difference_one_link;
      acoef.three_staple = 0;
      acoef.five_staple = 0;
      acoef.seven_staple = 0;
      acoef.lepage = 0;
      acoef.naik = hisq_coeff->difference_naik;
      QOPPC(asqtad_deriv)(&tinfo, Wgf, force_accum_1, &acoef,
			  force_accum_0, force_accum_0_naik);
      nflops += tinfo.final_flop;
    }

    QLA_Real coeff_mult;
    if( inaik==0 ) {
      coeff_mult = 1.0;
    } else {
      coeff_mult = hisq_coeff->eps_naik[inaik];
    }
    for(int dir=0; dir<4; dir++) {
      QDP_M_peq_r_times_M(force_accum_2[dir], &coeff_mult,
			  force_accum_1[dir], QDP_all);
    }
    nflops += 4*36*QDP_sites_on_node;

    n_naik_shift += n_orders_naik[inaik];
  }

  QOP_hisq_unitarize_method_t umethod = hisq_coeff->umethod;
  if ( umethod==QOP_UNITARIZE_NONE ){

    // smearing level 1
    QOP_asqtad_coeffs_t acoef;
    acoef.one_link = hisq_coeff->fat7_one_link;
    acoef.three_staple = hisq_coeff->fat7_three_staple;
    acoef.five_staple = hisq_coeff->fat7_five_staple;
    acoef.seven_staple = hisq_coeff->fat7_seven_staple;
    acoef.lepage = 0;
    acoef.naik = 0;
    for(int dir=0; dir<4; dir++)
      QDP_M_eq_zero(force_accum_1[dir], QDP_all);
    QOPPC(asqtad_deriv)(&tinfo, Ugf, force_accum_1, &acoef,
			force_accum_2, NULL);
    nflops += tinfo.final_flop;

  } else if ( umethod==QOP_UNITARIZE_RATIONAL ) {

    for(int mu=0; mu<4; mu++) QDP_M_eq_Ma(force_accum_1u[mu], force_accum_2[mu], QDP_all);
    // reunitarization
    QOPPC(hisq_force_multi_reunit)(info, Vgf, force_accum_2, force_accum_1u);
    for(int mu=0; mu<4; mu++) QDP_M_eq_Ma(force_accum_1u[mu], force_accum_2[mu], QDP_all);
    nflops += 1000*QDP_sites_on_node; // ??? need to count

    // smearing level 1
    QOP_asqtad_coeffs_t acoef;
    acoef.one_link = hisq_coeff->fat7_one_link;
    acoef.three_staple = hisq_coeff->fat7_three_staple;
    acoef.five_staple = hisq_coeff->fat7_five_staple;
    acoef.seven_staple = hisq_coeff->fat7_seven_staple;
    acoef.lepage = 0;
    acoef.naik = 0;
    for(int dir=XUP;dir<=TUP;dir++)
      QDP_M_eq_zero(force_accum_1[dir], QDP_all);
    QOPPC(asqtad_deriv)(&tinfo, Ugf, force_accum_1, &acoef,
			force_accum_1u, NULL);
    nflops += tinfo.final_flop;

  } else {
    QOP_printf0("Unknown or unsupported unitarization method\n");
    exit(1);
  }

  // contraction with the link in question should be done here,
  // after contributions from all levels of smearing are taken into account
  for(int dir=0; dir<4; dir++) {
    QDP_M_eq_M_times_Ma(force_final[dir], Ugf[dir], force_accum_1[dir], QDP_all);
  }
  nflops += 4*198*QDP_sites_on_node;

  // take into account even/odd parity (it is NOT done in "smearing" routine)
  //eps multiplication done outside QOP 
  for(int dir=0; dir<4; dir++) {
    QLA_Real treal = 2;
    QDP_M_eq_M(tmat, force_final[dir], QDP_all);
    QDP_M_eq_r_times_M(force_final[dir], &treal, tmat, QDP_even);
    QDP_M_eqm_r_times_M(force_final[dir], &treal, tmat, QDP_odd);
  }
  nflops += 4*18*QDP_sites_on_node;

  // Put antihermitian traceless part into momentum 
  // add force to momentum
  for(int dir=0; dir<4; dir++) {
    QDP_M_eq_antiherm_M(tmat, force_final[dir], QDP_all);
    QDP_M_peq_M(force[dir], tmat, QDP_all);
  }
  nflops += 4*(24+18)*QDP_sites_on_node;

  for(int i=0; i<4; i++) {
     QDP_destroy_M( force_accum_0[i] );
     QDP_destroy_M( force_accum_0_naik[i] );
     QDP_destroy_M( force_accum_1[i] );
     QDP_destroy_M( force_accum_1u[i] );
     QDP_destroy_M( force_accum_2[i] );
     QDP_destroy_M( force_final[i] );
  }
  QDP_destroy_M( tmat );

  info->final_sec = QDP_time() - dtime;
  info->final_flop = nflops;
  info->status = QOP_SUCCESS;
}
