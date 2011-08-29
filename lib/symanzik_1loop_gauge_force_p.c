/****** symanzik_1loop_gauge_force_p.c   ************/
/* Gauge action Symanzik improved 1x1 + 1x2 + 1x1x1 */
/* 9/2006 started by L.L. */

#include <qop_internal.h>
#define OPP_DIR(i) (7-i)

#define CHKSUM

static QDP_ColorMatrix *fblink[8];

#if 0
void Print(QDP_ColorMatrix * field){
  QLA_ColorMatrix *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);
  
  mom0 = QDP_expose_M(field);        
  printf("Field: %e %e %e\n\n", mom0[ind].e[0][0].real, mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);

  QDP_reset_M(field);
}
#endif

void 
QOPPC(symanzik_1loop_gauge_force1) (QOP_info_t *info, QOP_GaugeField *gauge, 
		   QOP_Force *force, QOP_gauge_coeffs_t *coeffs, REAL eps)
{
  REAL Plaq, Rect, Pgm ;
  QDP_ColorMatrix *tempmom_qdp[4];
  QDP_ColorMatrix *Amu[6]; // products of 2 links Unu(x)*Umu(x+nu)
  QDP_ColorMatrix *tmpmat;
  QDP_ColorMatrix *tmpmat1;
  QDP_ColorMatrix *tmpmat2;
  QDP_ColorMatrix *staples;
  QDP_ColorMatrix *tmpmat3;
  QDP_ColorMatrix *tmpmat4;

  int i, k;
  int mu, nu, sig;
  double dtime;
  //REAL eb3 = -eps*beta/3.0;
  REAL eb3 = -eps/3.0;
  int j[3][2] = {{1,2},
                 {0,2},
                 {0,1}};
  
  //  QOP_printf0("beta: %e, eb3: %e\n", beta, eb3);
  dtime = -QOP_time();

  for(mu=0; mu<4; mu++) {
    tempmom_qdp[mu] = QDP_create_M();
    QDP_M_eq_zero(tempmom_qdp[mu], QDP_all);
  }

  tmpmat = QDP_create_M();
  for(i=0; i<QOP_common.ndim; i++) {
    fblink[i] = gauge->links[i];
    fblink[OPP_DIR(i)] = QDP_create_M();
    QDP_M_eq_sM(tmpmat, fblink[i], QDP_neighbor[i], QDP_backward, QDP_all);
    QDP_M_eq_Ma(fblink[OPP_DIR(i)], tmpmat, QDP_all);
  }
  

  for(i=0; i<6; i++) {
    Amu[i] = QDP_create_M();
  }

  staples = QDP_create_M();
  tmpmat1 = QDP_create_M();
  tmpmat2 = QDP_create_M();
  tmpmat3 = QDP_create_M();
  tmpmat4 = QDP_create_M();

  Plaq = coeffs->plaquette;
  Rect = coeffs->rectangle;
  Pgm  = coeffs->parallelogram;

  //Construct 3-staples and rectangles
  for(mu=0; mu<4; mu++) {
    i=0;
    for(nu=0; nu<4; nu++) {
      if(nu!=mu){
	// tmpmat1 = Umu(x+nu)
	QDP_M_eq_sM(tmpmat1, fblink[mu], QDP_neighbor[nu], QDP_forward, QDP_all); 
        QDP_M_eq_M_times_M(Amu[i], fblink[nu], tmpmat1, QDP_all);

        //tmpmat2 = Umu(x-nu)
	QDP_M_eq_sM(tmpmat2, fblink[mu], QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_eq_M_times_M(Amu[i+3], fblink[OPP_DIR(nu)], tmpmat2, QDP_all);
       

 
	//tmpmat = U_{nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(staples, Amu[i], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Plaq, staples, QDP_all);
 
        //tmpmat = U_{-nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_Ma_times_M(tmpmat3, fblink[OPP_DIR(nu)], staples, QDP_all);
        QDP_M_eq_M_times_M(tmpmat4, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Rect, tmpmat, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat4, tmpmat2, tmpmat3, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[mu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat3, QDP_all);

        //tmpmat = U_{-nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat3, tmpmat2, tmpmat, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat, tmpmat3, staples, QDP_all);        
        QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat3, QDP_all);




        //tmpmat = U_{-nu}(x+mu) 
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(staples, Amu[i+3], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Plaq, staples, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat3, fblink[nu], staples, QDP_all);
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_M(tmpmat4, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat4, QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[mu], &Rect, tmpmat, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat, tmpmat3, tmpmat1, QDP_all);
        QDP_M_eq_sM(tmpmat4, tmpmat, QDP_neighbor[mu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat4, QDP_all);

        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_M(tmpmat3, staples, tmpmat, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat4, tmpmat3, tmpmat1, QDP_all);
        QDP_M_peq_r_times_M(tempmom_qdp[nu], &Rect, tmpmat4, QDP_all);
        i++;
      }
      
    }

    // Construct the  pgm staples and add them to force
    QDP_M_eq_zero(staples, QDP_all);
    i=0;
    for(nu=0; nu<4; nu++){
      if(nu!=mu){
        k=0;
	for(sig=0; sig<4;sig ++){
	  if(sig!=mu && nu!=sig){
	    
	    // the nu_sig_mu ... staple and 3 reflections
            //tmpmat = Amu["sig"](x+nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]], QDP_neighbor[nu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[nu], tmpmat, QDP_all);   
            //tmpmat3 = Unu(x+mu+sig)
            QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_forward, QDP_all); // HERE?
            //tmpmat2 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = Usig(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[sig], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);


            //tmpmat = Amu["sig"](x-nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]], QDP_neighbor[nu], QDP_backward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["sig"](x-nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(nu)], tmpmat, QDP_all);   
            //tmpmat3 = U_{-nu}(x+mu+sig)
            QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_forward, QDP_all); // HERE?
            //tmpmat2 = U_{-nu}nu(x)*Amu["sig"](x-nu)*adj(Unu(x+mu+sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = Usig(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[sig], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["sig"](x-nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);


            //tmpmat = Amu["-sig"](x-nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]+3], QDP_neighbor[nu], QDP_backward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["-sig"](x-nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(nu)], tmpmat, QDP_all);   
            //tmpmat = U_{-nu}(x+mu-sig)
            QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_backward, QDP_all); // HERE?
            //tmpmat2 = U_{-nu}nu(x)*Amu["-sig"](x-nu)*adj(Unu(x+mu-sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = U_{-sig}(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(sig)], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = U_{-nu}(x)*Amu["-sig"](x-nu)*adj(Unu(x+mu-sig))*adj(U_{-sig}(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);

            


            //tmpmat = Amu["-sig"](x+nu)
	    QDP_M_eq_sM(tmpmat, Amu[j[i][k]+3], QDP_neighbor[nu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["-sig"](x+nu)
            QDP_M_eq_M_times_M(tmpmat1, fblink[nu], tmpmat, QDP_all);   
            //tmpmat3 = Unu(x+mu-sig)
            QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
	    QDP_M_eq_sM(tmpmat3, tmpmat, QDP_neighbor[sig], QDP_backward, QDP_all); // HERE?
            //tmpmat2 = Unu(x)*Amu["-sig"](x+nu)*adj(Unu(x+mu-sig))
	    QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all);
            //tmpmat = U_{-sig}(x+mu)
	    QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(sig)], QDP_neighbor[mu], QDP_forward, QDP_all);
            //tmpmat1 = Unu(x)*Amu["sig"](x+nu)*adj(Unu(x+mu+sig))*adj(Usig(x+mu))
	    QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	    QDP_M_peq_M(staples, tmpmat1, QDP_all);

	    k++;
	  }//close if sig!=nu ...
	}//close sig loop
	i++;
      }// close if nu!=mu
    }//close the pgm nu loop

    QDP_M_peq_r_times_M(tempmom_qdp[mu], &Pgm, staples, QDP_all);
   

    
  }// closes the mu loop

#ifdef CHKSUM
  QLA_ColorMatrix qcm;
  QLA_Complex det, chk;
  QLA_c_eq_r(chk, 0);
#endif
  for(mu=0; mu<4; mu++){
    QDP_M_eq_M_times_Ma(tmpmat, fblink[mu], tempmom_qdp[mu], QDP_all); // HERE?
    QDP_M_eq_r_times_M_plus_M( tempmom_qdp[mu], &eb3, tmpmat, force->force[mu], QDP_all);// HERE?
    QDP_M_eq_antiherm_M(force->force[mu], tempmom_qdp[mu], QDP_all);// HERE
#ifdef CHKSUM
    QDP_m_eq_sum_M(&qcm, force->force[mu], QDP_all);
    QLA_C_eq_det_M(&det, &qcm);
    QLA_c_peq_c(chk, det);
#endif
  }
#ifdef CHKSUM
  QOP_printf0("chksum: %g %g\n", QLA_real(chk), QLA_imag(chk));
#endif

  //DESTROY various fields

  QDP_destroy_M(tmpmat);
  QDP_destroy_M(tmpmat1);
  QDP_destroy_M(tmpmat2);
  QDP_destroy_M(tmpmat3);
  QDP_destroy_M(staples);
  QDP_destroy_M(tmpmat4);

  for(mu=0; mu<4; mu++){
    QDP_destroy_M(tempmom_qdp[mu]);
  }
  for(i=0; i<6; i++) {
    QDP_destroy_M(Amu[i]);
  }

  for(i=4; i<8; i++) {
    QDP_destroy_M(fblink[i]);
  }

  dtime += QOP_time();

  double nflop = 96720;
  info->final_sec = dtime;
  info->final_flop = nflop*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
  //QOP_printf0("Time in slow g_force: %e\n", info->final_sec);
} 

#define QOP_F3_symanzik_1loop_gauge_force QOP_F3_symanzik_1loop_gauge_force2
#define QOP_D3_symanzik_1loop_gauge_force QOP_D3_symanzik_1loop_gauge_force2
#include "symanzik_1loop_gauge_force2_p.c"
#undef QOP_F3_symanzik_1loop_gauge_force
#undef QOP_D3_symanzik_1loop_gauge_force

void 
QOPPC(symanzik_1loop_gauge_force)(QOP_info_t *info, QOP_GaugeField *gauge, 
				  QOP_Force *force, QOP_gauge_coeffs_t *coeffs, REAL eps)
{
  static int n=0;
  if(QOP_common.verbosity==QOP_VERB_DEBUG) {
    QOPPC(symanzik_1loop_gauge_force2)(info, gauge, force, coeffs, eps);
  } else {
    if(((n/6)&1)==0) {
      QOPPC(symanzik_1loop_gauge_force1)(info, gauge, force, coeffs, eps);
    } else {
      QOPPC(symanzik_1loop_gauge_force2)(info, gauge, force, coeffs, eps);
    }
    n++;
  }
}
