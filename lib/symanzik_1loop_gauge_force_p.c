/****** symanzik_1loop__force_p.c   ************/
/* Gauge action Symanzik improved 1x1 + 1x2 + 1x1x1 */
/* 9/2006 started by L.L. */

#include <qop_internal.h>
#define OPP_DIR(i) (7-i)


static QDP_ColorMatrix *fblink[8];
extern REAL beta;


void Print(QDP_ColorMatrix * field){
  QLA_ColorMatrix *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);
  
  mom0 = QDP_expose_M(field);        
  printf("Field: %e %e %e\n\n", mom0[ind].e[0][0].real, mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);

  QDP_reset_M(field);
}

void 
QOPPC(symanzik_1loop_gauge_force) (QOP_info_t *info, QOP_GaugeField *gauge, QOP_Force *force,
			QOP_gauge_coeffs_t *coeffs, REAL eps)
{
  REAL Plaq, Rect, Pgm ;
  QDP_ColorMatrix *tempmom_qdp[4];
  QDP_ColorMatrix *Amu[6]; // products of 2 links Unu(x)*Umu(x+nu)
  //QDP_ColorMatrix *Smu[6]; // staples Unu(x)*Umu(x+nu)*Unu^+(x+mu)
  QDP_ColorMatrix *tmpmat;
  QDP_ColorMatrix *tmpmat1;
  QDP_ColorMatrix *tmpmat2;
  QDP_ColorMatrix *staples;
  QDP_ColorMatrix *tmpmat3;
  int i, k;
  int mu, nu, sig;
  double dtime;
  REAL eb3 = -eps*beta/3.0;

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
    //Smu[i] = QDP_create_M();
  }

  staples = QDP_create_M();
  tmpmat1 = QDP_create_M();
  tmpmat2 = QDP_create_M();
  tmpmat3 = QDP_create_M();

  Plaq = coeffs->plaquette;
  Rect = coeffs->rectangle;
  Pgm  = coeffs->parallelogram;

  //QOP_printf0("before main mu loop\n Plaq %e, Rect %e, Pgm %e\n", Plaq, Rect, Pgm);
  for(mu=0; mu<4; mu++) {
    QDP_M_eq_zero(staples, QDP_all);
    i=0;
    // Construct 2 link products to be reused for dir mu
    for(nu=0; nu<4; nu++) {
      if(nu!=mu){
	// tmpmat = Umu(x+nu)
	QDP_M_eq_sM(tmpmat, fblink[mu], QDP_neighbor[nu], QDP_forward, QDP_all); 
        QDP_M_eq_M_times_M(Amu[i], fblink[nu], tmpmat, QDP_all);

        //tmpmat = Umu(x-nu)
	QDP_M_eq_sM(tmpmat, fblink[mu], QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_eq_M_times_M(Amu[i+3], fblink[OPP_DIR(nu)], tmpmat, QDP_all);
        
	//tmpmat = U_{nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat1, Amu[i], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(staples, &Plaq, tmpmat1, QDP_all);
       
        //tmpmat = U_{-nu}(x+mu) 
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_Ma(tmpmat2, Amu[i+3], tmpmat, QDP_all);        
        QDP_M_peq_r_times_M(staples, &Plaq, tmpmat2, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat3, fblink[OPP_DIR(nu)], tmpmat1, QDP_all);
        QDP_M_eq_M_times_M(tmpmat1, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat1, QDP_neighbor[nu], QDP_forward, QDP_all);
        QDP_M_peq_r_times_M(staples, &Rect, tmpmat, QDP_all);

        QDP_M_eq_Ma_times_M(tmpmat3, fblink[nu], tmpmat2, QDP_all);
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_M_times_M(tmpmat1, tmpmat3, tmpmat, QDP_all);
        QDP_M_eq_sM(tmpmat, tmpmat1, QDP_neighbor[nu], QDP_backward, QDP_all);
        QDP_M_peq_r_times_M(staples, &Rect, tmpmat, QDP_all);

        i++;
      }
      
    }
    QDP_M_peq_M(tempmom_qdp[mu], staples, QDP_all); 
    //QOP_printf0("After initializing Amu and Smu at mu %i\n", mu);
    //Add 3 staples to force. Is the force initiallized at beginning?
        
    
    //Construct rect 5 staples   
    QDP_M_eq_zero(staples, QDP_all);
    //Print(staples);
    i=0;
    for(nu=0; nu<4; nu++) {
      if(nu!=mu){
        // Second, the nu_mu_mu_nu_mu staple and reflection
//tmpmat = Umu(x+mu), keep that one for the rest of this loop
	QDP_M_eq_sM(tmpmat, fblink[mu], QDP_neighbor[mu], QDP_forward, QDP_all);
        //tmpmat2 = Umu(x+mu+nu)
	QDP_M_eq_sM(tmpmat2, tmpmat, QDP_neighbor[nu], QDP_forward, QDP_all);	         
//tmpmat1 = Amu[i]*Umu(x+mu+nu)
        QDP_M_eq_M_times_M(tmpmat1, Amu[i], tmpmat2, QDP_all); // HERE3?
	//QOP_printf0("In rect loop. After HERE3\n");
        //tmpmat3 = Unu(x+2mu)
	QDP_M_eq_sM(tmpmat2, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        QDP_M_eq_sM(tmpmat3, tmpmat2, QDP_neighbor[mu], QDP_forward, QDP_all); // HERE4?
        //QOP_printf0("In rect loop. After HERE4\n");
        //tmpmat2 = Amu[i]*Umu(x+mu+nu)*adj(Unu(x+2mu))
        QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat3, QDP_all); // HERE5?
        //QOP_printf0("In rect loop. After HERE5\n");
        //tmpmat1 = Amu[i]*Umu(x+mu+nu)*adj(Unu(x+2mu))*adj(Umu(x+mu))
	QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);  // HERE6?
	//QOP_printf0("In rect loop. After HERE6\n");
	QDP_M_peq_M(staples, tmpmat1, QDP_all);
        //printf("Second term\n");
        //Print(staples);


        //tmpmat1 = Umu(x+mu-nu)
	QDP_M_eq_sM(tmpmat1, tmpmat, QDP_neighbor[nu], QDP_backward, QDP_all);
        //tmpmat2 = Amu[i+3]*Umu(x+mu-nu)
        QDP_M_eq_M_times_M(tmpmat2, Amu[i+3], tmpmat1, QDP_all); // HERE7?
	//QOP_printf0("In rect loop. After HERE7\n");
   	//tmpmat3 = U_{-nu}(x+2mu)
	QDP_M_eq_sM(tmpmat1, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);// HERE7a?
	QDP_M_eq_sM(tmpmat3, tmpmat1,             QDP_neighbor[mu], QDP_forward, QDP_all); // HERE8?
        //tmpmat1 = Amu[i+3]*Umu(x+mu-nu)*adj(U_{-nu}(x+2mu))
        QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat3, QDP_all); // HERE9?
	//QOP_printf0("In rect loop. After HERE9\n");
        //tmpmat2 = Amu[i]*Umu(x+mu+nu)*adj(Unu(x+2mu))*adj(Umu(x+mu))
	QDP_M_eq_M_times_Ma(tmpmat2, tmpmat1, tmpmat, QDP_all); // HERE10?
        //QOP_printf0("In rect loop. After HERE10\n");

	QDP_M_peq_M(staples, tmpmat2, QDP_all);
	//printf("Second term refl.\n");
        //Print(staples);

        //tmpmat = Amu["nu"](x-mu)
        QDP_M_eq_sM(tmpmat, Amu[i], QDP_neighbor[mu], QDP_backward, QDP_all);
	//tmpmat1 = U_{-mu}(x)*Amu["nu"](x-mu)
	QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(mu)], tmpmat, QDP_all);
        //tmpmat = Umu(x+nu)
        QDP_M_eq_sM(tmpmat, fblink[mu], QDP_neighbor[nu], QDP_forward, QDP_all);
        //tmpmat2 = U_{-mu}(x)*Amu["nu"](x-mu)*Umu(x+nu)
	QDP_M_eq_M_times_M(tmpmat2, tmpmat1, tmpmat, QDP_all);
        //tmpmat = Unu(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
        //tmpmat1 = U_{-mu}(x)*Amu["nu"](x-mu)*Umu(x+nu)*adj(Unu(x+mu))
	QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	QDP_M_peq_M(staples, tmpmat1, QDP_all);
	//printf("Third term\n");
        //Print(staples);

       //tmpmat = Amu["-nu"](x-mu)
        QDP_M_eq_sM(tmpmat, Amu[i+3], QDP_neighbor[mu], QDP_backward, QDP_all);
	//tmpmat1 = U_{-mu}(x)*Amu["-nu"](x-mu)
	QDP_M_eq_M_times_M(tmpmat1, fblink[OPP_DIR(mu)], tmpmat, QDP_all);
        //tmpmat = Umu(x-nu)
        QDP_M_eq_sM(tmpmat, fblink[mu], QDP_neighbor[nu], QDP_backward, QDP_all);
        //tmpmat2 = U_{-mu}(x)*Amu["-nu"](x-mu)*Umu(x-nu)
	QDP_M_eq_M_times_M(tmpmat2, tmpmat1, tmpmat, QDP_all);
        //tmpmat = U_{-nu}(x+mu)
        QDP_M_eq_sM(tmpmat, fblink[OPP_DIR(nu)], QDP_neighbor[mu], QDP_forward, QDP_all);
        //tmpmat1 = U_{-mu}(x)*Amu["-nu"](x-mu)*Umu(x-nu)*adj(U_{-nu)x+mu))
	QDP_M_eq_M_times_Ma(tmpmat1, tmpmat2, tmpmat, QDP_all);

	QDP_M_peq_M(staples, tmpmat1, QDP_all);
	//printf("Third term and refl\n");
        //Print(staples);


	i++;
	//QOP_printf0("In rect loop. End of if\n");
      }//if(nu!=mu);
    }
    //QOP_printf0("End rect loop. Before staples addition Rect coef: %e\n", Rect);
    //printf("Rect staples sum in mu = %i:\n", mu); 
    //Print(staples);
    QDP_M_peq_r_times_M(tempmom_qdp[mu], &Rect, staples, QDP_all);
    //QOP_printf0("End rect loop. After staples addition\n");

    // Construct the 5 pgm staples and add them to force
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

  for(mu=0; mu<4; mu++){
    QDP_M_eq_M_times_Ma(tmpmat, fblink[mu], tempmom_qdp[mu], QDP_all); // HERE?
    QDP_M_eq_r_times_M_plus_M( tempmom_qdp[mu], &eb3, tmpmat, force->force[mu], QDP_all);// HERE?
    QDP_M_eq_antiherm_M(force->force[mu], tempmom_qdp[mu], QDP_all);// HERE
   }
  /*
  for(mu=0; mu<4; mu++){
    QDP_M_eq_M_times_Ma(tmpmat, fblink[mu], tempmom_qdp[mu], QDP_all); // HERE?
    QDP_M_eq_antiherm_M(tmpmat1, tmpmat, QDP_all);// HERE?
    QDP_M_peq_r_times_M( force->force[mu], &eb3, tmpmat1, QDP_all);// HERE?
  }
  */


  //DESTROY various fields

  QDP_destroy_M(tmpmat);
  QDP_destroy_M(tmpmat1);
  QDP_destroy_M(tmpmat2);
  QDP_destroy_M(tmpmat3);
  QDP_destroy_M(staples);
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

  int nflop = 112624;
  info->final_sec = dtime;
  info->final_flop = nflop*QDP_sites_on_node; 
  info->status = QOP_SUCCESS;
  //QOP_printf0("Time in slow g_force: %e\n", info->final_sec);
} 



