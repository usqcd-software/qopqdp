/**************** hisq_force_fnmat_p.c ******************************/
/*
 * Alan Gray (EPCC) HISQ force routines based on MILC code
 * AG: HISQ force routines converted from MILC code 27 Feb 2009
 * CD: Support simplified HISQ links structure 23 Apr 2011
 */

#include <qop_internal.h>
#include <string.h>
#include <hisq_action.h>


/**********************************************************************/
/*   Utilities                                                        */
/**********************************************************************/


#if 0
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
PrintVsite(QDP_ColorVector * field, int *x){
  QLA_ColorVector *mom0;
  int ind = QDP_index(x);

  mom0 = QDP_expose_V(field);
  printf("Field: %e %e\n", mom0[ind].c[0].real, mom0[ind].c[0].imag);
  printf("       %e %e\n", mom0[ind].c[1].real, mom0[ind].c[1].imag);
  printf("       %e %e\n", mom0[ind].c[2].real, mom0[ind].c[2].imag);

  QDP_reset_V(field);
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

static void 
PrintVone(QLA_ColorVector V){
  printf("Field: %e %e\n", V.c[0].real,	V.c[0].imag);
  printf("       %e %e\n", V.c[1].real,	V.c[1].imag);
  printf("       %e %e\n", V.c[2].real,	V.c[2].imag);
}

#endif

#ifdef DEBUG_FNMAT

static void 
PrintM(QDP_ColorMatrix * field){
  QLA_ColorMatrix *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);

  mom0 = QDP_expose_M(field);
  printf("Field: %e %e %e\n\n", mom0[ind].e[0][0].real,
	 mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);
  printf("       %e %e %e\n\n", mom0[ind].e[1][0].real,
	 mom0[ind].e[1][1].real, mom0[ind].e[1][2].real);
  printf("       %e %e %e\n\n", mom0[ind].e[2][0].real,
	 mom0[ind].e[2][1].real, mom0[ind].e[2][2].real);

  QDP_reset_M(field);
}

static void 
PrintV(QDP_ColorVector * field){
  QLA_ColorVector *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);

  mom0 = QDP_expose_V(field);
  printf("Field: %e %e %e\n\n", mom0[ind].c[0].real, 
	 mom0[ind].c[1].real, mom0[ind].c[2].real);

  QDP_reset_V(field);
}

#endif


/********************************************************************/
/* Construct path table from the action definition above
   Originally quark_stuff_hisq.c                                    */
/********************************************************************/

typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  REAL coeff;	/* coefficient, including minus sign if backwards */
  REAL forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

static QDP_Shift shiftdirs[8];
static QDP_Shift neighbor3[4];

static int qop_is_path_equal( int *path1, int* path2, int length );

static int 
qop_is_path_equal( int *path1, int* path2, int length ){
   register int i;
   for(i=0;i<length;i++)if(path1[i]!=path2[i])return(0);
   return(1);
}


static REAL act_path_coeff_1[NUM_BASIC_PATHS_1];
static REAL act_path_coeff_2[NUM_BASIC_PATHS_2];
static REAL act_path_coeff_3[NUM_BASIC_PATHS_3];

static Q_path q_paths_1[MAX_NUM],q_paths_2[MAX_NUM],q_paths_3[MAX_NUM];

/* number of paths in dslash */
static int num_q_paths_1,num_q_paths_2,num_q_paths_3; 

static int 
qop_add_basic_path_hisq( Q_path *this_q_paths, int path_table_index, 
			 int *basic_vec, int length, REAL coeff, 
			 int max_paths );


static int 
make_path_table_hisq( char *action_desc, int npaths, int max_paths,
		      int *path_length, REAL *coeff,
		      int paths[][MAX_LENGTH], REAL *act_coeff, 
		      Q_path *this_q_paths, REAL naik_term_mass, 
		      int index_onelink, int index_naik );


/********************************************************************/
/* Make table of paths in action */
/********************************************************************/
//originally make_path_table in quark_stuff_hisq.c
static void 
QOP_make_paths_and_dirs_hisq(QOP_hisq_coeffs_t *coef, 
			     QOP_hisq_unitarize_method_t umethod) {


  int index_naik = -1, index_onelink = -1;
  
#ifdef INDEX_NAIK
  index_naik = INDEX_NAIK;
  index_onelink = INDEX_ONELINK;
#endif


  int i;
  int disp[4] = {0,0,0,0};

  for(i=0; i<4; i++) {
    disp[i] = 3;
    neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
  }

  for(i=0; i<4; ++i) {
    shiftdirs[i] = QDP_neighbor[i];
    shiftdirs[i+4] = neighbor3[i];
  }


  REAL path_coeff_1[NUM_BASIC_PATHS_1] = {coef->fat7_one_link, 
					  coef->fat7_three_staple, 
					  coef->fat7_five_staple,
					  coef->fat7_seven_staple};

  REAL path_coeff_2[NUM_BASIC_PATHS_2] = {coef->asqtad_one_link, 
					  coef->asqtad_naik, 
					  coef->asqtad_three_staple, 
					  coef->asqtad_five_staple,
					  coef->asqtad_seven_staple, 
					  coef->asqtad_lepage};

  REAL path_coeff_3[NUM_BASIC_PATHS_3] = {coef->difference_one_link, 
					  coef->difference_naik};


  QOP_printf0("QOP MAKING PATH TABLES\n");


  num_q_paths_1 = 
    make_path_table_hisq( QUARK_ACTION_DESCRIPTION_1, 
			  quark_action_npaths_1,
			  MAX_NUM_1, path_length_in_1, 
			  path_coeff_1, path_ind_1, 
			  act_path_coeff_1, q_paths_1, 0.0, -1, -1 );
    

  if ( umethod==QOP_UNITARIZE_NONE )
    {
      QOP_printf0("QOP Unitarization method = QOP_UNITARIZE_NONE\n");
    }
  else if ( umethod==QOP_UNITARIZE_RATIONAL )
    {
      QOP_printf0("QOP Unitarization method = QOP_UNITARIZE_RATIONAL\n");
    }
  else
    {
      QOP_printf0("QOP Unknown or unsupported unitarization method\n"); 
      exit(0);
    }
  

  num_q_paths_2 = 
    make_path_table_hisq( QUARK_ACTION_DESCRIPTION_2, quark_action_npaths_2,
			  MAX_NUM_2, 
			  path_length_in_2, path_coeff_2, path_ind_2, 
			  act_path_coeff_2, q_paths_2, 0.0, -1, -1 );
  
  num_q_paths_3 = 
    make_path_table_hisq( QUARK_ACTION_DESCRIPTION_3, quark_action_npaths_3,
			  MAX_NUM_3, 
			  path_length_in_3, path_coeff_3, path_ind_3, 
			  act_path_coeff_3, q_paths_3, 0.0, -1, -1 );
  
}




static int 
make_path_table_hisq( char *action_desc, int npaths, int max_paths,
		      int *path_length, REAL *coeff,
		      int paths[][MAX_LENGTH], REAL *act_coeff, 
		      Q_path *this_q_paths, REAL naik_term_mass, 
		      int index_onelink, int index_naik ) 
{

  int i,j;
  int n_basic_paths; // number of paths before rotation/reflection
  int n_q_paths; // total number of paths in table
#ifdef TADPOLE_IMPROVE
  int k;
#endif

 

  
  QOP_printf0("%s\n",action_desc);
  n_q_paths = 0;
  n_basic_paths = 0;
  if(MAX_LENGTH > MAX_PATH_LENGTH){
    printf("QOP Path length for this action is too long.  Recompile.\n");
    exit(1);
  }
  
  /* add rots. and reflects to table, print out the path coefficients */
  QOP_printf0("QOP path coefficients: npath  path_coeff  multiplicity\n");
  for(j=0;j<npaths;j++) {
    REAL this_coeff;
    this_coeff = coeff[j];
#ifdef TADPOLE_IMPROVE
        for(k=1;k< path_length[j];k++)this_coeff /= u0;
#endif
    act_coeff[j] = this_coeff ;

    i = qop_add_basic_path_hisq( this_q_paths, n_q_paths, paths[j],
			path_length[j], this_coeff, max_paths );
    n_q_paths += i;
    n_basic_paths++;
    QOP_printf0("                    %d      %e     %d\n", j,this_coeff,i);
  }
  return( n_q_paths );
} //make_path_table_hisq()


/********************************************************************/
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */
/********************************************************************/
static int 
qop_add_basic_path_hisq( Q_path *this_q_paths, int path_table_index, 
		int *basic_vec, int length, REAL coeff, int max_paths ) {
  // this_q_paths is array of paths we are building
    // path_table_index is starting index when called
    // basic_vec is list of directions in basic path
    // length is length of basic path
    // coeff is coefficient in action

    int perm[8],pp[8],ir[4];
    int j,path_num;
    int vec[MAX_LENGTH];
    int flag;
         //node0_printf("ADD BASIC PATH %d:  ",path_table_index);
	 //printpath( basic_vec, length );

    path_num = 0;  // number of paths made from this basic path so far
    /* now fill the long table with all rotations and reflections
	of the fundamental path.  The path presented to us is for
        the positive x component of dslash, so if the x coordinate
        is reflected it will appear with a negative sign. */
      /* permutations */
      for(perm[0]=0;perm[0]<4;perm[0]++)
      for(perm[1]=0;perm[1]<4;perm[1]++)
      for(perm[2]=0;perm[2]<4;perm[2]++)
      for(perm[3]=0;perm[3]<4;perm[3]++){
	if(perm[0] != perm[1] && perm[0] != perm[2] 
	  && perm[0] != perm[3] && perm[1] != perm[2]
	  && perm[1] != perm[3] && perm[2] != perm[3] ) {
	  /* reflections*/
	  for(ir[0]=0;ir[0]<2;ir[0]++)
	  for(ir[1]=0;ir[1]<2;ir[1]++)
	  for(ir[2]=0;ir[2]<2;ir[2]++)
	  for(ir[3]=0;ir[3]<2;ir[3]++){
	    for(j=0;j<4;j++){
	      pp[j]=perm[j];

	      if(ir[j] == 1) pp[j]=OPP_DIR(pp[j]);
	      pp[OPP_DIR(j)]=OPP_DIR(pp[j]);
	    }
	    /* create new vector*/
	    for(j=0;j<length;j++) vec[j]=pp[basic_vec[j]];
	    for(j=length;j<MAX_LENGTH;j++) vec[j]=NODIR;

            flag=0;
	    /* check if it's a new set: */
	    for(j=0;j<path_table_index;j++){
	      flag = qop_is_path_equal( vec, this_q_paths[j].dir, 
					MAX_LENGTH );
	      if(flag==1)break;
	    }
	    if(flag == 0 ){
	      if(path_table_index>=max_paths){
		QOP_printf0("OOPS: MAX_NUM too small\n");
		exit(0);
	      }
	      this_q_paths[path_table_index].length=length;
	      for(j=0;j<MAX_LENGTH;j++) {
		this_q_paths[path_table_index].dir[j]=vec[j];
	      }
		/* remember to copy NODIR's, or comparison will get confused */
	      if(ir[0]==0){
		this_q_paths[path_table_index].coeff =  coeff;
		this_q_paths[path_table_index].forwback =  +1;
	      }
	      else{
		this_q_paths[path_table_index].coeff = -coeff;
		this_q_paths[path_table_index].forwback = -1;
	      }
	      path_table_index++;
	      path_num++;
	    }

	  } /* end reflection*/
        } /* end permutation if block */
      } /* end permutation */

    return(path_num);
} /* add_basic_path_hiqs */



static int 
qop_get_num_q_paths_1(){
  return num_q_paths_1;
}

static int 
qop_get_num_q_paths_2(){
  return num_q_paths_2;
}

static int 
qop_get_num_q_paths_3(){
  return num_q_paths_3;
}


static Q_path *
qop_get_q_paths_1(){
  return q_paths_1;
}

static Q_path *
qop_get_q_paths_2(){
  return q_paths_2;
}

static Q_path *
qop_get_q_paths_3(){
  return q_paths_3;
}



/* Map netdir to the QDP shift dir */

/* Assumes the ordering
   XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN,
   X3UP, Y3UP, Z3UP, T3UP, T3DOWN, Z3DOWN, Y3DOWN, X3DOWN
*/


static QDP_Shift 
fnshift(int netdir){
  if      (netdir <  4) return shiftdirs[netdir];
  else if (netdir <  8) return shiftdirs[7-netdir];
  else if (netdir < 12) return shiftdirs[netdir-4];
  else                  return shiftdirs[19-netdir];
}

static QDP_ShiftDir 
fndir(int netdir){
  if      (netdir <  4) return QDP_forward;
  else if (netdir <  8) return QDP_backward;
  else if (netdir < 12) return QDP_forward;
  else                  return QDP_backward;
}

/* special case to transport a "connection" by one link, does both parities */
static void 
link_transport_connection_qdp( QDP_ColorMatrix *dest, QDP_ColorMatrix *src,
			       QDP_ColorMatrix *gf[4], QDP_ColorMatrix *work,
                               QDP_ColorMatrix *st[8], int dir ){
  if( GOES_FORWARDS(dir) ) {
    QDP_M_eq_M(work, src, QDP_all);
    QDP_M_eq_sM(st[dir], work, QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(dest, gf[dir], st[dir], QDP_all);
    QDP_discard_M(st[dir]);
  }
  else { /* GOES_BACKWARDS(dir) */
    QDP_M_eq_Ma_times_M(work, gf[OPP_DIR(dir)], src, QDP_all);
    QDP_M_eq_sM(st[dir], work, QDP_neighbor[OPP_DIR(dir)], 
		QDP_backward,QDP_all);
    QDP_M_eq_M(dest, st[dir], QDP_all);
    QDP_discard_M(st[dir]);
  }
} /* link_transport_connection_qdp */

static int 
find_backwards_gather( Q_path *path ){
  int disp[4], i;
  /* compute total displacement of path */
  for(i=XUP;i<=TUP;i++)disp[i]=0;
  for( i=0; i<path->length; i++){
    if( GOES_FORWARDS(path->dir[i]) )
      disp[        path->dir[i]  ]++;
    else
      disp[OPP_DIR(path->dir[i]) ]--;
  }
  
  // There must be an elegant way??
  if( disp[XUP]==+1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(XDOWN);
  if( disp[XUP]==-1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(XUP);
  if( disp[XUP]== 0 && disp[YUP]==+1 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(YDOWN);
  if( disp[XUP]== 0 && disp[YUP]==-1 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(YUP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+1 && disp[TUP]== 0 )
    return(ZDOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-1 && disp[TUP]== 0 )
    return(ZUP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+1 )
    return(TDOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-1 )
    return(TUP);
  
  if( disp[XUP]==+3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(X3DOWN);
  if( disp[XUP]==-3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(X3UP);
  if( disp[XUP]== 0 && disp[YUP]==+3 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(Y3DOWN);
  if( disp[XUP]== 0 && disp[YUP]==-3 && disp[ZUP]== 0 && disp[TUP]== 0 )
    return(Y3UP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+3 && disp[TUP]== 0 )
    return(Z3DOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-3 && disp[TUP]== 0 )
    return(Z3UP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+3 )
    return(T3DOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-3 )
    return(T3UP);
  QOP_printf0("OOOPS: NODIR\n"); exit(0);
  return( NODIR );
} //find_backwards_gather


// Make a new path table.  Sorted principally by total displacement of path.
// Below that, sort by direction of first link
// Below that, sort by direction of second link - note special case of one link paths
static int 
sort_quark_paths_hisq( Q_path *src_table, Q_path *dest_table, 
		       int npaths, int num_back_dirs )
{
    int netdir,dir0,dir1,dir1tmp,thislength,num_new,i,j;

    num_new=0; // number of paths in sorted table
    for( i=0; i<num_back_dirs; i++ ){ // loop over net_back_dirs
        netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
	for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
	  for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	    if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	    for( j=0; j<npaths; j++ ){ // pick out paths with right net displacement
	    thislength = src_table[j].length;
	        if( find_backwards_gather( &(src_table[j]) ) == netdir && 
			src_table[j].dir[0]==dir0 &&
			src_table[j].dir[1]==dir1tmp ){
		    dest_table[num_new] = src_table[j];
		    num_new++;
	        }
	    } // loop over paths
	  } //dir1
	} //dir0
    }
    if( num_new!=npaths){ QOP_printf0("OOPS: path table error, num_new=%d, npaths=%d\n",num_new,npaths); exit(0); }
    return 0;
} //sort_quark_paths_hisq


/**********************************************************************/
/*        HISQ force functions                                        */
/**********************************************************************/
// organised as smearing level 0, smearing level i, 
// plus any renormalisation
// and a wrapper
// originally in generic_ks/fermion_force_hisq_multi.c

/* Smearing level 0 Forward Declaration */

static void 
QOPPC(hisq_force_multi_smearing0_fnmat)(QOP_info_t *info,  
					REAL *residues,
					QDP_ColorVector *x[], 
					int nterms, 
					QDP_ColorMatrix *force_accum[4],
					QDP_ColorMatrix *force_accum_naik[4]);


/* Smearing level i Forward Declaration */
static void 
QOPPC(hisq_force_multi_smearing_fnmat)(QOP_info_t *info, 
				       QDP_ColorMatrix * gf[4],
				       REAL *residues,
				       QDP_ColorVector *x[], 
				       int nterms, 
				       QDP_ColorMatrix *force_accum[4],
				       QDP_ColorMatrix *force_accum_old[4],
				       QDP_ColorMatrix 
				        *force_accum_naik_old[4],
				       int internal_num_q_paths,
				       Q_path *internal_q_paths_sorted,
				       int *internal_netbackdir_table
				       );

/* HISQ Force: main wrapper*/

void 
QOPPC(hisq_force_multi_wrapper_fnmat)(QOP_info_t *info,  
				      QOPPC(FermionLinksHisq) *flh,
				      QOP_Force *Force, 
				      QOP_hisq_coeffs_t *hisq_coeff,
				      REAL *residues,
				      QOP_ColorVector *in_pt[], 
				      int *n_orders_naik)
  
{

  int i, ipath, dir;
  REAL coeff_mult;

  double *eps_naik = hisq_coeff->eps_naik;
  int n_naiks = hisq_coeff->n_naiks;
  QOP_hisq_unitarize_method_t umethod = hisq_coeff->umethod;

  // Quark paths sorted by net displacement and last directions
  static Q_path *q_paths_sorted_1 = NULL;
  static Q_path *q_paths_sorted_2 = NULL;
  static Q_path *q_paths_sorted_3 = NULL;

  static int *netbackdir_table_1 = NULL;
  static int *netbackdir_table_2 = NULL;
  static int *netbackdir_table_3 = NULL;

  static int first_force = 1;

  if(first_force == 1) 
    QOP_make_paths_and_dirs_hisq(hisq_coeff, umethod);

  int num_q_paths_1 = qop_get_num_q_paths_1();
  int num_q_paths_2 = qop_get_num_q_paths_2();
  int num_q_paths_3 = qop_get_num_q_paths_3();

  Q_path *q_paths_1 = qop_get_q_paths_1();
  Q_path *q_paths_2 = qop_get_q_paths_2();
  Q_path *q_paths_3 = qop_get_q_paths_3();

  Q_path *q_paths_sorted_current = NULL;
  int *netbackdir_table_current = NULL;

  int inaik;
  int n_naik_shift;

  QDP_ColorMatrix * force[4] =  {Force->force[0], Force->force[1], 
				 Force->force[2], Force->force[3]};

  int num_q_paths_current,n_orders_naik_current;//==nterms


  QDP_ColorMatrix *force_accum_0[4];
  QDP_ColorMatrix *force_accum_0_naik[4];
  QDP_ColorMatrix *force_accum_1[4];
  QDP_ColorMatrix *force_accum_1u[4];
  QDP_ColorMatrix *force_accum_2[4];
  QDP_ColorMatrix *force_final[4];


  QDP_ColorMatrix *Ugf[4], *Vgf[4], *Wgf[4];

  int nterms = 0, n_order_naik_total;

  for(inaik = 0; inaik < n_naiks; inaik++)
    nterms += n_orders_naik[inaik];
  n_order_naik_total = nterms;

  for(i=0;i<4;i++) {
    Ugf[i] = flh->U_links[i];
    Vgf[i] = flh->V_links[i];
    Wgf[i] = flh->W_unitlinks[i];
  }



  QDP_ColorVector * x[(const int) nterms] ;
  for(i=0; i<nterms; i++) x[i] = in_pt[i] -> cv;


  QDP_ColorMatrix *tmat;
  QDP_ColorMatrix *mat_tmp0;

  REAL treal;


  if( first_force==1 ){
    if( q_paths_sorted_1==NULL ) 
      q_paths_sorted_1 = (Q_path *)malloc( num_q_paths_1*sizeof(Q_path) );
    if(netbackdir_table_1==NULL ) 
      netbackdir_table_1 = (int *)malloc( num_q_paths_1*sizeof(int) );
    if( q_paths_sorted_2==NULL ) 
      q_paths_sorted_2 = (Q_path *)malloc( num_q_paths_2*sizeof(Q_path) );
    if(netbackdir_table_2==NULL ) 
      netbackdir_table_2 = (int *)malloc( num_q_paths_2*sizeof(int) );
    if( q_paths_sorted_3==NULL ) 
      q_paths_sorted_3 = (Q_path *)malloc( num_q_paths_3*sizeof(Q_path) );
    if(netbackdir_table_3==NULL ) 
      netbackdir_table_3 = (int *)malloc( num_q_paths_3*sizeof(int) );
    else{QOP_printf0("WARNING: remaking sorted path tables\n"); exit(0); }
    // make sorted tables
    sort_quark_paths_hisq( q_paths_1, q_paths_sorted_1, num_q_paths_1, 8 );

    for( ipath=0; ipath<num_q_paths_1; ipath++ )
      netbackdir_table_1[ipath] = 
	find_backwards_gather( &(q_paths_sorted_1[ipath]) );

    sort_quark_paths_hisq( q_paths_2, q_paths_sorted_2, num_q_paths_2, 16 );

    for( ipath=0; ipath<num_q_paths_2; ipath++ )
      netbackdir_table_2[ipath] = 
	find_backwards_gather( &(q_paths_sorted_2[ipath]) );

    sort_quark_paths_hisq( q_paths_3, q_paths_sorted_3, num_q_paths_3, 16 );

    for( ipath=0; ipath<num_q_paths_3; ipath++ )
      netbackdir_table_3[ipath] = 
	find_backwards_gather( &(q_paths_sorted_3[ipath]) );

    first_force=0;
  }

  tmat = QDP_create_M();
  mat_tmp0 = QDP_create_M();

  for(i=XUP;i<=TUP;i++){
     force_accum_0[i] = QDP_create_M();
     force_accum_0_naik[i] = QDP_create_M();
     force_accum_1[i] = QDP_create_M();
     force_accum_1u[i] = QDP_create_M();
     force_accum_2[i] = QDP_create_M();
     force_final[i] = QDP_create_M();
  }


  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum_2[dir], QDP_all);


  // loop on different naik masses
  n_naik_shift = 0;


  for( inaik=0; inaik<n_naiks; inaik++ ) {

    // smearing level 0
    if( 0==inaik ) {
      n_orders_naik_current = n_order_naik_total;
    }
    else {
      n_orders_naik_current = n_orders_naik[inaik];
    }
    

    QOPPC(hisq_force_multi_smearing0_fnmat)(info,residues+n_naik_shift, 
					    x+n_naik_shift, n_orders_naik_current,
					    force_accum_0, force_accum_0_naik);
 
    
    // smearing level 2
    if( 0==inaik ) {
      q_paths_sorted_current = q_paths_sorted_2;
      num_q_paths_current = num_q_paths_2;
      netbackdir_table_current = netbackdir_table_2;
    }
    else {
      q_paths_sorted_current = q_paths_sorted_3;
      num_q_paths_current = num_q_paths_3;
      netbackdir_table_current = netbackdir_table_3;
    }
    

    QOPPC(hisq_force_multi_smearing_fnmat)( info,Wgf,residues+n_naik_shift, 
					    x+n_naik_shift, 
					    n_orders_naik_current, 
					    force_accum_1, 
					    force_accum_0, force_accum_0_naik, 
					    num_q_paths_current, 
					    q_paths_sorted_current, 
					    netbackdir_table_current );
    

    if( 0==inaik ) {
      coeff_mult = 1.0;
    }
    else {
      coeff_mult = eps_naik[inaik];
    }
    
    
    for(dir=XUP;dir<=TUP;dir++) {
      QDP_M_peq_r_times_M(force_accum_2[dir],&coeff_mult,
			  force_accum_1[dir],QDP_all);
    }
    n_naik_shift += n_orders_naik[inaik];


  }

 

  if ( umethod==QOP_UNITARIZE_NONE ){

    // smearing level 1

    QOPPC(hisq_force_multi_smearing_fnmat)( info,Ugf,residues, 
					    x, 
					    nterms, force_accum_1, 
					    force_accum_2, NULL, 
					    num_q_paths_1, 
					    q_paths_sorted_1, 
					    netbackdir_table_1 );
    
  }
  else if ( umethod==QOP_UNITARIZE_RATIONAL ){

    
    // reunitarization
    QOPPC(hisq_force_multi_reunit)(info,Vgf,force_accum_1u,
					 force_accum_2);
    
    // smearing level 1
    QOPPC(hisq_force_multi_smearing_fnmat)( info,Ugf,residues, 
					    x, 
					    nterms, force_accum_1, 
					    force_accum_1u, NULL, 
					    num_q_paths_1, 
					    q_paths_sorted_1, 
					    netbackdir_table_1 );    

  }
  else
    {
      QOP_printf0("Unknown or unsupported unitarization method\n");
      exit(1);
      
    }


  // contraction with the link in question should be done here,
  // after contributions from all levels of smearing are taken into account

  for(dir=XUP;dir<=TUP;dir++){

    QDP_M_eq_M_times_M(force_final[dir],Ugf[dir],force_accum_1[dir],QDP_all);

  }



  // take into account even/odd parity (it is NOT done in "smearing" routine)
  //eps multiplication done outside QOP 

  for(dir=XUP;dir<=TUP;dir++){
    QDP_M_eq_M(tmat,force_final[dir],QDP_all);

    treal = 2.0;
    QDP_M_eq_r_times_M(force_final[dir],&treal,tmat,QDP_even);

    treal = -2.0;
    QDP_M_eq_r_times_M(force_final[dir],&treal,tmat,QDP_odd);

  }


  // Put antihermitian traceless part into momentum 
  // add force to momentum

  for(dir=XUP; dir<=TUP; dir++){

    QDP_M_eq_antiherm_M(mat_tmp0, force_final[dir], QDP_all);
    QDP_M_peq_M(force[dir], mat_tmp0, QDP_all);

  }



  for(i=XUP;i<=TUP;i++){
     QDP_destroy_M( force_accum_0[i] );
     QDP_destroy_M( force_accum_0_naik[i] );
     QDP_destroy_M( force_accum_1[i] );
     QDP_destroy_M( force_accum_1u[i] );
     QDP_destroy_M( force_accum_2[i] );
     QDP_destroy_M( force_final[i] );
  }

     QDP_destroy_M( tmat );
     QDP_destroy_M( mat_tmp0 );

} //hisq_force_multi_wrapper_fnmat


// like link_transport, except doesn't multiply by link matrices.  
// use this, for example,
// when storing the intermediate HISQ force (a connection) at the lattice site
// associated with a link
static void 
link_gather_connection_qdp( QDP_ColorMatrix *dest, 
			    QDP_ColorMatrix *src,
			    QDP_ColorMatrix *work,
			    int dir ){



  if (dir >= 8) //3 link shift needed
    {
      dir=dir-8;

      //do initial 2 shifts
      if( GOES_FORWARDS(dir) ) {
	
	QDP_M_eq_sM(dest, src, QDP_neighbor[dir], QDP_forward, QDP_all);
	QDP_M_eq_sM(work, dest, QDP_neighbor[dir], QDP_forward, QDP_all);

      }
      else { /* GOES_BACKWARDS(dir) */
	
	QDP_M_eq_sM(dest, src, QDP_neighbor[OPP_DIR(dir)], 
		    QDP_backward, QDP_all);
	QDP_M_eq_sM(work, dest, QDP_neighbor[OPP_DIR(dir)], 
		    QDP_backward, QDP_all);

      }
    }
  else{ //only 1 link shift needed

    QDP_M_eq_M(work, src,  QDP_all);

  }

 

  //do final shift
  if( GOES_FORWARDS(dir) ) {

    QDP_M_eq_sM(dest, work, QDP_neighbor[dir], QDP_forward, QDP_all);

  }
  else { /* GOES_BACKWARDS(dir) */

    QDP_M_eq_sM(dest, work, QDP_neighbor[OPP_DIR(dir)], QDP_backward, QDP_all);

  }


} /* link_gather_connection_qdp */

/* Smearing level 0 */
static void 
QOPPC(hisq_force_multi_smearing0_fnmat)(QOP_info_t *info,  
					REAL *residues,
					QDP_ColorVector *x[], 
					int nterms, 
					QDP_ColorMatrix *force_accum[4],
					QDP_ColorMatrix *force_accum_naik[4])
{


  int term;
  int i,k;
  int dir;
  REAL coeff;


  QDP_ColorMatrix *tmat;
  QDP_ColorMatrix *oprod_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *mat_tmp0;
  QDP_ColorVector *vec_tmp[2];


if( nterms==0 )return;

mat_tmp0   = QDP_create_M();
tmat       = QDP_create_M();
vec_tmp[0] = QDP_create_V();
vec_tmp[1] = QDP_create_V();

for(i=0;i<=MAX_PATH_LENGTH;i++){
  oprod_along_path[i] = QDP_create_M();
}


  // clear force accumulators
  
  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum[dir], QDP_all);


      for(dir=XUP;dir<=TUP;dir++){ //AB loop on directions, path table is not needed

	k=0; // which vec_tmp we are using (0 or 1)
   QDP_V_eq_sV(vec_tmp[k], x[0], 
     fnshift(OPP_DIR(dir)), fndir(OPP_DIR(dir)), QDP_all);
   QDP_M_eq_zero(oprod_along_path[0], QDP_all);


    for(term=0;term<nterms;term++){
      if(term<nterms-1)
	QDP_V_eq_sV(vec_tmp[1-k], x[term+1], 
		    fnshift(OPP_DIR(dir)), fndir(OPP_DIR(dir)), QDP_all);

      QDP_M_eq_V_times_Va(tmat, x[term], vec_tmp[k], QDP_all);
      QDP_discard_V(vec_tmp[k]);
      QDP_M_peq_r_times_M(oprod_along_path[0], &residues[term], tmat, 
			  QDP_all);
      
      k=1-k; // swap 0 and 1

      
    } // end loop over terms in rational function expansion 
    

 
    link_gather_connection_qdp(oprod_along_path[1], oprod_along_path[0], tmat,
			       dir );


    coeff = 1.;


    QDP_M_peq_r_times_M(force_accum[dir],&coeff,oprod_along_path[1],QDP_all);


  } // end of loop on directions //


  // *** Naik part *** /
  
  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum_naik[dir], QDP_all);


  for(dir=XUP;dir<=TUP;dir++){ //AB loop on directions, path table is not needed



    k=0; // which vec_tmp we are using (0 or 1)

    QDP_V_eq_sV(vec_tmp[k], x[0], fnshift(OPP_3_DIR( DIR3(dir) )), 
    	fndir(OPP_3_DIR( DIR3(dir) )), QDP_all);

 
    QDP_M_eq_zero(oprod_along_path[0], QDP_all);

    for(term=0;term<nterms;term++){
      if(term<nterms-1)

      QDP_V_eq_sV(vec_tmp[1-k], x[term+1], fnshift(OPP_3_DIR( DIR3(dir) )), 
		  fndir(OPP_3_DIR( DIR3(dir) )), QDP_all);

      QDP_M_eq_V_times_Va(tmat, x[term], vec_tmp[k], QDP_all);
      QDP_discard_V(vec_tmp[k]);
      QDP_M_peq_r_times_M(oprod_along_path[0], &residues[term], tmat, QDP_all);

      k=1-k; // swap 0 and 1
    } // end loop over terms in rational function expansion 


 

    link_gather_connection_qdp(oprod_along_path[1], oprod_along_path[0], tmat, 
			       DIR3(dir) );

    coeff = 1; // fermion_eps is outside this routine in "wrapper" routine


    QDP_M_peq_r_times_M(force_accum_naik[dir],&coeff,
			oprod_along_path[1],QDP_all);


  } // end of loop on directions 


QDP_destroy_V( vec_tmp[0] );
QDP_destroy_V( vec_tmp[1] );
QDP_destroy_M( mat_tmp0 );
QDP_destroy_M( tmat );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
    QDP_destroy_M( oprod_along_path[i] );
  }
  return;
} //hisq_force_multi_smearing0_fnmat

/* Smearing level i*/
static void 
QOPPC(hisq_force_multi_smearing_fnmat)(QOP_info_t *info, 
				       QDP_ColorMatrix * gf[4],
				       REAL *residues,
				       QDP_ColorVector *x[], 
					int nterms, 
				       QDP_ColorMatrix *force_accum[4],
				       QDP_ColorMatrix *force_accum_old[4],
				       QDP_ColorMatrix 
				        *force_accum_naik_old[4],
				       int internal_num_q_paths,
				       Q_path *internal_q_paths_sorted,
				       int *internal_netbackdir_table
					)

{

  int i,j,k,lastdir=-99,ipath,ilink;
  int length,dir,odir;
  REAL coeff;


  QDP_ColorMatrix *tmat;
  QDP_ColorMatrix *oprod_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *mats_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *mat_tmp0,*mat_tmp1, *stmp[8];;
  QDP_ColorVector *vec_tmp[2];

  int netbackdir, last_netbackdir;

// table of net path displacements (backwards from usual convention)

  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path


  /* Allocate fields */
 for(i=0;i<=MAX_PATH_LENGTH;i++){
   oprod_along_path[i] = QDP_create_M();
 }
 for(i=1;i<=MAX_PATH_LENGTH;i++){ 
    // 0 element is never used (it's unit matrix)
   mats_along_path[i] = QDP_create_M();
 }



 mat_tmp0   = QDP_create_M();
 mat_tmp1   = QDP_create_M();
 for(i=0; i<8; i++) stmp[i] = QDP_create_M();
 tmat       = QDP_create_M();
 vec_tmp[0] = QDP_create_V();
 vec_tmp[1] = QDP_create_V();
 


  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum[dir], QDP_all);

  // loop over paths, and loop over links in path 
  last_netbackdir = NODIR;
  last_path = NULL;



  for( ipath=0; ipath<internal_num_q_paths; ipath++ ){

   this_path = &(internal_q_paths_sorted[ipath]); 

   if(this_path->forwback== -1)continue;	// skip backwards dslash 


    length = this_path->length;
    netbackdir = internal_netbackdir_table[ipath];

    // move f(i-1) force from current site in positive direction,
    //  this corresponds to outer product |X><Y| calculated at the endpoint of the path 


    if( netbackdir<8) { // Not a Naik path
      link_gather_connection_qdp(oprod_along_path[0] , 
				 force_accum_old[OPP_DIR(netbackdir)],
				 tmat, 
				 netbackdir );
    }
    else { // Naik path
      if( NULL==force_accum_naik_old ) {
        QOP_printf0( "hisq_force_multi_smearing_fnmat:  mismatch:\n" );
        QOP_printf0( "force_accum_naik_old is NULL, but path table contains Naik paths(!)\n" );
        exit(0);
      }
      // CONVERSION FROM 3-LINK DIRECTION TO 1-LINK DIRECTION

      link_gather_connection_qdp(oprod_along_path[0] , 
      			 force_accum_naik_old[OPP_DIR(netbackdir-8)],
				 tmat,
      			 netbackdir );

    }


    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;


    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection_qdp( oprod_along_path[length-ilink], 
				     oprod_along_path[length-ilink-1], gf,
				     mat_tmp0, stmp, this_path->dir[ilink]  );
    }



    // maintain an array of transports "to this point" along the path.
    //	Don't recompute beginning parts of path if same as last path 
    ilink=0; // first link where new transport is needed
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){


	    QDP_M_eq_sM(tmat, gf[dir], QDP_neighbor[dir],
			QDP_backward, QDP_all);
	    QDP_M_eq_Ma(mats_along_path[1], tmat, QDP_all);
	    QDP_discard_M(tmat);

          }
          else{

	    QDP_M_eq_M(mats_along_path[1], gf[OPP_DIR(dir)], QDP_all);

          }


      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);

	link_transport_connection_qdp( mats_along_path[ilink+1], 
				       mats_along_path[ilink], gf,
				       mat_tmp0, stmp, dir );

      }
    } // end loop over links



    // A path has (length+1) points, counting the ends.  At first
    //	 point, no "down" direction links have their momenta "at this
    //	 point". At last, no "up" ... 
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;



      coeff = this_path->coeff;

      if( (ilink%2)==1 )coeff = -coeff;



      // add in contribution to the force 
      if( ilink<length && GOES_FORWARDS(dir) ){


	link_gather_connection_qdp(mat_tmp1, 
		       oprod_along_path[length-ilink-1], tmat, dir );

        if(ilink==0) 
	  {
	    QDP_M_eq_M(mat_tmp0,mat_tmp1,QDP_all);
	  }
        else
	  {

	    QDP_M_eq_M_times_Ma(mat_tmp0, mats_along_path[ilink], 
				mat_tmp1, QDP_all);
	    QDP_M_eq_Ma(mat_tmp1,mat_tmp0,QDP_all);


	  }

	QDP_M_peq_r_times_M(force_accum[dir],&coeff,mat_tmp1,QDP_all);
	

      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        if( ilink==1 ){
	  QDP_M_eq_M(mat_tmp0,oprod_along_path[length-ilink],QDP_all);
	  QDP_M_eq_Ma(mat_tmp1,mat_tmp0,QDP_all);
	}
        else{

	  link_gather_connection_qdp(mat_tmp1, mats_along_path[ilink-1], 
				     tmat, odir );


	  QDP_M_eq_M_times_Ma(mat_tmp0, oprod_along_path[length-ilink], 
			      mat_tmp1, QDP_all);
	  QDP_M_eq_Ma(mat_tmp1, mat_tmp0, QDP_all);


        }


	QDP_M_peq_r_times_M(force_accum[odir],&coeff,mat_tmp1,QDP_all);

      }
      lastdir = dir;
    } // end loop over links in path //

  } // end loop over paths //



  QDP_destroy_V( vec_tmp[0] );
  QDP_destroy_V( vec_tmp[1] );
  QDP_destroy_M( mat_tmp0 );
  QDP_destroy_M( mat_tmp1 );
  QDP_destroy_M( tmat );
  for(i=0; i<8; i++) QDP_destroy_M(stmp[i]);
  for(i=0;i<=MAX_PATH_LENGTH;i++){
    QDP_destroy_M( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
    QDP_destroy_M( mats_along_path[i] );
  }

  return;
}//hisq_force_multi_smearing_fnmat

