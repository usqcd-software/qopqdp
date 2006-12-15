/**************** fermion_force_fn_qdp.c ******************************/
/* MIMD version 7 */
/* General force for any FN-type action.
 * Optimized to transport only one set of SU(3) matrices.
 * D.T. 12/05 Version  3. created for improved fermion RHMC.
 * C.D. 10/06 Converted to QDP (original code in fermion_force_fn_multi.c)
 *      Added multi vector version.	
 */
#include <qop_internal.h>
#include <string.h>
#include <asqtad_action.h>

/** #define DEBUG_FNMAT **/

/********************************************************************/
/* Construct path table from the action definition above
   Originally quark_stuff.c                                         */
/********************************************************************/
						
typedef struct {
  int dir[MAX_PATH_LENGTH];	/* directions in path */
  int length;		/* length of path */
  REAL coeff;	/* coefficient, including minus sign if backwards */
  REAL forwback;	/* +1 if in forward Dslash, -1 if in backward */
} Q_path;

static QDP_Shift shiftdirs[8];
static QDP_Shift neighbor3[4];

static Q_path q_paths[MAX_NUM];
static int num_q_paths; /* number of paths in dslash */
static int num_basic_paths;	/* number of paths before rotation/reflection */

static int qop_is_path_equal( int *path1, int* path2, int length );
static int qop_add_basic_path( int *vec, int length, REAL coeff );

static void 
QOP_make_paths_and_dirs(QOP_asqtad_coeffs_t *coef) {

  int i,j;
  int disp[4] = {0,0,0,0};

  /* table of directions, 1 for each kind of path */
  /**int path_ind[MAX_BASIC_PATHS][MAX_LENGTH];**/
  /* table of coefficients in action, for each path */
  for(i=0; i<4; i++) {
    disp[i] = 3;
    neighbor3[i] = QDP_create_shift(disp);
    disp[i] = 0;
  }

  for(i=0; i<4; ++i) {
    shiftdirs[i] = QDP_neighbor[i];
    shiftdirs[i+4] = neighbor3[i];
  }
  
#ifdef DEBUG_FNMAT
  QOP_printf0("%s\n",QUARK_ACTION_DESCRIPTION);
#endif
  num_q_paths = 0;
  num_basic_paths = 0;

  if(MAX_LENGTH > MAX_PATH_LENGTH){
    printf("Path length for this action is too long.	Recompile.\n");
    exit(0);
  }
  
  REAL path_coeff[MAX_BASIC_PATHS] = {coef->one_link, coef->naik, 
				      coef->three_staple, coef->five_staple,
				      coef->seven_staple, coef->lepage};  
  /* add rots. and reflects to table, print out the path coefficients */

  for(j=0;j<quark_action_npaths;j++) {
    REAL this_coeff;
    this_coeff = path_coeff[j];
    //QOP_printf0("this_coeff = %e\n", this_coeff);
    
    i = qop_add_basic_path( path_ind[j], path_length_in[j],
			this_coeff );
  }
}

static int 
qop_get_num_q_paths(){
  return num_q_paths;
}

static Q_path *
qop_get_q_paths(){
  return q_paths;
}



/********************************************************************/
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */
/********************************************************************/
static int 
qop_add_basic_path( int *basic_vec, int length, REAL coeff ) {

  int perm[8],pp[8],ir[4];
  int j,path_num;
  int vec[MAX_LENGTH];
  int flag;
  
  path_num = 0;
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
		    for(j=0;j<num_q_paths;j++){
		      flag = qop_is_path_equal( vec, q_paths[j].dir, MAX_LENGTH );
		      if(flag==1)break;
		    }
		    if(flag == 0 ){
		      if(num_q_paths>=MAX_NUM){
			QOP_printf0("OOPS: MAX_NUM too small\n");
			exit(0);
		      }
		      q_paths[num_q_paths].length=length;
		      for(j=0;j<MAX_LENGTH;j++) q_paths[num_q_paths].dir[j]=vec[j];
		      /* remember to copy NODIR's, or comparison will get confused */
		      if(ir[0]==0){
			q_paths[num_q_paths].coeff =  coeff;

			q_paths[num_q_paths].forwback =	 +1;
		      }
		      else{
			q_paths[num_q_paths].coeff = -coeff;
			q_paths[num_q_paths].forwback = -1;
		      }
		      num_q_paths++;
		      path_num++;
		      /**node0_printf("ADD PATH %d:  rx=%d ",num_q_paths-1,ir[0]);
			 printpath( vec, length );**/
		    }
		    
		  } /* end reflection*/
	  } /* end permutation if block */
	} /* end permutation */
  num_basic_paths++;
  return(path_num);
} /* add_basic_path */

static int 
qop_is_path_equal( int *path1, int* path2, int length ){
   register int i;
   for(i=0;i<length;i++)if(path1[i]!=path2[i])return(0);
   return(1);
}


/**********************************************************************/
/*   Utilities                                                        */
/**********************************************************************/

#ifdef DEBUG_FNMAT

static void 
PrintM(QDP_ColorMatrix * field){
  QLA_ColorMatrix *mom0;
  const int x[4] = {0,0,0,0};
  int ind = QDP_index(x);

  mom0 = QDP_expose_M(field);
  printf("Field: %e %e %e\n\n", mom0[ind].e[0][0].real,
	 mom0[ind].e[0][1].real, mom0[ind].e[0][2].real);

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
  if( disp[XUP]==+1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XDOWN);
  if( disp[XUP]==-1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XUP);
  if( disp[XUP]== 0 && disp[YUP]==+1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YDOWN);
  if( disp[XUP]== 0 && disp[YUP]==-1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YUP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+1 && disp[TUP]== 0 )return(ZDOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-1 && disp[TUP]== 0 )return(ZUP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+1 )return(TDOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-1 )return(TUP);
  
  if( disp[XUP]==+3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3DOWN);
  if( disp[XUP]==-3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3UP);
  if( disp[XUP]== 0 && disp[YUP]==+3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3DOWN);
  if( disp[XUP]== 0 && disp[YUP]==-3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3UP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+3 && disp[TUP]== 0 )return(Z3DOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-3 && disp[TUP]== 0 )return(Z3UP);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+3 )return(T3DOWN);
  if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-3 )return(T3UP);
  QOP_printf0("OOOPS: NODIR\n"); exit(0);
  return( NODIR );
} //find_backwards_gather


// Make a new path table.  Sorted principally by total displacement of
// path.  Below that, sort by direction of first link 
// Below that, sort by direction of second link - note special case of
// one link paths
static int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths ){
  int netdir,dir0,dir1,dir1tmp,thislength,num_new,i,j;
  
  num_new=0; // number of paths in sorted table
  for( i=0; i<16; i++ ){ // loop over net_back_dirs
    netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
    for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
      for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	for( j=0; j<npaths; j++ ){
	  // pick out paths with right net displacement
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
  if( num_new!=npaths){ QOP_printf0("OOPS: path table error\n"); exit(0); }
  return 0;
}

/**********************************************************************/
/*   General QOP FN Version for "nterms" sources			*/
/**********************************************************************/

void 
QOPPC(asqtad_force_multi_fnmat)(QOP_info_t *info,  QOP_GaugeField *Gauge,
			QOP_Force *Force, QOP_asqtad_coeffs_t *asq_coeff,
			REAL eps[], QOP_ColorVector *in_pt[], int nterms)
{
#ifdef DEBUG_FNMAT
  printf("coeff: %e %e %e\n %e %e %e\n\n", asq_coeff->one_link, asq_coeff->naik,
	 asq_coeff->three_staple, asq_coeff->five_staple,
	 asq_coeff->seven_staple, asq_coeff->lepage);
#endif
  
  /* note CG_solution and Dslash * solution are combined in "x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x transported from both ends of path. */
  int term;
  int i,j,k,lastdir=-99,ipath,ilink;
  int length,dir,odir;

  QDP_ColorMatrix * force[4] =  {Force->force[0], Force->force[1], 
				 Force->force[2], Force->force[3]};
#ifdef DEBUG_FNMAT
  PrintM(force[0]);
#endif

  QDP_ColorMatrix * gf[4] = {Gauge->links[0], Gauge->links[1], 
			     Gauge->links[2], Gauge->links[3]};
#ifdef DEBUG_FNMAT
  PrintM(gf[0]);
#endif

  QDP_ColorVector * x[(const int) nterms] ;
  for(i=0; i<nterms; i++) x[i] = in_pt[i] -> cv;
#ifdef DEBUG_FNMAT
  PrintV(x[0]);  
#endif
   
  QDP_ColorMatrix *tmat;
  REAL coeff;

  QDP_ColorMatrix *mat_tmp0, *stmp[8];
  QDP_ColorMatrix *oprod_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *mats_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *force_accum[4];  // accumulate force
  QDP_ColorVector *vec_tmp[2];

  int netbackdir, last_netbackdir;	// backwards direction for entire path
  //int tempflops = 0; //TEMP
  // Quark paths sorted by net displacement and last directions
  static Q_path *q_paths_sorted = NULL; 
  // table of net path displacements (backwards from usual convention)
  static int *netbackdir_table = NULL;
  static int first_force = 1;  // 1 if force hasn't been called yet

  if(first_force == 1) QOP_make_paths_and_dirs(asq_coeff);
  int num_q_paths = qop_get_num_q_paths();
  Q_path *q_paths = qop_get_q_paths();
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path


  char myname[] = "QOP_asqtad_force_fnmat";


  double dtime;

  if( nterms==0 )return;

  dtime=-QOP_time();
	
  /* Allocate fields */
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = QDP_create_M();
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ 
    // 0 element is never used (it's unit matrix)
     mats_along_path[i] = QDP_create_M();
  }
  for(i=XUP;i<=TUP;i++){
     force_accum[i] = QDP_create_M();
  }
  mat_tmp0   = QDP_create_M();
  for(i=0; i<8; i++) stmp[i] = QDP_create_M();
  tmat       = QDP_create_M();
  vec_tmp[0] = QDP_create_V();
  vec_tmp[1] = QDP_create_V();
  if( vec_tmp[1] == NULL ){
    QOP_printf0("%s NO ROOM\n",myname); 
    exit(0);
  }
  
  /* Sort the paths */
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) 
      q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    if(netbackdir_table==NULL ) 
      netbackdir_table = (int *)malloc( num_q_paths*sizeof(int) );
    else{ QOP_printf0("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    for( ipath=0; ipath<num_q_paths; ipath++ ){
      netbackdir_table[ipath] = 
	find_backwards_gather( &(q_paths_sorted[ipath]) );
      //QOP_printf0 ("sortedpath coeff = %e\n", q_paths_sorted[ipath].coeff);
    }
    first_force=0;
  }
  
  //  QOP_printf0 ("\n\n sortedpath coeff = %e\n", q_paths_sorted[13].coeff);
  //  QOP_printf0 (" length = %i\n", q_paths_sorted[13].length);
  //  QOP_printf0 (" fwdbck = %e\n", q_paths_sorted[13].forwback);
  //  QOP_printf0 (" dir6 = %i\n", q_paths_sorted[13].dir[6]);
  //  QOP_printf0 (" netbackdir = %i\n\n", netbackdir_table[13]);
  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum[dir], QDP_all);

  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    this_path = &(q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */
    
    length = this_path->length;
    //QOP_printf0("ipath = %i length = %i\n", ipath, length);
    // find gather to bring x[term] from "this site" to end of path
    //netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
    netbackdir = netbackdir_table[ipath];
    // and bring x to end - no gauge transformation !!
    // The resulting outer product matrix has gauge transformation
    // properties of a connection from start to end of path
    if( netbackdir != last_netbackdir){ 
      // don't need to repeat this if same net disp. as last path
      k=0; // which vec_tmp we are using (0 or 1)
      QDP_V_eq_sV(vec_tmp[k], x[0], 
	  fnshift(netbackdir), fndir(netbackdir), QDP_all);
      //PrintV(vec_tmp[k]);
      // actually last site in path
      QDP_M_eq_zero(oprod_along_path[0], QDP_all);
      
      for(term=0;term<nterms;term++){
	if(term<nterms-1){
	  QDP_V_eq_sV(vec_tmp[1-k], x[term+1], 
		      fnshift(netbackdir), fndir(netbackdir), QDP_all);
	}
	QDP_M_eq_V_times_Va(tmat, x[term], vec_tmp[k], QDP_all);
	QDP_discard_V(vec_tmp[k]);
	QDP_M_peq_r_times_M(oprod_along_path[0], &eps[term], tmat, QDP_all);
	//PrintM(oprod_along_path[0]);
	k=1-k; // swap 0 and 1
      } /* end loop over terms in rational function expansion */
    }
    //tempflops+=54*nterms;
    //tempflops+=36*nterms;

    /* path transport the outer product, or projection matrix, of x[term]
       (EVEN sites)  and Dslash*x[term] (ODD sites) from far end.
       
       maintain a matrix of the outer product transported backwards
       along the path to all sites on the path.
       If new "net displacement", need to completely recreate it.
       Otherwise, use as much of the previous path as possible 
       
       Note this array is indexed backwards - the outer product transported
       to site number n along the path is in oprod_along_path[length-n].
       This makes reusing it for later paths easier.
       
       Sometimes we need this at the start point of the path, and sometimes
       one link into the path, so don't always have to do the last link. */
    
    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( netbackdir == last_netbackdir )
      while ( j>0 && 
	      this_path->dir[j] == 
	      last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;
    
    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection_qdp( oprod_along_path[length-ilink], 
				     oprod_along_path[length-ilink-1], gf,
				     mat_tmp0, stmp, this_path->dir[ilink]  );
      //tempflops+=9*22;
    }

    /* maintain an array of transports "to this point" along the path.
       Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    if( last_path != NULL )
      while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
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
	//tempflops+=9*22;
      }
    } // end loop over links

    /* A path has (length+1) points, counting the ends.  At first
       point, no "down" direction links have their momenta "at this
       point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;
      coeff = this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;
      //QOP_printf0("coeff = %e\n", coeff);

      if(ilink==0 && GOES_FORWARDS(dir) ){
	QDP_M_eq_M(mat_tmp0, oprod_along_path[length],QDP_all); 
      }
      else if( ilink>0){
	QDP_M_eq_M_times_Ma(mat_tmp0, oprod_along_path[length-ilink],  
			    mats_along_path[ilink], QDP_all);
      }
      //if(ilink>0)tempflops+=9*22;

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      if( ilink<length && GOES_FORWARDS(dir) ){
	QDP_M_peq_r_times_M(force_accum[dir], &coeff, mat_tmp0, QDP_even);
	QDP_M_meq_r_times_M(force_accum[dir], &coeff, mat_tmp0, QDP_odd);
	//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
	QDP_M_meq_r_times_M(force_accum[odir], &coeff, mat_tmp0, QDP_even);
	QDP_M_peq_r_times_M(force_accum[odir], &coeff, mat_tmp0, QDP_odd);
      }
      //tempflops+=36;
      
      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted[ipath]);
  } /* end loop over paths */

  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++){
    QDP_M_eq_antiherm_M(mat_tmp0, force_accum[dir], QDP_all);
    QDP_M_peq_M(force[dir], mat_tmp0, QDP_all);
  }
  //tempflops+=4*18;
  //tempflops+=4*18;
  
  QDP_destroy_V( vec_tmp[0] );
  QDP_destroy_V( vec_tmp[1] );
  QDP_destroy_M( mat_tmp0 );
  for(i=0; i<8; i++) QDP_destroy_M(stmp[i]);
  QDP_destroy_M( tmat );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     QDP_destroy_M( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     QDP_destroy_M( mats_along_path[i] );
  }
  for(i=XUP;i<=TUP;i++){
     QDP_destroy_M( force_accum[i] );
  }

  dtime += QOP_time();
  int nflop = 966456 + 1440*nterms;
  info->final_sec = dtime;
  info->final_flop = nflop*QDP_sites_on_node;
  info->status = QOP_SUCCESS;
}



