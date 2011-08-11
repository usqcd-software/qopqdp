/**********************************************************************/
/*   Originally MILC dirs.h                                           */
/**********************************************************************/

#define MAX_PATH_LENGTH 16
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15

#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir) ((dir)+8)
#define FORALL3UPDIR(dir) for(dir=X3UP; dir<=T3UP; dir++)


#define NODIR -1  /* not a direction */

#define OPP_DIR(dir)	(7-(dir))	/* Opposite direction */
#define NDIRS 8				/* number of directions */


static int net_back_dirs[16] = { XDOWN, YDOWN, ZDOWN, TDOWN, 
				 XUP, YUP, ZUP, TUP, 
				 X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, 
				 X3UP, Y3UP, Z3UP, T3UP };

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)



/**********************************************************************/
/*   Originally MILC asqtad_action.h                                  */
/**********************************************************************/

    /* The fat link action with seven link paths designed to zero
	couplings at momentum pi in any direction. The term introduced
	by Lepage to cancel the additional O(a^2) errors introduced
	by the fattening is added.  The Naik term corrects the dispersion
	relation	 */
    /* Specify paths in orientation in which they appear in the
	forward part of the x component of dslash().  Rotations and
	reflections will be automatically included. Be careful
	about signs of coefficients.  See long comment at bottom
	of quark_stuff.c. */


//AG undefine this just now to match MILC code
//#define TADPOLE_IMPROVE /* use tadpole improvement in quark action */
//static REAL u0 = 0.890;


#define ASQ_OPTIMIZED_FATTENING
#define ASQ_OPTIMIZED_FORCE

/* We don't include the explicit path coefficients, since they are
   being passed into QOP by the caller */


/**********************************************************************/
/*   Originally MILC hisq_su3_action.h                                */
/**********************************************************************/


//#define HISQ_NAIK_ADJUSTABLE // allow for adjustable epsilon in Naik term
#define HISQ_NAIK_2ND_ORDER (-0.675)

#define MAX_NUM 688  // should be obsolete, for now max of MAX_NUM_[12]
#define QUARK_ACTION_DESCRIPTION "\"HISQ action version 1\""

#define MAX_LENGTH 7	// Maximum length of path in any path table
#define MAX_BASIC_PATHS 6  // Max. no. of basic paths in any path table

// Smearing for first level
// This is Fat7
#define NUM_BASIC_PATHS_1 4
#define MAX_NUM_1 632
//#define ASQ_OPTIMIZED_FATTENING_1
//#define ASQ_OPTIMIZED_FORCE_1
#define QUARK_ACTION_DESCRIPTION_1 "\"Fat 7 (level 1)\""
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    };
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1,3,5,7};
/* We don't include the explicit path coefficients, since they are
   being passed into QOP by the caller */

// Unitarization algorithm
// Choices are:
//   UNITARIZE_NONE

// This option is set with QOP_hisq_coeffs_t
//#define UNITARIZE_NONE 0
//#define UNITARIZE_RATIONAL 1

// This option is set with QOP_hisq_links_set_opts.
//#define UNITARIZATION_METHOD UNITARIZE_NONE
//#define UNITARIZATION_METHOD UNITARIZE_RATIONAL

// Smearing for second level
#define NUM_BASIC_PATHS_2 6
#define MAX_NUM_2 688
//#define ASQ_OPTIMIZED_FATTENING_2
//#define ASQ_OPTIMIZED_FORCE_2
#define QUARK_ACTION_DESCRIPTION_2 "\"Fat7 + 2xLepage\""
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },      /* 5-link compensation    */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1,3,3,5,7,5};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
/* We don't include the explicit path coefficients, since they are
   being passed into QOP by the caller */
#define INDEX_ONELINK 0
#define INDEX_NAIK 1


// Smearing for second level -- only difference in 1-link and Naik
#define NUM_BASIC_PATHS_3 2
#define MAX_NUM_3 16
#define QUARK_ACTION_DESCRIPTION_3 "\"1-link + Naik\""
    static int path_ind_3[NUM_BASIC_PATHS_3][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    };
    static int path_length_in_3[NUM_BASIC_PATHS_3] = {1,3};
    static int quark_action_npaths_3 = NUM_BASIC_PATHS_3 ;

/* We don't include the explicit path coefficients, since they are
   being passed into QOP by the caller */


