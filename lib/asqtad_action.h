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
#define MAX_BASIC_PATHS 6
#define MAX_LENGTH 7
#define MAX_NUM 688
#define TADPOLE_IMPROVE /* use tadpole improvement in quark action */
#define ASQ_OPTIMIZED_FATTENING
#define ASQ_OPTIMIZED_FORCE
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation	  */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5,7,5};
    static int quark_action_npaths = MAX_BASIC_PATHS ;

    static char quark_action_description[] =
	"\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\"";

/* We don't include the explicit path coefficients, since they are
   being passed into QOP by the caller */
