#include <test_common.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <qop_internal.h>


// amplitude of the random matrix
#define MAX_AMPLITUDE 10
// number of repetitions
//#define MAX_REPEAT 10000
#define MAX_REPEAT 2



static void print_mat( QLA_ColorMatrix *A ) {
  int i,j;
  for( i=0; i<3; i++) {
    for( j=0; j<3; j++) {
      printf("(%+.4e,%+.4e)", QLA_real(QLA_elem_M(*A,i,j)), QLA_imag(QLA_elem_M(*A,i,j)));
      if( j!=2 ) printf("  ");
    }
    printf("\n");
  }
  printf("\n");
}


int
main(int argc, char *argv[]) {
  int irpt, i, j;
  QLA_RandomState *s;
  QLA_ColorMatrix A;
  QLA_Complex det;

  // initialize QDP -- SO FAR NOT NEEDED?

  // allocate memory for random number generator
  s=(QLA_RandomState*)malloc(sizeof(QLA_RandomState));

  /* initialize random number generator */
  QLA_seed_random( s, 0, 0 );

  // unit matrix
//  unity_su3mat( &Unit );

#if ( QOP_Precision==1 )
  printf( "Testing in SINGLE precision\n" );
#else
  printf( "Testing in DOUBLE precision\n" );
#endif

  /* MAIN REPETITION LOOP */
  for( irpt=0; irpt<MAX_REPEAT; irpt++ ) {
    /* randomize matrix */
    for( i=0; i<3; i++) {
      for( j=0; j<3; j++) {
#if 0
	A.e[i][j].real
            = MAX_AMPLITUDE*( 2 * QLA_random(s) - 1 );
        A.e[i][j].imag 
            = MAX_AMPLITUDE*( 2 * QLA_random(s) - 1 );
#endif
	QLA_c_eq_r_plus_ir(QLA_elem_M(A,i,j),
			   MAX_AMPLITUDE*( 2 * QLA_random(s) - 1 ),
			   MAX_AMPLITUDE*( 2 * QLA_random(s) - 1 ));
      }
    }
    printf( "*** Iteration %d\n", irpt );
    printf( "Random 3x3 complex matrix A:\n" );
    print_mat( &A );
    det=QOPPC(su3_mat_det)( &A );
    printf( "det(A): (%18.10g,%18.10g)\n\n", det.real, det.imag );
  }


  // free random number generator
  free(s);
}
