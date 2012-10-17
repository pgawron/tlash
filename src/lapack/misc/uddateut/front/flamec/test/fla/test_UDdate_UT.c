/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2011, The University of Texas

   libflame is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   libflame is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with libflame; if you did not receive a copy, see
   http://www.gnu.org/licenses/.

   For more information, please contact us at flame@cs.utexas.edu or
   send mail to:

   Field G. Van Zee and/or
   Robert A. van de Geijn
   The University of Texas at Austin
   Department of Computer Sciences
   1 University Station C0500
   Austin TX 78712
*/

#include "FLAME.h"

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

void time_UDdate_UT(
                 int variant, int type, int n_repeats, int mB, int mC, int mD, int n,
                 FLA_Obj B, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj R, FLA_Obj E,
                 double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_input,
    mB_input, mC_input, mD_input,
    mB, mC, mD, n,
    p_first, p_last, p_inc,
    p,
    b_alg,
    variant,
    n_repeats,
    i,
    n_variants = 1;
  
  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    B, C, D, T, R, E;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter n (-1 means bind to problem size): ", '%' );
  scanf( "%d", &n_input );
  fprintf( stdout, "%c %d\n", '%', n_input );

  fprintf( stdout, "%c enter mB mC mD (-1 means bind to problem size): ", '%' );
  scanf( "%d %d %d", &mB_input, &mC_input, &mD_input );
  fprintf( stdout, "%c %d %d %d\n", '%', mB_input, mC_input, mD_input );


  fprintf( stdout, "\nclear all;\n\n" );



  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    mB = mB_input;
    mC = mC_input;
    mD = mD_input;
    n  = n_input;

    if( mB < 0 ) mB = p / abs(mB_input);
    if( mC < 0 ) mC = p / abs(mC_input);
    if( mD < 0 ) mD = p / abs(mD_input);
    if( n  < 0 ) n  = p / abs(n_input);

    for ( variant = 0; variant < n_variants; variant++ ){
      
      FLA_Obj_create( datatype, mB, n, 0, 0, &B );
      FLA_Obj_create( datatype, mC, n, 0, 0, &C );
      FLA_Obj_create( datatype, mD, n, 0, 0, &D );
      FLA_Obj_create( datatype, b_alg, n, 0, 0, &T );
      FLA_Obj_create( datatype, n,  n, 0, 0, &R );
      FLA_Obj_create( datatype, n,  n, 0, 0, &E );

      FLA_Random_matrix( B );
      FLA_Random_matrix( C );
      FLA_Random_matrix( D );

      FLA_Set( FLA_ZERO, R );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D, FLA_ONE, R );
      FLA_Chol( FLA_UPPER_TRIANGULAR, R );

      FLA_Set( FLA_ZERO, E );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, E );
      FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C, FLA_ONE, E );
      FLA_Chol( FLA_UPPER_TRIANGULAR, E );

      fprintf( stdout, "data_uddate_ut( %d, 1:5 ) = [ %d  ", i, p );
      fflush( stdout );

      time_UDdate_UT( variant, FLA_ALG_FRONT, n_repeats, mB, mC, mD, n,
                      B, C, D, T, R, E, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &B );
      FLA_Obj_free( &C );
      FLA_Obj_free( &D );
      FLA_Obj_free( &T );
      FLA_Obj_free( &R );
      FLA_Obj_free( &E );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_variants; i++ ) {
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 2 ), '%c:%c' ); \n",
            colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_qr_ut( :,1 ), data_qr_ut( :, 4 ), '%c-.%c' ); \n",
            colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_variants; i++ )
    fprintf( stdout, "'ref\\_qr\\_ut', 'fla\\_qr\\_ut', ... \n" );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME UDdate_UT front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc qr_ut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

