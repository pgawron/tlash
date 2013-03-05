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

#include "plasma.h"
#include "FLAME.h" 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 


#define OUTPUT_PATH "./results"
#define OUTPUT_FILE "chol_plasma"


int main( int argc, char *argv[] ) 
{ 
   int
      i, j,
      size,
      n_threads,
      n_repeats,
      n_trials,
      nb_alg,
      increment,
      begin;

   FLA_Datatype
      datatype = FLA_DOUBLE;
   
   FLA_Obj 
      A, x, b, b_norm;
   
   double 
      length,
      b_norm_value,
      dtime,
      *dtimes,
      *flops;

   char
      output_file_m[100];
   
   FILE
      *fpp;

   fprintf( stdout, "%c Enter number of repeats: ", '%' );
   scanf( "%d", &n_repeats );
   fprintf( stdout, "%c %d\n", '%', n_repeats );

   fprintf( stdout, "%c Enter blocksize: ", '%' );
   scanf( "%d", &nb_alg );
   fprintf( stdout, "%c %d\n", '%', nb_alg );

   fprintf( stdout, "%c Enter problem size parameters: first, inc, num: ", '%' );
   scanf( "%d%d%d", &begin, &increment, &n_trials );
   fprintf( stdout, "%c %d %d %d\n", '%', begin, increment, n_trials );

   fprintf( stdout, "%c Enter number of threads: ", '%' );
   scanf( "%d", &n_threads );
   fprintf( stdout, "%c %d\n\n", '%', n_threads );

   sprintf( output_file_m, "%s/%s_output.m", OUTPUT_PATH, OUTPUT_FILE );
   fpp = fopen( output_file_m, "a" );

   fprintf( fpp, "%%\n" );
   fprintf( fpp, "%% | Matrix Size |    PLASMA   |\n" );
   fprintf( fpp, "%% |    n x n    |    GFlops   |\n" );
   fprintf( fpp, "%% -----------------------------\n" );
   
   FLA_Init();
   PLASMA_Init( n_threads );

   PLASMA_Disable( PLASMA_AUTOTUNING );
   PLASMA_Set( PLASMA_TILE_SIZE, nb_alg );
   PLASMA_Set( PLASMA_INNER_BLOCK_SIZE, nb_alg / 4 );

   dtimes = ( double * ) FLA_malloc( n_repeats * sizeof( double ) );
   flops  = ( double * ) FLA_malloc( n_trials  * sizeof( double ) );

   fprintf( fpp, "%s = [\n", OUTPUT_FILE );

   for ( i = 0; i < n_trials; i++ )
   {
      size = begin + i * increment;
      
      FLA_Obj_create( datatype, size, size, 0, 0, &A ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &x ); 
      FLA_Obj_create( datatype, size, 1,    0, 0, &b ); 
      FLA_Obj_create( datatype, 1,    1,    0, 0, &b_norm ); 
      
      for ( j = 0; j < n_repeats; j++ )
      {
         FLA_Random_matrix( A );
         FLA_Random_matrix( x );
         FLA_Random_matrix( b );

         length = ( double ) FLA_Obj_length( A );
         FLA_Add_to_diag( &length, A );

         FLA_Symv_external( FLA_LOWER_TRIANGULAR, FLA_ONE, A, x, FLA_ZERO, b );

         dtime = FLA_Clock();

         PLASMA_dpotrf( PlasmaLower, size, FLA_Obj_buffer_at_view( A ), size );

         dtime = FLA_Clock() - dtime;
         dtimes[j] = dtime;      
      }
      
      dtime = dtimes[0];
      for ( j = 1; j < n_repeats; j++ )
         dtime = min( dtime, dtimes[j] );
      flops[i] = 1.0 / 3.0 * size * size * size / dtime / 1e9;
      
      fprintf( fpp, "   %d   %6.3f\n", size, flops[i] );
      
      FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, 
                         FLA_NONUNIT_DIAG, A, b );
      FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, 
                         FLA_NONUNIT_DIAG, A, b ); 
      FLA_Axpy_external( FLA_MINUS_ONE, x, b );
      FLA_Nrm2_external( b, b_norm );
      FLA_Obj_extract_real_scalar( b_norm, &b_norm_value );
      
      printf( "Time: %e  |  GFlops: %6.3f\n", 
              dtime, flops[i] );
      printf( "Matrix size: %d x %d  |  nb_alg: %d\n", 
              size, size, nb_alg ); 
      printf( "Norm of difference: %le\n\n", b_norm_value ); 
      
      FLA_Obj_free( &A ); 
      FLA_Obj_free( &x ); 
      FLA_Obj_free( &b ); 
      FLA_Obj_free( &b_norm ); 
   }

   fprintf( fpp, "];\n" );
   
   fflush( fpp );
   fclose( fpp );

   FLA_free( dtimes );
   FLA_free( flops );

   PLASMA_Finalize();
   FLA_Finalize(); 
   
   return 0; 
}
