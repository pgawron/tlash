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

#define N_PARAM_COMBOS    2

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PARAM_COMBOS] = { "l", "u" };

void time_Chol(
                int param_combo, int type, int n_repeats, int n,
                FLA_Obj A, FLA_Obj A_ref,
                double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_threads,
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;
  
  dim_t
    nb_alg;

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, A_ref;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &nb_alg );
  fprintf( stdout, "%c %u\n", '%', nb_alg );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%d", &n_threads );
  fprintf( stdout, "%c %d\n", '%', n_threads );


  fprintf( stdout, "\nclear all;\n\n" );


  if     ( m_input >  0 ) {
    sprintf( m_dim_desc, "m = %d", m_input );
    sprintf( m_dim_tag,  "m%dc", m_input);
  }
  else if( m_input <  -1 ) {
    sprintf( m_dim_desc, "m = p/%d", -m_input );
    sprintf( m_dim_tag,  "m%dp", -m_input );
  }
  else if( m_input == -1 ) {
    sprintf( m_dim_desc, "m = p" );
    sprintf( m_dim_tag,  "m%dp", 1 );
  }


  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  FLASH_Queue_set_num_threads( n_threads );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;

    if( m < 0 ) m = p / abs(m_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){

      FLASH_Obj_create( datatype, m, m, 1, &nb_alg, &A );
      FLASH_Obj_create( datatype, m, m, 1, &nb_alg, &A_ref );

      if ( pc_str[param_combo][0] == 'l' )
        FLASH_Random_spd_matrix( FLA_LOWER_TRIANGULAR, A );
      else
        FLASH_Random_spd_matrix( FLA_UPPER_TRIANGULAR, A );
      
      FLASH_Copy( A, A_ref );

      fprintf( stdout, "data_chol_%s( %d, 1:5 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );


      time_Chol( param_combo, FLA_ALG_REFERENCE, n_repeats, m,
                 A, A_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Chol( param_combo, FLA_ALG_FRONT, n_repeats, m,
                 A, A_ref, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &A_ref );
    }

    fprintf( stdout, "\n" );
  }


  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_param_combos; i++ ) {
    fprintf( stdout, "plot( data_chol_%s( :,1 ), data_chol_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_chol_%s( :,1 ), data_chol_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_param_combos; i++ )
    fprintf( stdout, "'ref\\_chol\\_%s', 'fla\\_chol\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME chol front-end performance (%s)' );\n", m_dim_desc );
  fprintf( stdout, "print -depsc chol_front_%s.eps\n", m_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );

  FLA_Finalize( );

  return 0;
}

