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

#define N_PIVOT_COMBOS    1

#define FLA_ALG_REFERENCE 0
#define FLA_ALG_FRONT     1

char* pc_str[N_PIVOT_COMBOS] = { "nopiv" };

void time_LU(
              int is_pivoting, int type, int n_repeats, int m, int n,
              FLA_Obj C, FLA_Obj b, FLA_Obj b_orig, FLA_Obj b_norm, 
              double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    m_input, n_input,
    m, n, min_m_n,
    p_first, p_last, p_inc,
    pp,
    pivot_combo,
    n_repeats,
    i,
    n_pivot_combos = N_PIVOT_COMBOS;
  
  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char n_dim_desc[14];
  char m_dim_tag[10];
  char n_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    C, b, b_orig, b_norm;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );


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
  if     ( n_input >  0 ) {
    sprintf( n_dim_desc, "n = %d", n_input );
    sprintf( n_dim_tag,  "n%dc", n_input);
  }
  else if( n_input <  -1 ) {
    sprintf( n_dim_desc, "n = p/%d", -n_input );
    sprintf( n_dim_tag,  "n%dp", -n_input );
  }
  else if( n_input == -1 ) {
    sprintf( n_dim_desc, "n = p" );
    sprintf( n_dim_tag,  "n%dp", 1 );
  }

  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( pp = p_first, i = 1; pp <= p_last; pp += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if( m < 0 ) m = pp / abs(m_input);
    if( n < 0 ) n = pp / abs(n_input);

    min_m_n = min( m, n );

    for ( pivot_combo = 0; pivot_combo < n_pivot_combos; pivot_combo++ ){
      
      FLA_Obj_create( datatype, m, n, 0, 0, &C );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b_orig );

      if ( FLA_Obj_is_single_precision( C ) )
        FLA_Obj_create( FLA_FLOAT, 1, 1, 0, 0, &b_norm );
      else
        FLA_Obj_create( FLA_DOUBLE, 1, 1, 0, 0, &b_norm );

      FLA_Random_matrix( C );
      FLA_Random_matrix( b );
      FLA_Copy_external( b, b_orig );


      fprintf( stdout, "data_lu_%s( %d, 1:5 ) = [ %d  ", pc_str[pivot_combo], i, pp );
      fflush( stdout );

      //time_LU( pivot_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
      //         C, b, b_orig, b_norm, &dtime, &diff, &gflops );

      //fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      //fflush( stdout );

      time_LU( pivot_combo, FLA_ALG_FRONT, n_repeats, m, n,
               C, b, b_orig, b_norm, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &C );
      FLA_Obj_free( &b );
      FLA_Obj_free( &b_orig );
      FLA_Obj_free( &b_norm );
    }

    fprintf( stdout, "\n" );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 0; i < n_pivot_combos; i++ ) {
    fprintf( stdout, "plot( data_lu_%s( :,1 ), data_lu_%s( :, 2 ), '%c:%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
    fprintf( stdout, "plot( data_lu_%s( :,1 ), data_lu_%s( :, 4 ), '%c-.%c' ); \n",
            pc_str[i], pc_str[i], colors[ i ], ticks[ i ] );
  }

  fprintf( stdout, "legend( ... \n" );

  for ( i = 0; i < n_pivot_combos; i++ )
    fprintf( stdout, "'ref\\_lu\\_%s', 'fla\\_lu\\_%s', ... \n", pc_str[i], pc_str[i] );

  fprintf( stdout, "'Location', 'SouthEast' ); \n" );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME LU front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc lu_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

