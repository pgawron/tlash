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

void time_Bidiag_UT(
                int param_combo, int type, int n_repeats, int m, int n,
                FLA_Obj A, FLA_Obj tu, FLA_Obj tv, FLA_Obj TU, FLA_Obj TV,
                double *dtime, double *diff, double *gflops );

int main(int argc, char *argv[])
{
  int 
    m_input, n_input,
    m, n,
    min_m_n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    //n_param_combos = N_PARAM_COMBOS;
    n_param_combos = 1;

  FLA_Datatype datatype;

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  dim_t b_alg;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj
    A, tu, tv, TU, TV;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c Enter blocking size:", '%' );
  scanf( "%u", &b_alg );
  fprintf( stdout, "%c %u\n", '%', b_alg );

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


  //datatype = FLA_FLOAT;
  //datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  datatype = FLA_DOUBLE_COMPLEX;

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if ( m < 0 ) m = p / abs(m_input);
    if ( n < 0 ) n = p / abs(n_input);

    min_m_n = min( m, n );

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ ){
      
      FLA_Obj_create( datatype, m,       n, 0, 0, &A );

      FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &tu );
      FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &tv );

      if ( b_alg == 0 )
      {
        FLA_Bidiag_UT_create_T( A, &TU, &TV );
      }
      else
      {
        FLA_Obj_create( datatype, b_alg, min_m_n, 0, 0, &TU );
        FLA_Obj_create( datatype, b_alg, min_m_n, 0, 0, &TV );
      }

      FLA_Random_matrix( A );

      fprintf( stdout, "data_bidiag_ut( %d, 1:3 ) = [ %d  ", i, p );
      fflush( stdout );

      time_Bidiag_UT( param_combo, FLA_ALG_REFERENCE, n_repeats, m, n,
                       A, tu, tv, TU, TV, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      time_Bidiag_UT( param_combo, FLA_ALG_FRONT, n_repeats, m, n,
                       A, tu, tv, TU, TV, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );

      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A );
      FLA_Obj_free( &tu );
      FLA_Obj_free( &tv );
      FLA_Obj_free( &TU );
      FLA_Obj_free( &TV );
    }

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
  fprintf( stdout, "title( 'FLAME Bidiag_UT front-end performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc bidiag_ut_front_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

