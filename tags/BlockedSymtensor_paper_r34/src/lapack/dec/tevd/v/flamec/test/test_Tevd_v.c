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

#define FLA_ALG_REFERENCE     0
#define FLA_ALG_UNBLOCKED     1
#define FLA_ALG_UNB_OPT       2
#define FLA_ALG_BLOCKED       3


void time_Tevd_v(
               int variant, int type, int n_repeats, int m, int k_accum, int b_alg, int n_iter_max,
               FLA_Obj A_orig, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W2, FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops );


int main(int argc, char *argv[])
{
  int 
    m_input,
    m,
    p_first, p_last, p_inc,
    p,
    k_accum,
    b_alg,
    n_iter_max,
    variant,
    n_repeats,
    i,
    n_variants = 2;

  char *colors = "brkgmcbrkg";
  char *ticks  = "o+*xso+*xs";
  char m_dim_desc[14];
  char m_dim_tag[10];

  double max_gflops=6.0;

  double
    dtime,
    gflops,
    diff1, diff2;

  FLA_Datatype datatype, dt_real;

  FLA_Obj
    A, l, Q, Ql, TT, r, d, e, A_orig, G, R, W2, de, alpha;

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter n_iter_max (per eigenvalue): ", '%' );
  scanf( "%d", &n_iter_max );
  fprintf( stdout, "%c %d\n", '%', n_iter_max );

  fprintf( stdout, "%c enter number of sets of Givens rotations to accumulate:", '%' );
  scanf( "%d", &k_accum );
  fprintf( stdout, "%c %d\n", '%', k_accum );

  fprintf( stdout, "%c enter blocking size for application of G:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );


  fprintf( stdout, "\n" );


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


  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;

    if( m < 0 ) m = p / abs(m_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype,  m,         m, 0, 0, &A );
    FLA_Obj_create( datatype,  m,         m, 0, 0, &A_orig );
    FLA_Obj_create( datatype,  m,         m, 0, 0, &Q );
    FLA_Obj_create( datatype,  m,         m, 0, 0, &Ql );
    FLA_Obj_create( datatype,  m,         1, 0, 0, &r );
    FLA_Obj_create( datatype,  m,         m, 0, 0, &W2 );
    FLA_Obj_create( datatype,  m-1, k_accum, 0, 0, &G );

	dt_real = FLA_Obj_datatype_proj_to_real( A );

    FLA_Obj_create( dt_real, m,      1, 0, 0, &l );
    FLA_Obj_create( dt_real, m,      1, 0, 0, &d );
    FLA_Obj_create( dt_real, m-1,    1, 0, 0, &e );
    FLA_Obj_create( dt_real, m,      m, 0, 0, &R );

	FLA_Obj_create( dt_real, 1,      1, 0, 0, &alpha );

	*FLA_DOUBLE_PTR( alpha ) = 1.0 / ( sqrt( sqrt( (double) m ) ) );

	FLA_Random_unitary_matrix( Q );

	//FLA_Fill_with_uniform_dist( FLA_ONE,   l );
	//FLA_Fill_with_inverse_dist( FLA_ONE,   l );
	FLA_Fill_with_geometric_dist( alpha,   l );


    {
      FLA_Copy( Q, Ql );
      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, Ql );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, Ql, Q, FLA_ZERO, A );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
      FLA_Copy( A, A_orig );
    }

    FLA_Set( FLA_ZERO, l );
    FLA_Set( FLA_ZERO, Q );

	FLA_Tridiag_UT_create_T( A, &TT );
	FLA_Tridiag_UT( FLA_LOWER_TRIANGULAR, A, TT );
	FLA_Tridiag_UT_realify( FLA_LOWER_TRIANGULAR, A, r );
	FLA_Tridiag_UT_extract_diagonals( FLA_LOWER_TRIANGULAR, A, d, e );
	FLA_Tridiag_UT_form_Q( FLA_LOWER_TRIANGULAR, A, TT );
	FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, r, A );
    FLA_Obj_free( &TT );

    time_Tevd_v( 0, FLA_ALG_REFERENCE, n_repeats, m, k_accum, b_alg, n_iter_max,
                 A_orig, d, e, G, R, W2, A, l, &dtime, &diff1, &diff2, &gflops );

    fprintf( stdout, "data_REFq( %d, 1:3 ) = [ %d %6.3lf %9.2e %6.2le %6.2le ]; \n", i, p, gflops, dtime, diff1, diff2 );
    fflush( stdout );

    for ( variant = 1; variant <= n_variants; variant++ ){
      
      fprintf( stdout, "data_var%d( %d, 1:3 ) = [ %d ", variant, i, p );
      fflush( stdout );

      time_Tevd_v( variant, FLA_ALG_UNB_OPT, n_repeats, m, k_accum, b_alg, n_iter_max,
                   A_orig, d, e, G, R, W2, A, l, &dtime, &diff1, &diff2, &gflops );

      fprintf( stdout, "%6.3lf %9.2e %6.2le %6.2le ", gflops, dtime, diff1, diff2 );
      fflush( stdout );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &A_orig );
    FLA_Obj_free( &Q );
    FLA_Obj_free( &Ql );
    FLA_Obj_free( &G );
    FLA_Obj_free( &W2 );
    FLA_Obj_free( &r );
    FLA_Obj_free( &l );
    FLA_Obj_free( &d );
    FLA_Obj_free( &e );
    FLA_Obj_free( &R );
    FLA_Obj_free( &alpha );
  }

/*
  fprintf( stdout, "figure;\n" );

  fprintf( stdout, "plot( data_REF( :,1 ), data_REF( :, 2 ), '-' ); \n" );

  fprintf( stdout, "hold on;\n" );

  for ( i = 1; i <= n_variants; i++ ) {
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 2 ), '%c:%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
    fprintf( stdout, "plot( data_var%d( :,1 ), data_var%d( :, 4 ), '%c-.%c' ); \n",
            i, i, colors[ i-1 ], ticks[ i-1 ] );
  }

  fprintf( stdout, "legend( ... \n" );
  fprintf( stdout, "'Reference', ... \n" );

  for ( i = 1; i < n_variants; i++ )
    fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d', ... \n", i, i );
  fprintf( stdout, "'unb\\_var%d', 'blk\\_var%d' ); \n", i, i );

  fprintf( stdout, "xlabel( 'problem size p' );\n" );
  fprintf( stdout, "ylabel( 'GFLOPS/sec.' );\n" );
  fprintf( stdout, "axis( [ 0 %d 0 %.2f ] ); \n", p_last, max_gflops );
  fprintf( stdout, "title( 'FLAME Hevd_lv performance (%s, %s)' );\n", 
           m_dim_desc, n_dim_desc );
  fprintf( stdout, "print -depsc tridiag_%s_%s.eps\n", m_dim_tag, n_dim_tag );
  fprintf( stdout, "hold off;\n");
  fflush( stdout );
*/

  FLA_Finalize( );

  return 0;
}

