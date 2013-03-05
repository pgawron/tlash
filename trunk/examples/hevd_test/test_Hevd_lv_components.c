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

void fill_eigenvalues( FLA_Obj l );

void time_Hevd_lv_components(
               int variant, int type, int n_repeats, int m, int n_iter_max, int k_accum, int b_alg,
               FLA_Obj A, FLA_Obj l,
               double* dtime, double* diff1, double* diff2, double* gflops,
               double* dtime_tred, double* gflops_tred,
               double* dtime_tevd, double* gflops_tevd,
               double* dtime_appq, double* gflops_appq, int* k_perf );


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
    dist_type,
    variant,
    n_repeats,
    i,
    k_perf,
    first_var = 1,
    last_var  = 4;

  double
    dtime,
    gflops,
    diff1,
    diff2;
  double
    dtime_tred,
    dtime_tevd,
    dtime_appq,
    gflops_tred,
    gflops_tevd,
    gflops_appq;

  FLA_Datatype datatype, dt_real;

  FLA_Obj
    A, l, Q, Ql, alpha, shift, nc;
  

  FLA_Init();


  fprintf( stdout, "%c number of repeats:", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter n_iter_max:", '%' );
  scanf( "%d", &n_iter_max );
  fprintf( stdout, "%c %d\n", '%', n_iter_max );

  fprintf( stdout, "%c enter number of Givens rotations to accumulate:", '%' );
  scanf( "%d", &k_accum );
  fprintf( stdout, "%c %d\n", '%', k_accum );

  fprintf( stdout, "%c enter blocking size:", '%' );
  scanf( "%d", &b_alg );
  fprintf( stdout, "%c %d\n", '%', b_alg );

  fprintf( stdout, "%c enter problem size first, last, inc:", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m (-1 means bind to problem size): ", '%' );
  scanf( "%d", &m_input );
  fprintf( stdout, "%c %d\n", '%', m_input );

  fprintf( stdout, "%c enter distribution type (0=u, 1=i, 2=g, 3=l, 4=r, 5=c): ", '%' );
  scanf( "%d", &dist_type );
  fprintf( stdout, "%c %d\n", '%', dist_type );


  fprintf( stdout, "\n" );



  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {

    m = m_input;

    if( m < 0 ) m = p / abs(m_input);

    //datatype = FLA_FLOAT;
    //datatype = FLA_DOUBLE;
    //datatype = FLA_COMPLEX;
    datatype = FLA_DOUBLE_COMPLEX;

    FLA_Obj_create( datatype, m,  m, 0, 0, &A );
    FLA_Obj_create( datatype, m,  m, 0, 0, &Q );
    FLA_Obj_create( datatype, m,  m, 0, 0, &Ql );

    dt_real = FLA_Obj_datatype_proj_to_real( A );

    FLA_Obj_create( dt_real,  m,  1, 0, 0, &l );
    FLA_Obj_create( dt_real,  1,  1, 0, 0, &alpha );
    FLA_Obj_create( dt_real,  1,  1, 0, 0, &shift );
    FLA_Obj_create( FLA_INT,  1,  1, 0, 0, &nc );

    FLA_Random_unitary_matrix( Q );

    if ( dist_type == 0 )
    {
      // Linear
      *FLA_DOUBLE_PTR( shift ) = 0.0;
      *FLA_DOUBLE_PTR( alpha ) = 1.0;
      fprintf( stdout, "%c using linear dist.\n", '%' );
      fprintf( stdout, "%c delta = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_linear_dist( shift, alpha, l );
    }
    else if ( dist_type == 1 )
    {
      // Inverse
      *FLA_DOUBLE_PTR( alpha ) = 1.0;
      fprintf( stdout, "%c using inverse dist.\n", '%' );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_inverse_dist( alpha, l );
    }
    else if ( dist_type == 2 )
    {
      // Geometric
      *FLA_DOUBLE_PTR( alpha ) = 1.0 / (double)m;
      fprintf( stdout, "%c using geometric dist.\n", '%' );
      fprintf( stdout, "%c alpha = %10.4e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_geometric_dist( alpha, l );
    }
    else if ( dist_type == 3 )
    {
      // Logarithmic
      *FLA_DOUBLE_PTR( alpha ) = 10;
      fprintf( stdout, "%c using logarithmic dist.\n", '%' );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_logarithmic_dist( alpha, l );
    }
    else if ( dist_type == 4 )
    {
      // Random
      *FLA_DOUBLE_PTR( shift ) = 0.0;
      *FLA_DOUBLE_PTR( alpha ) = (double)m;
      fprintf( stdout, "%c using random dist.\n", '%' );
      fprintf( stdout, "%c shift = %13.8e\n", '%', *FLA_DOUBLE_PTR( shift ) );
      fprintf( stdout, "%c alpha = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_random_dist( shift, alpha, l );
    }
    else if ( dist_type == 5 )
    {
      // Cluster
      //*FLA_INT_PTR( nc )       = 10;
      //*FLA_DOUBLE_PTR( alpha ) = 1.0e-9;
      *FLA_INT_PTR( nc )       = 10;
      *FLA_DOUBLE_PTR( alpha ) = 1.0e-9;
      fprintf( stdout, "%c using cluster dist.\n", '%' );
      fprintf( stdout, "%c num clusters  = %d\n", '%', *FLA_INT_PTR( nc ) );
      fprintf( stdout, "%c cluster width = %9.3e\n", '%', *FLA_DOUBLE_PTR( alpha ) );
      FLA_Fill_with_cluster_dist( nc, alpha, l );
    }

    {
      FLA_Copy( Q, Ql );
      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, Ql );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                FLA_ONE, Ql, Q, FLA_ZERO, A );
      FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
    }

    FLA_Set( FLA_ZERO, l );

    fprintf( stdout, "%c                               total----------                       reduction------    tri evd--------   form/app Q-----\n", '%' );
    fprintf( stdout, "%c                         %4s  gflops dtime      resid    |I-QQ'|    gflops dtime       gflops dtime      gflops dtime    niter\n", '%', "p" );

    time_Hevd_lv_components( -3, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_refq( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );

    time_Hevd_lv_components( -4, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_refd( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );

    time_Hevd_lv_components( -5, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_refr( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );

    time_Hevd_lv_components( 0, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_REFq( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );

    time_Hevd_lv_components( -1, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_REFd( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );

    time_Hevd_lv_components( -2, FLA_ALG_REFERENCE, n_repeats,
                             m, n_iter_max, k_accum, b_alg,
                             A, l, &dtime, &diff1, &diff2, &gflops,
                             &dtime_tred, &gflops_tred,
                             &dtime_tevd, &gflops_tevd,
                             &dtime_appq, &gflops_appq, &k_perf );
    if ( dtime_tred == 1.0 ) gflops_tred = gflops_tevd = gflops_appq =
                             dtime_tred  =  dtime_tevd =  dtime_appq = 0.0;

    fprintf( stdout, "data_REFr( %2d, 1:12 ) = [ %4d %6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ];\n",
             i, p, gflops, dtime, diff1, diff2,
             gflops_tred, dtime_tred,
             gflops_tevd, dtime_tevd,
             gflops_appq, dtime_appq, k_perf );
    fflush( stdout );


    for ( variant = first_var; variant <= last_var; variant++ ){
      
      fprintf( stdout, "data_var%d( %2d, 1:12 ) = [ %4d ", variant, i, p );
      fflush( stdout );

      time_Hevd_lv_components( variant, FLA_ALG_UNBLOCKED, n_repeats,
                               m, n_iter_max, k_accum, b_alg,
                               A, l, &dtime, &diff1, &diff2, &gflops,
                               &dtime_tred, &gflops_tred,
                               &dtime_tevd, &gflops_tevd,
                               &dtime_appq, &gflops_appq, &k_perf );

      fprintf( stdout, "%6.3lf %9.2e   %6.2le %6.2le  %6.3lf %9.2e  %7.3lf %9.2e  %6.3lf %9.2e %5d ",
               gflops, dtime, diff1, diff2,
               gflops_tred, dtime_tred,
               gflops_tevd, dtime_tevd,
               gflops_appq, dtime_appq, k_perf );

      fprintf( stdout, "];\n" );
      fflush( stdout );
    }

    fprintf( stdout, "\n" );

    FLA_Obj_free( &A );
    FLA_Obj_free( &l );
    FLA_Obj_free( &Q );
    FLA_Obj_free( &Ql );
    FLA_Obj_free( &alpha );
    FLA_Obj_free( &shift );
    FLA_Obj_free( &nc );
  }

  FLA_Finalize( );

  return 0;
}

