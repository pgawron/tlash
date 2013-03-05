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


#define N_PARAM_COMBOS    1

#define FLA_ALG_FRONT     0
#define FLA_ALG_FRONT_ALT 1

char* pc_str[N_PARAM_COMBOS] = { "" };

void time_CAQR_UT_inc(
               int param_combo, int type, int nrepeats, int m, int n, dim_t n_panels,
               FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj b, FLA_Obj x,
               double *dtime, double *diff, double *gflops );


int main(int argc, char *argv[])
{
  int 
    datatype,
    n_threads,
    m_input,
    m,
    n_input,
    n,
    p_first, p_last, p_inc,
    p,
    n_repeats,
    param_combo,
    i,
    n_param_combos = N_PARAM_COMBOS;

  dim_t
    n_panels,
    nb_flash,
    nb_alg;

  double
    dtime,
    gflops,
    diff;

  FLA_Obj A, ATW, R, RTW, b, x;
  FLA_Obj A_flat, b_flat, x_flat;
  
  FLA_Init( );


  fprintf( stdout, "%c number of repeats: ", '%' );
  scanf( "%d", &n_repeats );
  fprintf( stdout, "%c %d\n", '%', n_repeats );

  fprintf( stdout, "%c enter algorithmic blocksize: ", '%' );
  scanf( "%u", &nb_alg );
  fprintf( stdout, "%c %u\n", '%', nb_alg );

  fprintf( stdout, "%c enter FLASH blocksize: ", '%' );
  scanf( "%u", &nb_flash );
  fprintf( stdout, "%c %u\n", '%', nb_flash );

  fprintf( stdout, "%c enter problem size first, last, inc: ", '%' );
  scanf( "%d%d%d", &p_first, &p_last, &p_inc );
  fprintf( stdout, "%c %d %d %d\n", '%', p_first, p_last, p_inc );

  fprintf( stdout, "%c enter m n (-1 means bind to problem size): ", '%' );
  scanf( "%d %d", &m_input, &n_input );
  fprintf( stdout, "%c %d %d\n", '%', m_input, n_input );

  fprintf( stdout, "%c enter the number of QR subproblem panels: ", '%' );
  scanf( "%u", &n_panels );
  fprintf( stdout, "%c %u\n", '%', n_panels );

  fprintf( stdout, "%c enter the number of SuperMatrix threads: ", '%' );
  scanf( "%d", &n_threads );
  fprintf( stdout, "%c %d\n", '%', n_threads );



  //datatype = FLA_FLOAT;
  datatype = FLA_DOUBLE;
  //datatype = FLA_COMPLEX;
  //datatype = FLA_DOUBLE_COMPLEX;

  //FLASH_Queue_disable();
  FLASH_Queue_set_num_threads( n_threads );
  //FLASH_Queue_set_verbose_output( TRUE );
  // FLA_Check_error_level_set( FLA_NO_ERROR_CHECKING );

  for ( p = p_first, i = 1; p <= p_last; p += p_inc, i += 1 )
  {
    m = m_input;
    n = n_input;

    if ( m < 0 ) m = p * abs(m_input);
    if ( n < 0 ) n = p * abs(n_input);

    for ( param_combo = 0; param_combo < n_param_combos; param_combo++ )
    {
      FLA_Obj_create( datatype, m, n, 0, 0, &A_flat );
      FLA_Obj_create( datatype, n, 1, 0, 0, &x_flat );
      FLA_Obj_create( datatype, m, 1, 0, 0, &b_flat );

      FLA_Random_matrix( A_flat );
      FLA_Random_matrix( b_flat );

      FLASH_CAQR_UT_inc_create_hier_matrices( n_panels, A_flat, 1, &nb_flash, nb_alg,
                                              &A, &ATW, &R, &RTW );
      FLASH_Obj_create_hier_copy_of_flat( b_flat, 1, &nb_flash, &b );
      FLASH_Obj_create_hier_copy_of_flat( x_flat, 1, &nb_flash, &x );


      fprintf( stdout, "data_caqrutinc_%s( %d, 1:3 ) = [ %d  ", pc_str[param_combo], i, p );
      fflush( stdout );

      time_CAQR_UT_inc( param_combo, FLA_ALG_FRONT, n_repeats, m, n, n_panels,
                        A, ATW, R, RTW, b, x, &dtime, &diff, &gflops );

      fprintf( stdout, "%6.3lf %6.2le ", gflops, diff );
      fflush( stdout );


      fprintf( stdout, " ]; \n" );
      fflush( stdout );

      FLA_Obj_free( &A_flat );
      FLA_Obj_free( &b_flat );
      FLA_Obj_free( &x_flat );

      FLASH_Obj_free( &A );
      FLASH_Obj_free( &ATW );
      FLASH_Obj_free( &R );
      FLASH_Obj_free( &RTW );
      FLASH_Obj_free( &b );
      FLASH_Obj_free( &x );
    }

  }


  FLA_Finalize( );

  return 0;
}

