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

FLA_Error REF_Tevd_v( FLA_Obj d, FLA_Obj e, FLA_Obj U );
void time_Tevd_v(
               int variant, int type, int n_repeats, int m, int k_accum, int b_alg, int n_iter_max,
               FLA_Obj A_orig, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops );


void time_Tevd_v(
               int variant, int type, int n_repeats, int m, int k_accum, int b_alg, int n_iter_max,
               FLA_Obj A_orig, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj A, FLA_Obj l,
               double *dtime, double *diff1, double* diff2, double *gflops )
{
  int irep;

  double
    k, dtime_old = 1.0e9;

  FLA_Obj
    A_save, G_save, d_save, e_save;

  if (
       //( variant == 0 ) ||
       //( variant == 1 && type == FLA_ALG_UNB_OPT ) ||
       //( variant == 2 && type == FLA_ALG_UNB_OPT ) ||
       FALSE
     )
  {
    *dtime  = 0.0;
    *gflops = 0.0;
    *diff1  = 0.0;
    *diff2  = 0.0;
    return;
  }

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, G, &G_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, d, &d_save );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, e, &e_save );

  FLA_Copy_external( A, A_save );
  FLA_Copy_external( G, G_save );
  FLA_Copy_external( d, d_save );
  FLA_Copy_external( e, e_save );

  for ( irep = 0 ; irep < n_repeats; irep++ ){

    FLA_Copy_external( A_save, A );
    FLA_Copy_external( G_save, G );
    FLA_Copy_external( d_save, d );
    FLA_Copy_external( e_save, e );

    *dtime = FLA_Clock();

    switch( variant ){

    case 0:
      REF_Tevd_v( d, e, A );
      break;

    // Time variant 1
    case 1:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Tevd_v_opt_var1( n_iter_max, d, e, G, A, b_alg );
        break;
      }
      break;
    }

    // Time variant 2
    case 2:
    {
      switch( type ){
      case FLA_ALG_UNB_OPT:
        FLA_Tevd_v_opt_var2( n_iter_max, d, e, G, R, W, A, b_alg );
        break;
      }
      break;
    }

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );

  }
  {
    FLA_Obj V, A_rev_evd, norm, eye;

	FLA_Copy( d, l );

//FLA_Obj_show( "A_save", A_save, "%9.2e + %9.2e ", "" );
//FLA_Obj_show( "A_evd", A, "%9.2e + %9.2e ", "" );
	FLA_Sort_evd( FLA_FORWARD, l, A );

    FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &V ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_rev_evd ); 
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &eye ); 
    FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );


    FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, l, A );

    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, A, V, FLA_ZERO, A_rev_evd );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A_rev_evd );

/*
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, A, D, FLA_ZERO, A_rev_evd );
    FLA_Copy( A_rev_evd, D );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, D, V, FLA_ZERO, A_rev_evd );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A_rev_evd );
*/
//FLA_Obj_show( "A_rev_evd", A_rev_evd, "%9.2e + %9.2e ", "" );
 
    FLA_Axpy( FLA_MINUS_ONE, A_orig, A_rev_evd );
    FLA_Norm_frob( A_rev_evd, norm );
    FLA_Obj_extract_real_scalar( norm, diff1 );
    //*diff = FLA_Max_elemwise_diff( A_orig, A_rev_evd );

    FLA_Set_to_identity( eye );
	FLA_Copy( V, A_rev_evd );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
              FLA_ONE, V, A_rev_evd, FLA_MINUS_ONE, eye );
    FLA_Norm_frob( eye, norm );
    FLA_Obj_extract_real_scalar( norm, diff2 );

/*
FLA_Obj_free( &EL );
FLA_Obj_free( &EU );
FLA_Obj_free( &D );
FLA_Obj_free( &dc );
FLA_Obj_free( &ec );
*/

    FLA_Obj_free( &V );
    FLA_Obj_free( &A_rev_evd );
    FLA_Obj_free( &eye );
    FLA_Obj_free( &norm );
  }

  k = 2.00;

  if ( FLA_Obj_is_complex( A ) )
  {
    *gflops = (
                      (       4.5 * k * m * m     ) +
                2.0 * (       3.0 * k * m * m * m ) ) / 
              dtime_old / 1e9;
  }
  else 
  {
    *gflops = (
                      (       4.5 * k * m * m     ) +
                1.0 * (       3.0 * k * m * m * m ) ) / 
              dtime_old / 1e9;
  }

  *dtime = dtime_old;

  FLA_Copy_external( A_save, A );
  FLA_Copy_external( G_save, G );
  FLA_Copy_external( d_save, d );
  FLA_Copy_external( e_save, e );

  FLA_Obj_free( &A_save );
  FLA_Obj_free( &G_save );
  FLA_Obj_free( &d_save );
  FLA_Obj_free( &e_save );
}

