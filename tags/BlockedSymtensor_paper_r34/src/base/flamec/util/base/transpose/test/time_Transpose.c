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
#define FLA_ALG_BLOCKED   1
#define FLA_ALG_UNBLOCKED 2


void time_Transpose(
                  int variant, int type, int nrepeats, int n, int nb_alg,
                  FLA_Obj A, FLA_Obj A_ref,
                  double *dtime, double *diff, double *gflops );


void time_Transpose(
                  int variant, int type, int nrepeats, int n, int nb_alg,
                  FLA_Obj A, FLA_Obj A_ref,
                  double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9;

  FLA_Obj
    A_old, A_tmp;

  fla_blocksize_t*
    bp;
  fla_transpose_t*
    cntl_trans_var_unb;
  fla_transpose_t*
    cntl_trans_var_blk;
  fla_swap_t*
    cntl_swap_var_blk;
  fla_swap_t*
    cntl_swap_blas;


  bp                 = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_swap_blas     = FLA_Cntl_swap_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_swap_var_blk  = FLA_Cntl_swap_obj_create( FLA_FLAT, FLA_UNBLOCKED_VARIANT1, bp, cntl_swap_blas );
  cntl_trans_var_unb = FLA_Cntl_transpose_obj_create( FLA_FLAT, FLA_UNBLOCKED_VARIANT1, NULL, NULL, NULL );
  cntl_trans_var_blk = FLA_Cntl_transpose_obj_create( FLA_FLAT, variant, bp, cntl_trans_var_unb, cntl_swap_var_blk );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_old );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &A_tmp );

  FLA_Copy_external( A, A_old );


  for ( irep = 0 ; irep < nrepeats; irep++ ){
    FLA_Copy_external( A_old, A );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:

      //FLA_Copyt_external( FLA_TRANSPOSE, A, A_tmp );
      //FLA_Set( FLA_ZERO, A );
      //FLA_Copyt_external( FLA_NO_TRANSPOSE, A_tmp, A );
      FLA_Transpose( A );

      break;

    case 1:{

      /* Time variant 1 */
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Transpose_unb_var1( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Transpose_blk_var1( A, cntl_trans_var_blk );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{

      /* Time variant 2 */
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Transpose_unb_var2( A );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Transpose_blk_var2( A, cntl_trans_var_blk );
        break;
      default:
        printf("trouble\n");
      }

      break;
    } 

    }

    *dtime = FLA_Clock() - *dtime;
    dtime_old = min( *dtime, dtime_old );
  }

  FLA_Cntl_obj_free( cntl_trans_var_blk );
  FLA_Cntl_obj_free( cntl_trans_var_unb );
  FLA_Cntl_obj_free( cntl_swap_var_blk );
  FLA_Cntl_obj_free( cntl_swap_blas );
  FLA_Blocksize_free( bp );

  if ( variant == 0 ){
    FLA_Copy_external( A, A_ref );
    *diff = 0.0;
  }
  else{
    *diff = FLA_Max_elemwise_diff( A, A_ref );
  }

  *gflops = 4 * n * n /
            dtime_old / 1e9;

  *dtime = dtime_old;

  FLA_Copy_external( A_old, A );

  FLA_Obj_free( &A_old );
  FLA_Obj_free( &A_tmp );
}

