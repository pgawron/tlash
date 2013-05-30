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
#define FLA_ALG_UNBLOCKED 1
#define FLA_ALG_BLOCKED   2


void time_Symm_ll(
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops );


void time_Symm_ll( 
               int variant, int type, int nrepeats, int n, int nb_alg,
               FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj C_ref,
               double *dtime, double *diff, double *gflops )
{
  int
    irep;

  double
    dtime_old = 1.0e9; 

  FLA_Obj
    C_old;

  fla_blocksize_t*
    bp;
  fla_gemm_t*
    cntl_gemm_blas;
  fla_symm_t*
    cntl_symm_blas;
  fla_symm_t*
    cntl_symm_var;

  if ( type == FLA_ALG_UNBLOCKED && n > 300 )
  {
    *diff = 0.0;
    *gflops = 0.0;
    return;
  }

  bp             = FLA_Blocksize_create( nb_alg, nb_alg, nb_alg, nb_alg );
  cntl_gemm_blas = FLA_Cntl_gemm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL );
  cntl_symm_blas = FLA_Cntl_symm_obj_create( FLA_FLAT, FLA_SUBPROBLEM, NULL, NULL, NULL, NULL );
  cntl_symm_var  = FLA_Cntl_symm_obj_create( FLA_FLAT, variant, bp, cntl_symm_blas, cntl_gemm_blas, cntl_gemm_blas );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &C_old );

  FLA_Copy_external( C, C_old );


  for ( irep = 0 ; irep < nrepeats; irep++ )
  {
    FLA_Copy_external( C_old, C );

    *dtime = FLA_Clock();

    switch( variant ){
    case 0:
      // Time reference implementation
      REF_Symm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_ONE, A, B, FLA_ONE, C );
      break;

    case 1:{
      // Time variant 1
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var1( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var1( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 2:{
      // Time variant 2
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var2( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var2( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 3:{
      // Time variant 3
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var3( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var3( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 4:{
      // Time variant 4
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var4( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var4( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 5:{
      // Time variant 5
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var5( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var5( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 6:{
      // Time variant 6
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var6( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var6( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 7:{
      // Time variant 7
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var7( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var7( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 8:{
      // Time variant 8
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var8( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var8( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 9:{
      // Time variant 9
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var9( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var9( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
        break;
      default:
        printf("trouble\n");
      }

      break;
    }

    case 10:{
      // Time variant 10
      switch( type ){
      case FLA_ALG_UNBLOCKED:
        FLA_Symm_ll_unb_var10( FLA_ONE, A, B, FLA_ONE, C );
        break;
      case FLA_ALG_BLOCKED:
        FLA_Symm_ll_blk_var10( FLA_ONE, A, B, FLA_ONE, C, cntl_symm_var );
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

  FLA_Cntl_obj_free( cntl_symm_var );
  FLA_Cntl_obj_free( cntl_symm_blas );
  FLA_Cntl_obj_free( cntl_gemm_blas );
  FLA_Blocksize_free( bp );

  if ( variant == 0 )
  {
    FLA_Copy_external( C, C_ref );
    *diff = 0.0;
  }
  else
  {
    *diff = FLA_Max_elemwise_diff( C, C_ref );
  }

  *gflops = 2.0 * 
            FLA_Obj_length( C ) * 
            FLA_Obj_width( C ) * 
            FLA_Obj_width( A ) / 
            dtime_old / 
            1e9;

  *dtime = dtime_old;

  FLA_Copy_external( C_old, C );

  FLA_Obj_free( &C_old );
}

