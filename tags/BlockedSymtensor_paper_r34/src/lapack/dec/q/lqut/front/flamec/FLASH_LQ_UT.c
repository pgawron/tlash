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

extern fla_lqut_t*   flash_lqut_cntl;
extern fla_lqut_t*   fla_lqut_cntl_leaf;

FLA_Error FLASH_LQ_UT( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;
  dim_t     b_alg, b_flash;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LQ_UT_check( A, TW );

  // *** The current hierarchical LQ_UT algorithm assumes that the matrix
  // has a hierarchical depth of 1. We check for that here, because we
  //  anticipate that we'll use a more general algorithm in the future, and
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_LQ_UT() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Inspect the length of TTL to get the blocksize used by the LQ
  // factorization, which will be our inner blocksize for Apply_Q_UT.
  b_alg   = FLASH_Obj_scalar_length_tl( TW );
  b_flash = FLASH_Obj_scalar_width_tl( TW );

  // The traditional (non-incremental) LQ_UT algorithm-by-blocks requires
  // that the algorithmic blocksize be equal to the storage blocksize.
  if ( b_alg != b_flash )
  {
    FLA_Print_message( "FLASH_LQ_UT() requires that b_alg == b_store",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // The traditional (non-incremental) LQ_UT algorithm-by-blocks requires
  // that min_dim(A) % b_flash == 0.
  if ( FLASH_Obj_scalar_min_dim( A ) % b_flash != 0 )
  {
    FLA_Print_message( "FLASH_LQ_UT() requires that min_dim( A ) %% b_store == 0",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Invoke FLA_LQ_UT_internal() with hierarchical control tree.
  r_val = FLA_LQ_UT_internal( A, TW, flash_lqut_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

