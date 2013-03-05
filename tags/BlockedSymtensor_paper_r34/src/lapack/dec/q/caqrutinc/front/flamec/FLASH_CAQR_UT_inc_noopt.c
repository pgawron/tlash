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

extern fla_caqrutinc_t* flash_caqrutinc_cntl;

FLA_Error FLASH_CAQR_UT_inc_noopt( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW )
{
  FLA_Error r_val = FLA_SUCCESS;
  dim_t     nb_part;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_CAQR_UT_inc_check( p, A, ATW, R, RTW );

  // Compute the partition length from the number of partitions.
  nb_part = FLA_CAQR_UT_inc_compute_blocks_per_part( p, A );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Perform incremental QR's on each of the p partitions.
  FLA_CAQR_UT_inc_factorize_panels( nb_part, A, ATW );

  // Copy the triangles of A into R.
  FLA_CAQR_UT_inc_copy_triangles( nb_part, A, R );

  // Perform an incremental CAQR on the resulting upper triangular R's in A.
  FLA_CAQR_UT_inc_blk_var1( R, RTW, flash_caqrutinc_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

