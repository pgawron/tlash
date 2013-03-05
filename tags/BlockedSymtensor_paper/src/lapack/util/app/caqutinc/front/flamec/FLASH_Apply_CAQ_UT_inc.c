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

extern fla_apcaqutinc_t* flash_apcaqutinc_cntl;

FLA_Error FLASH_Apply_CAQ_UT_inc( dim_t p,
                                  FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                  FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj W, FLA_Obj B )
{
  FLA_Error r_val;
  dim_t     nb_part;
  FLA_Obj   WT, WB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_CAQ_UT_inc_check( side, trans, direct, storev, A, ATW, R, RTW, W, B );

  // Compute the partition length from the number of partitions.
  nb_part = FLA_CAQR_UT_inc_compute_blocks_per_part( p, R );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Apply the individual Q's from the incremental QR factorizations.
  FLA_Apply_CAQ_UT_inc_apply_panels( nb_part, A, ATW, W, B );

  FLA_Part_2x1( W,   &WT,
                     &WB,    1, FLA_TOP );

  // Apply the Q from the factorization of the upper triangular R's.
  r_val = FLA_Apply_CAQ_UT_inc_internal( side, trans, direct, storev,
                                         R, RTW, WT, B, flash_apcaqutinc_cntl );


  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

