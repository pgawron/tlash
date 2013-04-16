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

FLA_Error FLASH_FS_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b )
{
  dim_t     nb_alg;
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_FS_incpiv_check( A, p, L, b );

  // *** The current forward substitution algorithm implemented assumes that 
  // the matrix has a hierarchical depth of 1. We check for that here, because
  // we anticipate that we'll use a more general algorithm in the future, and 
  // we don't want to forget to remove the constraint. ***
  if ( FLASH_Obj_depth( A ) != 1 )
  {
    FLA_Print_message( "FLASH_FS_incpiv() currently only supports matrices of depth 1",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }
  
  // Inspect the width of a the top-left element of L to get the algorithmic
  // blocksize we'll use throughout the LU_incpiv algorithm.
  nb_alg = FLASH_Obj_scalar_width_tl( L );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();
  
  // Execute tasks.
  r_val = FLASH_FS_incpiv_aux1( A, p, L, b, nb_alg );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();

  return r_val;
}
