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

extern fla_lu_t* fla_lu_nopiv_cntl;
extern fla_lu_t* fla_lu_nopiv_cntl2;

FLA_Error FLA_LU_nopiv( FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_nopiv_check( A );

  // Invoke FLA_LU_nopiv_internal() with large control tree.
  r_val = FLA_LU_nopiv_internal( A, fla_lu_nopiv_cntl2 );

  // Check for singularity.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    r_val = FLA_LU_find_zero_on_diagonal( A );

  return r_val;
}

