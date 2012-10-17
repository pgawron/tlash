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

FLA_Error FLA_Form_perm_matrix( FLA_Obj p, FLA_Obj A )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Form_perm_matrix_check( p, A );

  // We assume that A is correctly sized, m x m, where m is the row
  // dimension of the matrix given to FLA_LU_piv() or similar function.
  FLA_Set_to_identity( A );

  // We assume that p contains pivots in native FLAME format. That is,
  // we assume the pivot type is FLA_NATIVE_PIVOTS. This is not a huge
  // assumption since the user has to go out of his way to shift the
  // pivots into LAPACK-indexed pivots.
  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

  return FLA_SUCCESS;
}

