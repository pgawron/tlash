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

FLA_Error FLA_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Obj R;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Random_spd_matrix_check( uplo, A );

  // Create a temporary object R conformal to A.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &R );

  // Randomize R to be an uplo-triangular matrix. Note that the diagonal of R
  // needs to be positive to ensure that R * R' is SPD/HPD.
  FLA_Random_tri_matrix( uplo, FLA_NONUNIT_DIAG, R );
    
  if ( uplo == FLA_LOWER_TRIANGULAR )
  {
    // A = R * R';
    FLA_Herk_external( uplo, FLA_NO_TRANSPOSE, FLA_ONE, R, FLA_ZERO, A );
  }
  else // if ( uplo == FLA_UPPER_TRIANGULAR )
  {
    // A = R' * R;
    FLA_Herk_external( uplo, FLA_CONJ_TRANSPOSE, FLA_ONE, R, FLA_ZERO, A );
  }

  // Free R.
  FLA_Obj_free( &R );

  return FLA_SUCCESS;
}

