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

FLA_Error FLASH_UDdate_UT_inc_solve( FLA_Obj R, FLA_Obj bR, FLA_Obj x )
{
	// Check parameters.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_UDdate_UT_inc_solve_check( R, bR, x );

	// Copy the contents of bR to x so that after the triangular solve, the
	// solution resides in x (and bR is preserved).
	FLASH_Copy( bR, x );
	
	// Perform a triangular solve with R the right-hand side.
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
	            FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	            FLA_ONE, R, x );

	return FLA_SUCCESS;
}

