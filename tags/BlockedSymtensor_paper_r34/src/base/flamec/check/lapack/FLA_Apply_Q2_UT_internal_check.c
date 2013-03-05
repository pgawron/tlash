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

FLA_Error FLA_Apply_Q2_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apq2ut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( D, T );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, W );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, C );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, E );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects.
	if ( side == FLA_LEFT )
	{
		if ( FLA_Obj_elemtype( D ) == FLA_MATRIX )
		{
			e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, D, T );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, C, W );
			FLA_Check_error_code( e_val );
		}
		else // if ( FLA_Obj_elemtype( D ) == FLA_SCALAR )
		{
			//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, E, C, D );
			//FLA_Check_error_code( e_val );

			e_val = FLA_Check_object_width_equals( C, FLA_Obj_width( E ) );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_object_length_equals( D, FLA_Obj_length( E ) );
			FLA_Check_error_code( e_val );
		}

		//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, C, E );
		//FLA_Check_error_code( e_val );

	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return FLA_SUCCESS;
}

