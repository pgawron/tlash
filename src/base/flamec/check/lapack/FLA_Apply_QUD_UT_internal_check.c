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

FLA_Error FLA_Apply_QUD_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( R, T );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, W );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, U );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, C );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, V );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, D );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects.
	if ( side == FLA_LEFT )
	{
		if ( FLA_Obj_elemtype( R ) == FLA_MATRIX )
		{
			e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( U ) );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_object_length_equals( T, max( FLA_Obj_length( U ),
			                                                FLA_Obj_length( V ) ) );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, W, R );
			FLA_Check_error_code( e_val );
		}
		else // FLA_Obj_elemtype( R ) == FLA_SCALAR
		{
		}

		e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, U, R, C );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, V, R, D );
		FLA_Check_error_code( e_val );

	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return FLA_SUCCESS;
}

