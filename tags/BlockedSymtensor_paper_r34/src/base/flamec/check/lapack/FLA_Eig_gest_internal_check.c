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

FLA_Error FLA_Eig_gest_internal_check( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, Y );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( A, B );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
	FLA_Check_error_code( e_val );

	if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 ||
	     FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		e_val = FLA_Check_object_width_equals( Y, FLA_Obj_width( A ) );
		FLA_Check_error_code( e_val );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 ||
	          FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		e_val = FLA_Check_object_length_equals( Y, FLA_Obj_length( A ) );
		FLA_Check_error_code( e_val );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, Y, A );
		FLA_Check_error_code( e_val );
	}

	return FLA_SUCCESS;
}

