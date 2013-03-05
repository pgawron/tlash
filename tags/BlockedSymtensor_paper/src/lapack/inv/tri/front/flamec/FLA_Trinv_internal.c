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

extern fla_trinv_t* flash_trinv_cntl;
extern fla_trinv_t* fla_trinv_cntl_leaf;

FLA_Error FLA_Trinv_internal( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Trinv_internal_check( uplo, diag, A, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Trinv_internal( uplo,
		                            diag,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            flash_trinv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Trinv( uplo, diag, A, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_trinv_cntl_leaf;
		}

		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			if      ( diag == FLA_NONUNIT_DIAG )
			{
				r_val = FLA_Trinv_ln( A, cntl );
			}
			else if ( diag == FLA_UNIT_DIAG )
			{
				r_val = FLA_Trinv_lu( A, cntl );
			}
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			if      ( diag == FLA_NONUNIT_DIAG )
			{
				r_val = FLA_Trinv_un( A, cntl );
			}
			else if ( diag == FLA_UNIT_DIAG )
			{
				r_val = FLA_Trinv_uu( A, cntl );
			}
		}
	}

	return r_val;
}

