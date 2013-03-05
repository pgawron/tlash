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

extern fla_sylv_t* flash_sylv_cntl;
extern fla_sylv_t* fla_sylv_cntl_leaf;

FLA_Error FLA_Sylv_internal( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Sylv_internal_check( transa, transb, isgn, A, B, C, scale, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Sylv_internal( transa,
		                           transb,
		                           isgn,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           *FLASH_OBJ_PTR_AT( C ),
		                           scale,
		                           flash_sylv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Sylv( transa, transb, isgn, A, B, C, scale, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_sylv_cntl_leaf;
		}

		// Parameter combinations
		if      ( transa == FLA_NO_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Sylv_nn( isgn, A, B, C, scale, cntl );
			else if ( transb == FLA_TRANSPOSE || transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Sylv_nh( isgn, A, B, C, scale, cntl );
		}
		else if ( transa == FLA_TRANSPOSE || transa == FLA_CONJ_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Sylv_hn( isgn, A, B, C, scale, cntl );
			else if ( transb == FLA_TRANSPOSE || transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Sylv_hh( isgn, A, B, C, scale, cntl );
		}
	}

	return r_val;
}

