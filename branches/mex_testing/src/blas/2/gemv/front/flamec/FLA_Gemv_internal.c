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

extern fla_gemv_t* flash_gemv_cntl_blas;
extern fla_gemv_t* flash_gemv_cntl_fm_rp;

FLA_Error FLA_Gemv_internal( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Gemv_internal_check( transa, alpha, A, x, beta, y, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Gemv_internal( transa, 
		                           alpha, 
		                           *FLASH_OBJ_PTR_AT( A ), 
		                           *FLASH_OBJ_PTR_AT( x ), 
		                           beta, 
		                           *FLASH_OBJ_PTR_AT( y ), 
		                           flash_gemv_cntl_fm_rp );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Gemv( transa, alpha, A, x, beta, y, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_gemv_cntl_blas;
		}

		// Parameter combinations
		if      ( transa == FLA_NO_TRANSPOSE )
		{
			r_val = FLA_Gemv_n( alpha, A, x, beta, y, cntl );
		}
		else if ( transa == FLA_TRANSPOSE )
		{
			r_val = FLA_Gemv_t( alpha, A, x, beta, y, cntl );
		}
		else if ( transa == FLA_CONJ_TRANSPOSE )
		{
			r_val = FLA_Gemv_h( alpha, A, x, beta, y, cntl );
		}
	}

	return r_val;
}

