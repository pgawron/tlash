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

extern fla_hemm_t* flash_hemm_cntl_blas;
extern fla_hemm_t* flash_hemm_cntl_mm;

FLA_Error FLA_Hemm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_hemm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Hemm_internal_check( side, uplo, alpha, A, B, beta, C, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Hemm_internal( side,
		                           uplo,
		                           alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           beta,
		                           *FLASH_OBJ_PTR_AT( C ),
		                           flash_hemm_cntl_mm );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Hemm( side, uplo, alpha, A, B, beta, C, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_hemm_cntl_blas;
		}

		// Parameter combinations
		if      ( side == FLA_LEFT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
				r_val = FLA_Hemm_ll( alpha, A, B, beta, C, cntl );
			else if ( uplo == FLA_UPPER_TRIANGULAR )
				r_val = FLA_Hemm_lu( alpha, A, B, beta, C, cntl );
		}
		else if ( side == FLA_RIGHT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
				r_val = FLA_Hemm_rl( alpha, A, B, beta, C, cntl );
			else if ( uplo == FLA_UPPER_TRIANGULAR )
				r_val = FLA_Hemm_ru( alpha, A, B, beta, C, cntl );
		}
	}

	return r_val;
}

