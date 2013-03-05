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

extern fla_trmm_t* flash_trmm_cntl_blas;
extern fla_trmm_t* flash_trmm_cntl_mm;

FLA_Error FLA_Trmm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Trmm_internal_check( side, uplo, transa, diag, alpha, A, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Trmm_internal( side,
		                           uplo,
		                           transa,
		                           diag,
		                           alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           flash_trmm_cntl_mm );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Trmm( side, uplo, transa, diag, alpha, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_trmm_cntl_blas;
		}

		// Parameter combinations
		if      ( side == FLA_LEFT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trmm_lln( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trmm_llt( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trmm_llc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trmm_llh( diag, alpha, A, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trmm_lun( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trmm_lut( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trmm_luc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trmm_luh( diag, alpha, A, B, cntl );
			}
		}
		else if ( side == FLA_RIGHT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trmm_rln( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trmm_rlt( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trmm_rlc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trmm_rlh( diag, alpha, A, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trmm_run( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trmm_rut( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trmm_ruc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trmm_ruh( diag, alpha, A, B, cntl );
			}
		}
	}

	return r_val;
}

