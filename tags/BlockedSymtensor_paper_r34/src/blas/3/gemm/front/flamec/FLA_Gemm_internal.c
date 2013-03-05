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

extern fla_gemm_t* flash_gemm_cntl_blas;
extern fla_gemm_t* flash_gemm_cntl_mm_op;

FLA_Error FLA_Gemm_internal( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Gemm_internal_check( transa, transb, alpha, A, B, beta, C, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          (FLA_Obj_elemtype( A ) == FLA_TENSOR || FLA_Obj_elemtype( A ) == FLA_MATRIX) &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Gemm_internal( transa, 
		                           transb, 
		                           alpha, 
		                           *FLASH_OBJ_PTR_AT( A ), 
		                           *FLASH_OBJ_PTR_AT( B ), 
		                           beta, 
		                           *FLASH_OBJ_PTR_AT( C ), 
		                           flash_gemm_cntl_mm_op );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Gemm( transa, transb, alpha, A, B, beta, C, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_gemm_cntl_blas;
		}

		// Parameter combinations
		if      ( transa == FLA_NO_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Gemm_nn( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_TRANSPOSE )
				r_val = FLA_Gemm_nt( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_NO_TRANSPOSE )
				r_val = FLA_Gemm_nc( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Gemm_nh( alpha, A, B, beta, C, cntl );
		}
		else if ( transa == FLA_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Gemm_tn( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_TRANSPOSE )
				r_val = FLA_Gemm_tt( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_NO_TRANSPOSE )
				r_val = FLA_Gemm_tc( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Gemm_th( alpha, A, B, beta, C, cntl );
		}
		else if ( transa == FLA_CONJ_NO_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Gemm_cn( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_TRANSPOSE )
				r_val = FLA_Gemm_ct( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_NO_TRANSPOSE )
				r_val = FLA_Gemm_cc( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Gemm_ch( alpha, A, B, beta, C, cntl );
		}
		else if ( transa == FLA_CONJ_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Gemm_hn( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_TRANSPOSE )
				r_val = FLA_Gemm_ht( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_NO_TRANSPOSE )
				r_val = FLA_Gemm_hc( alpha, A, B, beta, C, cntl );
			else if ( transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Gemm_hh( alpha, A, B, beta, C, cntl );
		}
	}

	return r_val;
}

