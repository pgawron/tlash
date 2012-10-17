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

extern fla_caqr2ut_t* flash_caqr2ut_cntl;
extern fla_caqr2ut_t* fla_caqr2ut_cntl_leaf;

FLA_Error FLA_CAQR2_UT_internal( FLA_Obj B,
                                 FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_CAQR2_UT_internal_check( B, D, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( B ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_CAQR2_UT_internal( *FLASH_OBJ_PTR_AT( B ),
		                              *FLASH_OBJ_PTR_AT( D ),
		                              *FLASH_OBJ_PTR_AT( T ),
		                              flash_caqr2ut_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( B ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		if      ( FLA_Obj_structure( D ) == FLA_FULL_MATRIX )
		{
			ENQUEUE_FLASH_QR2_UT( B, D, T, cntl );
		}
		else if ( FLA_Obj_structure( D ) == FLA_UPPER_TRIANGULAR )
		{
			ENQUEUE_FLASH_CAQR2_UT( B, D, T, cntl );
		}
		else if ( FLA_Obj_structure( D ) == FLA_ZERO_MATRIX )
		{
			// Don't enqueue any tasks for zero blocks.
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( B ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf.
			if      ( FLA_Obj_structure( D ) == FLA_FULL_MATRIX )
			{
				FLA_QR2_UT_task( B, D, T, NULL );
				return FLA_SUCCESS;
			}
			else if ( FLA_Obj_structure( D ) == FLA_UPPER_TRIANGULAR )
				cntl = fla_caqr2ut_cntl_leaf;
			else if ( FLA_Obj_structure( D ) == FLA_ZERO_MATRIX )
				return FLA_SUCCESS;
			else
				FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}

		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_unb_var1( B, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_opt_var1( B, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_blk_var1( B, D, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_CAQR2_UT_blk_var2( B, D, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

