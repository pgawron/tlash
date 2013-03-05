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

extern fla_qrut_t* flash_qrut_cntl;
extern fla_qrut_t* flash_qrut_cntl_leaf;
extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_internal( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_QR_UT_internal_check( A, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		if ( FLASH_Queue_get_enabled( ) )
		{
			// Enqueue
			ENQUEUE_FLASH_QR_UT_macro( A, *FLASH_OBJ_PTR_AT( T ), cntl );
		}
		else
		{
			// Execute
			r_val = FLA_QR_UT_macro_task( A, *FLASH_OBJ_PTR_AT( T ), cntl );
		}
	}
	else
	{
		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_QR_UT_unb_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_QR_UT_opt_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_QR_UT_blk_var1( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
		{
			r_val = FLA_QR_UT_unb_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
		{
			r_val = FLA_QR_UT_opt_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_QR_UT_blk_var2( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_QR_UT_blk_var3( A, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

