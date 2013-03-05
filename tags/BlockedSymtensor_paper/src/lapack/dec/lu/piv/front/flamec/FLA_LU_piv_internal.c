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

extern fla_lu_t* fla_lu_piv_cntl_leaf;

FLA_Error FLA_LU_piv_internal( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		FLA_Error e_val = FLA_Check_null_pointer( ( void* ) cntl );
		FLA_Check_error_code( e_val );
	}

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_LU_piv_macro( A, *FLASH_OBJ_PTR_AT( p ), cntl );
	}
	else
	{
/*
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_lu_piv_cntl_leaf;
		}
*/
		// Choose implementation.
		if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_EXTERN ||
		          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
		{
			if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER )
			{
				// When performing an unblocked operation on a hierarhical
				// object, we assume it's because we need to perform the
				// the operation on a macro block.
				r_val = FLA_LU_piv_macro_task( A, *FLASH_OBJ_PTR_AT( p ), cntl );
			}
			else
			{
				r_val = FLA_LU_piv_blk_ext( A, p );
			}
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_EXTERN )
		{
			if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER )
			{
				// When performing an unblocked operation on a hierarhical
				// object, we assume it's because we need to perform the
				// the operation on a macro block.
				r_val = FLA_LU_piv_macro_task( A, *FLASH_OBJ_PTR_AT( p ), cntl );
			}
			else
			{
				r_val = FLA_LU_piv_unb_ext( A, p );
			}
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
		{
			r_val = FLA_LU_piv_unb_var3( A, p );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
		{
			r_val = FLA_LU_piv_unb_var4( A, p );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
		{
			r_val = FLA_LU_piv_unb_var5( A, p );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
		{
			r_val = FLA_LU_piv_opt_var3( A, p );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT4 )
		{
			r_val = FLA_LU_piv_opt_var4( A, p );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT5 )
		{
			r_val = FLA_LU_piv_opt_var5( A, p );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_LU_piv_blk_var3( A, p, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
		{
			r_val = FLA_LU_piv_blk_var4( A, p, cntl );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
		{
			r_val = FLA_LU_piv_blk_var5( A, p, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

