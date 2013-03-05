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

FLA_Error FLASH_UDdate_UT_inc_update_rhs( FLA_Obj T, FLA_Obj bR,
                                          FLA_Obj C, FLA_Obj bC,
                                          FLA_Obj D, FLA_Obj bD )
{
	FLA_Obj W;
	FLA_Obj bC_copy;
	FLA_Obj bD_copy;

	// Check parameters.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_UDdate_UT_inc_update_rhs_check( T, bR, C, bC, D, bD );

	// Create hierarchical workspace.
	FLASH_Apply_QUD_UT_inc_create_workspace( T, bR, &W );
	
	// Make temporary copies of the bC and bD right-hand side objects so we
	// don't destory their original contents.
	FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bC, &bC_copy );
	FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bD, &bD_copy );

	// Apply the updowndating Q' incrementally to the right-hand sides.
	FLASH_Apply_QUD_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
	                        T, W,
	                           bR,
	                        C, bC_copy,
	                        D, bD_copy );

	// Free the temporary objects.
	FLASH_Obj_free( &bC_copy );
	FLASH_Obj_free( &bD_copy );

	// Free the workspace object.
	FLASH_Obj_free( &W );

	return FLA_SUCCESS;
}

