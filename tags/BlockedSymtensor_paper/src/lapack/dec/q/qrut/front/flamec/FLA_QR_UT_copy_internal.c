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

FLA_Error FLA_QR_UT_copy_internal( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_QR_UT_copy_internal_check( A, T, U, cntl );

	if ( FLASH_Queue_get_enabled() )
	{
		// Enqueue task.
		ENQUEUE_FLASH_QR_UT_copy( *FLASH_OBJ_PTR_AT( A ),
		                          *FLASH_OBJ_PTR_AT( T ),
		                          *FLASH_OBJ_PTR_AT( U ),
		                          NULL );
	}
	else
	{
		// Execute task immediately.
		FLA_QR_UT_copy_task( *FLASH_OBJ_PTR_AT( A ),
		                     *FLASH_OBJ_PTR_AT( T ),
		                     *FLASH_OBJ_PTR_AT( U ),
                             NULL );
	}

	return r_val;
}

