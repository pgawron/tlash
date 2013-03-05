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

FLA_Error FLASH_LQ_UT_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, FLA_Obj* A, FLA_Obj* TW )
{
	FLA_Datatype datatype;
	dim_t        m, n;
	dim_t        min_m_n;
	
	// *** The current LQ_UT algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1. We check for that here
	// because we anticipate that we'll use a more general algorithm in the
	// future, and we don't want to forget to remove the constraint. ***
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_LQ_UT() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Create hierarchical copy of matrix A_flat.
	FLASH_Obj_create_hier_copy_of_flat( A_flat, depth, b_flash, A );

	// Query the datatype of matrix A_flat.
	datatype = FLA_Obj_datatype( A_flat );
	
	// Query the minimum dimension of A_flat.
	min_m_n = FLA_Obj_min_dim( A_flat );

	// Set the m and n dimensions of TW to be min_m_n.
	m = min_m_n;
	n = min_m_n;

	// Create hierarchical matrices T and W.
	FLASH_Obj_create_ext( datatype, m, n, 
	                      depth, b_flash, b_flash, 
	                      TW );
	   
	return FLA_SUCCESS;
}

