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

FLA_Error FLASH_Apply_Q_UT_create_workspace( FLA_Obj TW, FLA_Obj B, FLA_Obj* W )
{
	FLA_Datatype datatype;
	dim_t        depth;
	dim_t        b_alg;
	dim_t        b_flash;
	dim_t        m, n;

	// Query the depth.
	depth = FLASH_Obj_depth( TW );

	// *** The current Apply_Q_UT algorithm implemented assumes that
	// the matrix has a hierarchical depth of 1. We check for that here
	// because we anticipate that we'll use a more general algorithm in the
	// future, and we don't want to forget to remove the constraint. ***
	if ( depth != 1 )
	{
	   FLA_Print_message( "FLASH_Apply_Q_UT() currently only supports matrices of depth 1",
	                      __FILE__, __LINE__ );
	   FLA_Abort();
	}

	// Query the datatype of matrix TW.
	datatype = FLA_Obj_datatype( TW );

	// Inspect the dimensions of a the top-left element of TW to get the
	// algorithmic/storage blocksize we'll use throughout the Apply_Q_UT
	// algorithm.
	b_alg   = FLASH_Obj_scalar_length_tl( TW );
	b_flash = FLASH_Obj_scalar_width_tl( TW );

	// The traditional (non-incremental) Apply_Q_UT algorithm-by-blocks
	// requires that the algorithmic blocksize be equal to the storage
	// blocksize.
	if ( b_alg != b_flash )
	{
		FLA_Print_message( "FLASH_Apply_Q_UT() requires that b_alg == b_store",
		                   __FILE__, __LINE__ );
		FLA_Abort();
	}

	// The scalar length of W should be the algorithmc/storage blocksize
	// encoded in TW.
	m = b_alg;

	// Query the scalar (not element) width of the right-hand side
	// matrix B.
	n = FLASH_Obj_scalar_width( B );

	// Create hierarchical matrix W.
	FLASH_Obj_create_ext( datatype, m, n, 
	                      depth, &b_alg, &b_flash, 
	                      W );
   
	return FLA_SUCCESS;
}

