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

FLA_Error FLASH_Obj_create_diag_panel( FLA_Obj A, FLA_Obj* U )
{
	FLA_Datatype datatype;
	dim_t        b_flash;
	dim_t        b_flash_last;
	dim_t        n_blocks_min_dim;
	dim_t        m_U, n_U;

	// Acquire the datatype of the matrix to be factored.
	datatype = FLA_Obj_datatype( A );

	// Acquire the storage blocksize of the top-left element.
	b_flash = FLASH_Obj_scalar_length_tl( A );

	// Get the number of storage blocks in the minimum dimension of A.
	n_blocks_min_dim = FLA_Obj_min_dim( A );

	// Compute the scalar length and width of U.
	m_U = b_flash;
	n_U = n_blocks_min_dim * b_flash;

	// Create U with storage blocksize of b_flash.
	FLASH_Obj_create( datatype, m_U, n_U, 1, &b_flash, U );

	// The last, bottom-right-most diagonal block of A might be smaller
	// than the other diagonal blocks. Compute the size of this block.
	b_flash_last = FLASH_Obj_scalar_min_dim( A ) % b_flash;

	// If the remainder is zero, then A does not need its last block
	// shrunk and thus it is ready as-is. However, if b_flash_last is
	// non-zero, then we must manually adjust the size of the last block of
	// U. Note that we are not freeing and re-allocating memory, just
	// changing the size of the view into the last block.

	if ( b_flash_last > 0 )
	{
		FLA_Obj  UL,    UR;
		FLA_Obj  URTL,  URTR,
	    	     URBL,  URBR;
		FLA_Obj* UR_p;

		// Repartition U so we can access the last block object.
		FLA_Part_1x2( *U,   &UL, &UR,    1, FLA_RIGHT );

		// Dereference the 1x1 object reference to get a pointer to
		// the actual block object in U.
		UR_p = FLASH_OBJ_PTR_AT( UR );

		// Repartition the last block object so that URTL is the
		// correct size.
		FLA_Part_2x2( *UR_p,    &URTL, &URTR,
		                        &URBL, &URBR,     b_flash_last,
		                                          b_flash_last, FLA_TL );

		// Overwrite the original object pointed to by UR_p with the
		// corrected object URTL.
		*UR_p = URTL;
	}

	return FLA_SUCCESS;
}
