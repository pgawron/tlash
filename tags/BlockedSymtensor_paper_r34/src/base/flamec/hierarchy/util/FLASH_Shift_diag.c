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

FLA_Error FLASH_Shift_diag( FLA_Conj conj, FLA_Obj sigma, FLA_Obj H )
{
	FLA_Obj F;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( H ) ) return FLA_SUCCESS;

	// Create a temporary flat copy of the hierarchical object.
	FLASH_Obj_create_flat_copy_of_hier( H, &F );

	// Shift the diagonal of the flat matrix object by sigma.
	FLA_Shift_diag( conj, sigma, F );

	// Copy the flat object's contents back to the hierarchical object.
	FLASH_Obj_hierarchify( F, H );
	
	// Free the temporary flat object.
	FLA_Obj_free( &F );
	
	return FLA_SUCCESS;
}
