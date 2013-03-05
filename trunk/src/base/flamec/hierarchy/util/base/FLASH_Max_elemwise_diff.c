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

double FLASH_Max_elemwise_diff( FLA_Obj A, FLA_Obj B )
{
	FLA_Obj A_flat, B_flat;
	double  max_diff;

	// Exit early if one dimension is zero.
	if ( FLA_Obj_has_zero_dim( A ) ) return -1.0;

	// Create a temporary flat copy of the hierarchical objects.
	FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
	FLASH_Obj_create_flat_copy_of_hier( B, &B_flat );

	// Get the maximum element-wise diff.
	max_diff = FLA_Max_elemwise_diff( A_flat, B_flat );
	
	// Free the temporary flat objects.
	FLA_Obj_free( &A_flat );
	FLA_Obj_free( &B_flat );
	
	return max_diff;
}

