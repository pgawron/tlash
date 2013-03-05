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

dim_t FLA_CAQR_UT_inc_compute_blocks_per_part( dim_t p, FLA_Obj A )
{
	dim_t nb_part;
	dim_t nb_left;
	dim_t num_blocks;

	// Query the element (not scalar) length of A.
	num_blocks = FLA_Obj_length( A );
	
	// Compute the number of blocks per partitions.
	nb_part = num_blocks / p;
	nb_left = num_blocks % p;

	// If there are leftover blocks, increase nb_part by one.
	if ( nb_left > 0 ) nb_part += 1;

	return nb_part;
}

