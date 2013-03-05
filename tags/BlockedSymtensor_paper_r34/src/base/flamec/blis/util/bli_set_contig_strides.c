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

#include "blis.h"

void bli_set_contig_strides( int m, int n, int* rs, int* cs )
{
	// Default to column-major order.
	*rs = 1;
	*cs = m;

	// Handle special cases first.
	// Check the strides, and modify them if needed.
	if ( *rs == 1 && *cs == 1 )
	{
		// If both strides are unit, we are probably trying to create a
		// 1-by-n matrix in column-major order, or an m-by-1 matrix in
		// row-major order. We have decided to "reserve" the case where
		// rs == cs == 1 for scalars only, as having unit strides can
		// upset the BLAS error checking when attempting to induce a
		// row-major operation.
		if ( m > 1 && n == 1 )
		{
			// Set the column stride to indicate that this is an m-by-1
			// matrix (or vector) stored in column-major order. This is
			// necessary because, in some cases, we have to satisfy error
			// checking in the underlying BLAS library, which expects the
			// leading dimension to be set to at least m, even if it will
			// never be used for indexing since there is only one column
			// of data. Note that rs is already set to 1.
			*cs = m;
		}
		else if ( m == 1 && 1 < n )
		{
			// Set the row stride to indicate that this is a 1-by-n matrix
			// stored in row-major order. Note that cs is already set to 1.
			*rs = n;
		}
		else
		{
			// If m == n == 1, then we are dealing with a scalar. Since rs
			// and cs do not exceed m and n, we don't have to do anything.
		}
	}
}

