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

void bli_srandm( int m, int n, float* a, int a_rs, int a_cs )
{
	float*    a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem );
		bli_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_srandv( n_elem,
		            a_begin, inca );
	}
}

void bli_drandm( int m, int n, double* a, int a_rs, int a_cs )
{
	double*   a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem );
		bli_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_drandv( n_elem,
		            a_begin, inca );
	}
}

void bli_crandm( int m, int n, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem );
		bli_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_crandv( n_elem,
		            a_begin, inca );
	}
}

void bli_zrandm( int m, int n, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	int       inca, lda;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem );
		bli_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_zrandv( n_elem,
		            a_begin, inca );
	}
}

