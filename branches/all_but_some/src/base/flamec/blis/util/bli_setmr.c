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

void bli_ssetmr( uplo_t uplo, int m, int n, float* sigma, float* a, int a_rs, int a_cs )
{
	float*    a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
	}
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			bli_ssetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bli_ssetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bli_dsetmr( uplo_t uplo, int m, int n, double* sigma, double* a, int a_rs, int a_cs )
{
	double*   a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
	}
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			bli_dsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bli_dsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bli_csetmr( uplo_t uplo, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
	}
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			bli_csetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bli_csetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bli_zsetmr( uplo_t uplo, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
	}
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			bli_zsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bli_zsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

