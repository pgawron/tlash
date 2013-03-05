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

void bli_sinvscalm( conj_t conj, int m, int n, float* alpha, float* a, int a_rs, int a_cs )
{
	float     alpha_inv;
	float*    a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_seq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_sinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_sscal( n_elem,
		           &alpha_inv,
		           a_begin, inca );
	}
}

void bli_dinvscalm( conj_t conj, int m, int n, double* alpha, double* a, int a_rs, int a_cs )
{
	double    alpha_inv;
	double*   a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_deq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_dinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_dscal( n_elem,
		           &alpha_inv,
		           a_begin, inca );
	}
}

void bli_csinvscalm( conj_t conj, int m, int n, float* alpha, scomplex* a, int a_rs, int a_cs )
{
	float     alpha_inv;
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_seq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_sinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_csscal( n_elem,
		            &alpha_inv,
		            a_begin, inca );
	}
}

void bli_cinvscalm( conj_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs )
{
	scomplex  alpha_inv;
	scomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_ceq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_cinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_cscal( n_elem,
		           &alpha_inv,
		           a_begin, inca );
	}
}

void bli_zdinvscalm( conj_t conj, int m, int n, double* alpha, dcomplex* a, int a_rs, int a_cs )
{
	double    alpha_inv;
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_deq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_dinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_zdscal( n_elem,
		            &alpha_inv,
		            a_begin, inca );
	}
}

void bli_zinvscalm( conj_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex  alpha_inv;
	dcomplex* a_begin;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;
	if ( bli_zeq1( alpha ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	bli_zinvert2s( conj, alpha, &alpha_inv );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bli_zscal( n_elem,
		           &alpha_inv,
		           a_begin, inca );
	}
}

