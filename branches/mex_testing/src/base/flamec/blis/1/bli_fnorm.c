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

void bli_sfnorm( int m, int n, float* a, int a_rs, int a_cs, float* norm )
{
	float*    a_ij;
	float     sum;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       i, j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
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
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0F;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += (*a_ij) * (*a_ij);
		}
	}
	
	// Compute the norm and store the result.
	*norm = ( float ) sqrt( sum );
}

void bli_dfnorm( int m, int n, double* a, int a_rs, int a_cs, double* norm )
{
	double*   a_ij;
	double    sum;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       i, j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
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
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += (*a_ij) * (*a_ij);
		}
	}
	
	// Compute the norm and store the result.
	*norm = sqrt( sum );
}

void bli_cfnorm( int m, int n, scomplex* a, int a_rs, int a_cs, float* norm )
{
	scomplex* a_ij;
	float     sum;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       i, j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
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
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0F;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += a_ij->real * a_ij->real + a_ij->imag * a_ij->imag;
		}
	}
	
	// Compute the norm and store the result.
	*norm = ( float ) sqrt( sum );
}

void bli_zfnorm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* norm )
{
	dcomplex* a_ij;
	double    sum;
	int       lda, inca;
	int       n_iter;
	int       n_elem;
	int       i, j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector separately.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
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
		
		// An optimization: if A is row-major, then let's access the matrix by
		// rows instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( a_rs, a_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
		}
	}

	// Initialize the accumulator variable.
	sum = 0.0;

	for ( j = 0; j < n_iter; j++ )
	{
		for ( i = 0; i < n_elem; i++ )
		{
			a_ij = a + i*inca + j*lda;
			sum += a_ij->real * a_ij->real + a_ij->imag * a_ij->imag;
		}
	}
	
	// Compute the norm and store the result.
	*norm = sqrt( sum );
}

