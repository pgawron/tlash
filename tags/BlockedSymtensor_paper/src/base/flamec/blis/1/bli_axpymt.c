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

void bli_saxpymt( trans_t trans, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bli_does_trans( trans ) )
		{
			bli_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bli_is_col_storage( a_rs, a_cs ) && bli_does_trans( trans ) ) ||
			     ( bli_is_row_storage( a_rs, a_cs ) && bli_does_notrans( trans ) ) )
			{
				bli_swap_ints( n_iter, n_elem );
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_saxpy( n_elem,
		           alpha,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bli_daxpymt( trans_t trans, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bli_does_trans( trans ) )
		{
			bli_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bli_is_col_storage( a_rs, a_cs ) && bli_does_trans( trans ) ) ||
			     ( bli_is_row_storage( a_rs, a_cs ) && bli_does_notrans( trans ) ) )
			{
				bli_swap_ints( n_iter, n_elem );
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_daxpy( n_elem,
		           alpha,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bli_caxpymt( trans_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	scomplex* a_temp;
	int       inca_temp;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bli_does_trans( trans ) )
		{
			bli_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bli_is_col_storage( a_rs, a_cs ) && bli_does_trans( trans ) ) ||
			     ( bli_is_row_storage( a_rs, a_cs ) && bli_does_notrans( trans ) ) )
			{
				bli_swap_ints( n_iter, n_elem );
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
			}
		}
	}

	if ( bli_does_conj( trans ) )
	{
		conj_t conj = bli_proj_trans_to_conj( trans );

		a_temp = bli_callocv( n_elem );
		inca_temp = 1;

		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            a_temp,  inca_temp );

			bli_caxpy( n_elem,
			           alpha,
			           a_temp,  inca_temp, 
			           b_begin, incb );
		}

		bli_cfree( a_temp );
	}
	else // if ( !bli_does_conj( trans ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_caxpy( n_elem,
			           alpha,
			           a_begin, inca, 
			           b_begin, incb );
		}
	
	}
}

void bli_zaxpymt( trans_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	dcomplex* a_temp;
	int       inca_temp;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying axpy
	// gets invoked only once.
	if ( bli_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bli_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bli_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bli_vector_inc( BLIS_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bli_does_trans( trans ) )
		{
			bli_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrices by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bli_is_col_storage( a_rs, a_cs ) && bli_does_trans( trans ) ) ||
			     ( bli_is_row_storage( a_rs, a_cs ) && bli_does_notrans( trans ) ) )
			{
				bli_swap_ints( n_iter, n_elem );
				bli_swap_ints( lda, inca );
				bli_swap_ints( ldb, incb );
			}
		}
	}

	if ( bli_does_conj( trans ) )
	{
		conj_t conj = bli_proj_trans_to_conj( trans );

		a_temp = bli_zallocv( n_elem );
		inca_temp = 1;

		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            a_temp,  inca_temp );

			bli_zaxpy( n_elem,
			           alpha,
			           a_temp,  inca_temp, 
			           b_begin, incb );
		}

		bli_zfree( a_temp );
	}
	else // if ( !bli_does_conj( trans ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zaxpy( n_elem,
			           alpha,
			           a_begin, inca, 
			           b_begin, incb );
		}
	
	}
}

