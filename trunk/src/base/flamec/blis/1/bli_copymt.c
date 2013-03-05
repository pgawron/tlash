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

void bli_icopymt( trans_t trans, int m, int n, int* a, int a_rs, int a_cs, int* b, int b_rs, int b_cs )
{
	int*      a_begin;
	int*      b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// after a possible transposition, then let's access the matrix by rows
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
		
		bli_icopyv( trans,
		            n_elem,
		            a_begin, inca, 
		            b_begin, incb );
	}
}

void bli_scopymt( trans_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
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

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// after a possible transposition, then let's access the matrix by rows
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
		
		bli_scopy( n_elem,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bli_dcopymt( trans_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
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

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// after a possible transposition, then let's access the matrix by rows
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
		
		bli_dcopy( n_elem,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bli_ccopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// after a possible transposition, then let's access the matrix by rows
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
		
		bli_ccopy( n_elem,
		           a_begin, inca, 
		           b_begin, incb );

		if ( bli_does_conj( trans ) )
			bli_cconjv( n_elem,
			            b_begin, incb );
	}
}

void bli_zcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// after a possible transposition, then let's access the matrix by rows
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
		
		bli_zcopy( n_elem,
		           a_begin, inca, 
		           b_begin, incb );

		if ( bli_does_conj( trans ) )
			bli_zconjv( n_elem,
			            b_begin, incb );
	}
}

// --- Mixed-datatype and general stride copy routines---------------

// ss
void bli_sscopymt( trans_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_scopyv( conj,
		            n_elem,
		            a_begin, inca,
		            b_begin, incb );
	}
}

// sd ds
void bli_sdcopymt( trans_t trans, int m, int n, float* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	float*    a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_sdcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_dscopymt( trans_t trans, int m, int n, double* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	double*   a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_dscopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// sc cs
void bli_sccopymt( trans_t trans, int m, int n, float* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_sccopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_cscopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_cscopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// sz zs
void bli_szcopymt( trans_t trans, int m, int n, float* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_szcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_zscopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_zscopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// dd
void bli_ddcopymt( trans_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_dcopyv( conj,
		            n_elem,
		            a_begin, inca,
		            b_begin, incb );
	}
}

// dc cd
void bli_dccopymt( trans_t trans, int m, int n, double* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_dccopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_cdcopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_cdcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// dz zd
void bli_dzcopymt( trans_t trans, int m, int n, double* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_dzcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_zdcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_zdcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// cc
void bli_cccopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_ccopyv( conj,
		            n_elem,
		            a_begin, inca,
		            b_begin, incb );
	}
}

// cz zc
void bli_czcopymt( trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_czcopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}
void bli_zccopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_zccopyv( conj,
		             n_elem,
		             a_begin, inca,
		             b_begin, incb );
	}
}

// zz
void bli_zzcopymt( trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
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
		// Initialize with optimal values for column-major storage of B.
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

		// An optimization: if B is row-major, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bli_is_row_storage( b_rs, b_cs ) )
		{
			bli_swap_ints( n_iter, n_elem );
			bli_swap_ints( lda, inca );
			bli_swap_ints( ldb, incb );
		}
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	for ( j = 0; j < n_iter; ++j )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bli_zcopyv( conj,
		            n_elem,
		            a_begin, inca,
		            b_begin, incb );
	}
}

