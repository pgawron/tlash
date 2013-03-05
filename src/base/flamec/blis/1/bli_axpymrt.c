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

void bli_saxpymrt( uplo_t uplo, trans_t trans, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bli_is_col_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bli_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bli_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bli_does_trans( trans ) )
	{
		bli_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bli_saxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bli_saxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_daxpymrt( uplo_t uplo, trans_t trans, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bli_is_col_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bli_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bli_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bli_does_trans( trans ) )
	{
		bli_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bli_daxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bli_daxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_caxpymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bli_is_col_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bli_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bli_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bli_does_trans( trans ) )
	{
		bli_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bli_caxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bli_caxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_zaxpymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
	conj_t    conj;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bli_is_col_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bli_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		if ( bli_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bli_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bli_is_upper( uplo ) )
		{
			n_iter     = bli_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bli_does_trans( trans ) )
	{
		bli_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bli_proj_trans_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bli_zaxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bli_zaxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

