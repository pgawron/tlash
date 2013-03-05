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

void bli_scopymrt( uplo_t uplo, trans_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
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
		
			bli_scopyv( conj,
			            n_elem,
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
		
			bli_scopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_dcopymrt( uplo_t uplo, trans_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
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
		
			bli_dcopyv( conj,
			            n_elem,
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
		
			bli_dcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_ccopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
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
		
			bli_ccopyv( conj,
			            n_elem,
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
		
			bli_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bli_zcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
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
		
			bli_zcopyv( conj,
			            n_elem,
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
		
			bli_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// --- Mixed-datatype and general stride copy routines---------------

// ss
void bli_sscopymrt( uplo_t uplo, trans_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
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
		
			bli_scopyv( conj,
			            n_elem,
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
		
			bli_scopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// sd
void bli_sdcopymrt( uplo_t uplo, trans_t trans, int m, int n, float* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	float*    a_begin;
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
		
			bli_sdcopyv( conj,
			             n_elem,
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
		
			bli_sdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// sc
void bli_sccopymrt( uplo_t uplo, trans_t trans, int m, int n, float* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
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
		
			bli_sccopyv( conj,
			             n_elem,
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
		
			bli_sccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// sz
void bli_szcopymrt( uplo_t uplo, trans_t trans, int m, int n, float* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
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
		
			bli_szcopyv( conj,
			             n_elem,
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
		
			bli_szcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// ds
void bli_dscopymrt( uplo_t uplo, trans_t trans, int m, int n, double* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	double*   a_begin;
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
		
			bli_dscopyv( conj,
			             n_elem,
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
		
			bli_dscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// dd
void bli_ddcopymrt( uplo_t uplo, trans_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
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
		
			bli_dcopyv( conj,
			            n_elem,
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
		
			bli_dcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// dc
void bli_dccopymrt( uplo_t uplo, trans_t trans, int m, int n, double* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
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
		
			bli_dccopyv( conj,
			             n_elem,
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
		
			bli_dccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// dz
void bli_dzcopymrt( uplo_t uplo, trans_t trans, int m, int n, double* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
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
		
			bli_dzcopyv( conj,
			             n_elem,
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
		
			bli_dzcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cs
void bli_cscopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
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
		
			bli_cscopyv( conj,
			             n_elem,
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
		
			bli_cscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cd
void bli_cdcopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
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
		
			bli_cdcopyv( conj,
			             n_elem,
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
		
			bli_cdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cc
void bli_cccopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
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
		
			bli_ccopyv( conj,
			            n_elem,
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
		
			bli_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// cz
void bli_czcopymrt( uplo_t uplo, trans_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
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
		
			bli_czcopyv( conj,
			             n_elem,
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
		
			bli_czcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zs
void bli_zscopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
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
		
			bli_zscopyv( conj,
			             n_elem,
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
		
			bli_zscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zd
void bli_zdcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
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
		
			bli_zdcopyv( conj,
			             n_elem,
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
		
			bli_zdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zc
void bli_zccopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
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
		
			bli_zccopyv( conj,
			             n_elem,
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
		
			bli_zccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zz
void bli_zzcopymrt( uplo_t uplo, trans_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
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
		
			bli_zcopyv( conj,
			            n_elem,
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
		
			bli_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

