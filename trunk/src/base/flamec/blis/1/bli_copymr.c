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

void bli_scopymr( uplo_t uplo, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) && bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_scopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_scopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bli_dcopymr( uplo_t uplo, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) && bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_dcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_dcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bli_ccopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) && bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_ccopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_ccopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

void bli_zcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if A and B are both row-major, then let's access the
	// matrices by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) && bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_zcopy( n_elem,
			           a_begin, inca, 
			           b_begin, incb );
		}
	}
}

// --- Mixed-datatype and general stride copy routines---------------

// ss
void bli_sscopymr( uplo_t uplo, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_scopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_scopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// sd ds
void bli_sdcopymr( uplo_t uplo, int m, int n, float* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	float*    a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_sdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_sdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_dscopymr( uplo_t uplo, int m, int n, double* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	double*   a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_dscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_dscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// sc cs
void bli_sccopymr( uplo_t uplo, int m, int n, float* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_sccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_sccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_cscopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_cscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_cscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// sz zs
void bli_szcopymr( uplo_t uplo, int m, int n, float* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_szcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_szcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_zscopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_zscopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// dd
void bli_ddcopymr( uplo_t uplo, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_dcopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_dcopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// dc cd
void bli_dccopymr( uplo_t uplo, int m, int n, double* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_dccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_dccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_cdcopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_cdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_cdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// dz zd
void bli_dzcopymr( uplo_t uplo, int m, int n, double* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_dzcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_dzcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_zdcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_zdcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// cc
void bli_cccopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_ccopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_ccopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

// cz zc
void bli_czcopymr( uplo_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_czcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_czcopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}
void bli_zccopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_zccopyv( BLIS_NO_TRANSPOSE,
			             n_elem,
			             a_begin, inca, 
			             b_begin, incb );
		}
	}
}

// zz
void bli_zzcopymr( uplo_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem_max;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// We initialize for column-major.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;
	ldb        = b_cs;
	incb       = b_rs;

	// An optimization: if B is row-major, then let's access the matrix
	// by rows instead of by columns for increased spatial locality.
	if ( bli_is_row_storage( b_rs, b_cs ) )
	{
		bli_swap_ints( n_iter, n_elem_max );
		bli_swap_ints( lda, inca );
		bli_swap_ints( ldb, incb );
		bli_toggle_uplo( uplo );
	}
	
	
	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;

			bli_zcopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;

			if ( n_elem <= 0 ) break;

			bli_zcopyv( BLIS_NO_TRANSPOSE,
			            n_elem,
			            a_begin, inca, 
			            b_begin, incb );
		}
	}
}

