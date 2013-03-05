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

void bli_srandmr( uplo_t uplo, diag_t diag, int m, int n, float* a, int a_rs, int a_cs )
{
	float*    a_begin;
	float*    ajj;
	float     one;
	float     zero;
	float     ord;
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

	// Initialize some scalars.
	one      = bli_s1();
	zero     = bli_s0();
	ord      = ( float ) bli_max( m, n );

	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bli_srandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bli_sinvscalv( BLIS_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_srands( ajj );
					bli_sabsval2( ajj, ajj );
					bli_sadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bli_ssetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bli_ssetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_srands( ajj );
					bli_sabsval2( ajj, ajj );
					bli_sadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bli_srandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bli_sinvscalv( BLIS_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bli_drandmr( uplo_t uplo, diag_t diag, int m, int n, double* a, int a_rs, int a_cs )
{
	double*   a_begin;
	double*   ajj;
	double    one;
	double    zero;
	double    ord;
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

	// Initialize some scalars.
	one      = bli_d1();
	zero     = bli_d0();
	ord      = ( double ) bli_max( m, n );

	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bli_drandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bli_dinvscalv( BLIS_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_drands( ajj );
					bli_dabsval2( ajj, ajj );
					bli_dadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bli_dsetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bli_dsetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_drands( ajj );
					bli_dabsval2( ajj, ajj );
					bli_dadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bli_drandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bli_dinvscalv( BLIS_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bli_crandmr( uplo_t uplo, diag_t diag, int m, int n, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_begin;
	scomplex* ajj;
	scomplex  one;
	scomplex  zero;
	scomplex  ord;
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

	// Initialize some scalars.
	one      = bli_c1();
	zero     = bli_c0();
	ord      = bli_c0();
	ord.real = ( float ) bli_max( m, n );

	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bli_crandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bli_cinvscalv( BLIS_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_crands( ajj );
					bli_cabsval2( ajj, ajj );
					bli_cadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bli_csetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bli_csetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_crands( ajj );
					bli_cabsval2( ajj, ajj );
					bli_cadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bli_crandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bli_cinvscalv( BLIS_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

void bli_zrandmr( uplo_t uplo, diag_t diag, int m, int n, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_begin;
	dcomplex* ajj;
	dcomplex  one;
	dcomplex  zero;
	dcomplex  ord;
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

	// Initialize some scalars.
	one      = bli_z1();
	zero     = bli_z0();
	ord      = bli_z0();
	ord.real = ( double ) bli_max( m, n );

	if ( bli_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Randomize super-diagonal elements.
			bli_zrandv( n_elem,
			            a_begin, inca );

			// Normalize super-diagonal elements by order of the matrix.
			bli_zinvscalv( BLIS_NO_CONJUGATE,
			               n_elem,
			               &ord,
			               a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_zrands( ajj );
					bli_zabsval2( ajj, ajj );
					bli_zadd3( ajj, &one, ajj );
				}

				// Initialize sub-diagonal elements to zero.
				bli_zsetv( n_elem_max - j - 1,
				           &zero,
				           ajj + inca, inca );
			}
		}
	}
	else // if ( bli_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bli_min( j, n_elem_max );
			a_begin = a + j*lda;

			// Initialize super-diagonal to zero.
			bli_zsetv( n_elem,
			           &zero,
			           a_begin, inca );

			// Initialize diagonal and sub-diagonal elements only if there are
			// elements left in the column (ie: j < n_elem_max).
			if ( j < n_elem_max )
			{
				ajj = a_begin + j*inca;

				// Initialize diagonal element.
				if      ( bli_is_unit_diag( diag ) )    *ajj = one;
				else if ( bli_is_zero_diag( diag ) )    *ajj = zero;
				else if ( bli_is_nonunit_diag( diag ) )
				{
					// We want positive diagonal elements between 1 and 2.
					bli_zrands( ajj );
					bli_zabsval2( ajj, ajj );
					bli_zadd3( ajj, &one, ajj );
				}

				// Randomize sub-diagonal elements.
				bli_zrandv( n_elem_max - j - 1,
				            ajj + inca, inca );

				// Normalize sub-diagonal elements by order of the matrix.
				bli_zinvscalv( BLIS_NO_CONJUGATE,
				               n_elem_max - j - 1,
				               &ord,
				               ajj + inca, inca );

			}
		}
	}
}

