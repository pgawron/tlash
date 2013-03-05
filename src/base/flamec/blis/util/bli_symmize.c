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

void bli_ssymmize( conj_t conj, uplo_t uplo, int m, float* a, int a_rs, int a_cs )
{
	float*    a_src;
	float*    a_dst;
	int       rs_src, cs_src, inc_src;
	int       rs_dst, cs_dst, inc_dst;
	int       n_iter;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bli_is_col_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bli_is_col_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bli_scopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );
	}
}

void bli_dsymmize( conj_t conj, uplo_t uplo, int m, double* a, int a_rs, int a_cs )
{
	double*   a_src;
	double*   a_dst;
	int       rs_src, cs_src, inc_src;
	int       rs_dst, cs_dst, inc_dst;
	int       n_iter;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bli_is_col_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bli_is_col_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bli_dcopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );
	}
}

void bli_csymmize( conj_t conj, uplo_t uplo, int m, scomplex* a, int a_rs, int a_cs )
{
	scomplex* a_src;
	scomplex* a_dst;
	scomplex* a_jj;
	int       rs_src, cs_src, inc_src;
	int       rs_dst, cs_dst, inc_dst;
	int       n_iter;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bli_is_col_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bli_is_col_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bli_ccopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );

		if ( bli_is_conj( conj ) )
		{
			a_jj = a + j*a_rs + j*a_cs;
			a_jj->imag = bli_s0();
		}
	}
}

void bli_zsymmize( conj_t conj, uplo_t uplo, int m, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* a_src;
	dcomplex* a_dst;
	dcomplex* a_jj;
	int       rs_src, cs_src, inc_src;
	int       rs_dst, cs_dst, inc_dst;
	int       n_iter;
	int       j;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Assume A is square.
	n_iter = m;

	// Initialize with appropriate values based on storage.
	if      ( bli_is_col_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 1;
		rs_src  = 0;
		inc_src = a_cs;
		cs_dst  = a_cs;
		rs_dst  = 0;
		inc_dst = 1;
	}
	else if ( bli_is_col_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = a_cs;
		rs_src  = 0;
		inc_src = 1;
		cs_dst  = 1;
		rs_dst  = 0;
		inc_dst = a_cs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		cs_src  = 0;
		rs_src  = a_rs;
		inc_src = 1;
		cs_dst  = 0;
		rs_dst  = 1;
		inc_dst = a_rs;
	}
	else if ( bli_is_row_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		cs_src  = 0;
		rs_src  = 1;
		inc_src = a_rs;
		cs_dst  = 0;
		rs_dst  = a_rs;
		inc_dst = 1;
	}
	else if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_lower( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = 1 * a_rs;
			rs_src  = 0;
			inc_src = a_cs;
			cs_dst  = a_cs;
			rs_dst  = 0;
			inc_dst = 1 * a_rs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = a_rs;
			inc_src = 1 * a_cs;
			cs_dst  = 0;
			rs_dst  = 1 * a_cs;
			inc_dst = a_rs;
		}
	}
	else // if ( bli_is_gen_storage( a_rs, a_cs ) && bli_is_upper( uplo ) )
	{
		// General stride with column-major tilt looks similar to column-major.
		// General stride with row-major tilt looks similar to row-major.
		if ( a_rs < a_cs )
		{
			cs_src  = a_cs;
			rs_src  = 0;
			inc_src = 1 * a_rs;
			cs_dst  = 1 * a_rs;
			rs_dst  = 0;
			inc_dst = a_cs;
		}
		else // if ( a_rs > a_cs )
		{
			cs_src  = 0;
			rs_src  = 1 * a_cs;
			inc_src = a_rs;
			cs_dst  = 0;
			rs_dst  = a_rs;
			inc_dst = 1 * a_cs;
		}
	}
	
	for ( j = 0; j < n_iter; j++ )
	{
		a_src = a + j*cs_src + j*rs_src;
		a_dst = a + j*cs_dst + j*rs_dst;

		bli_zcopyv( conj,
		            j,
		            a_src, inc_src,
		            a_dst, inc_dst );

		if ( bli_is_conj( conj ) )
		{
			a_jj = a + j*a_rs + j*a_cs;
			a_jj->imag = bli_d0();
		}
	}
}

