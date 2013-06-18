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

void bli_sger( conj_t conjx, conj_t conjy, int m, int n, float* alpha, float* x, int incx, float* y, int incy, float* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_screate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( m, n );
		bli_swap_ints( lda, inca );
		bli_swap_ints( incx, incy );
		bli_swap_conj( conjx, conjy );
		bli_sswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	bli_sger_blas( m,
	               n,
	               alpha,
	               x, incx,
	               y, incy,
	               a, lda );

	// Free the temporary contiguous matrix.
	bli_sfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_dger( conj_t conjx, conj_t conjy, int m, int n, double* alpha, double* x, int incx, double* y, int incy, double* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_dcreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( m, n );
		bli_swap_ints( lda, inca );
		bli_swap_ints( incx, incy );
		bli_swap_conj( conjx, conjy );
		bli_dswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	bli_dger_blas( m,
	               n,
	               alpha,
	               x, incx,
	               y, incy,
	               a, lda );

	// Free the temporary contiguous matrix.
	bli_dfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_cger( conj_t conjx, conj_t conjy, int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_conj;
	int       incx_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_ccreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( m, n );
		bli_swap_ints( lda, inca );
		bli_swap_ints( incx, incy );
		bli_swap_conj( conjx, conjy );
		bli_cswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated.
	if ( bli_is_conj( conjx ) )
	{
		x_conj    = bli_callocv( m );
		incx_conj = 1;

		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	// Conjugation of y is supported in the BLAS.
	if ( bli_is_conj( conjy ) )
	{
		bli_cgerc_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}
	else
	{
		bli_cgeru_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conjx ) )
		bli_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_cfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_zger( conj_t conjx, conj_t conjy, int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_conj;
	int       incx_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_zcreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( m, n );
		bli_swap_ints( lda, inca );
		bli_swap_ints( incx, incy );
		bli_swap_conj( conjx, conjy );
		bli_zswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated.
	if ( bli_is_conj( conjx ) )
	{
		x_conj    = bli_zallocv( m );
		incx_conj = 1;

		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	// Conjugation of y is supported in the BLAS.
	if ( bli_is_conj( conjy ) )
	{
		bli_zgerc_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}
	else
	{
		bli_zgeru_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conjx ) )
		bli_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_zfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_sger_blas( int m, int n, float* alpha, float* x, int incx, float* y, int incy, float* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_sger( cblas_order,
	            m,
	            n,
	            *alpha,
	            x, incx,
	            y, incy,
	            a, lda );
#else
	F77_sger( &m,
	          &n,
	          alpha,
	          x, &incx,
	          y, &incy,
	          a, &lda );
#endif
}

void bli_dger_blas( int m, int n, double* alpha, double* x, int incx, double* y, int incy, double* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_dger( cblas_order,
	            m,
	            n,
	            *alpha,
	            x, incx,
	            y, incy,
	            a, lda );
#else
	F77_dger( &m,
	          &n,
	          alpha,
	          x, &incx,
	          y, &incy,
	          a, &lda );
#endif
}

void bli_cgerc_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_cgerc( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_cgerc ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bli_cgeru_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_cgeru( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_cgeru ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bli_zgerc_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_zgerc( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_zgerc ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bli_zgeru_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_zgeru( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_zgeru ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

