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

void bli_sher2( uplo_t uplo, conj_t conj, int m, float* alpha, float* x, int incx, float* y, int incy, float* a, int a_rs, int a_cs )
{
	bli_ssyr2( uplo,
	           m,
	           alpha,
	           x, incx,
	           y, incy,
	           a, a_rs, a_cs );
}

void bli_dher2( uplo_t uplo, conj_t conj, int m, double* alpha, double* x, int incx, double* y, int incy, double* a, int a_rs, int a_cs )
{
	bli_dsyr2( uplo,
	           m,
	           alpha,
	           x, incx,
	           y, incy,
	           a, a_rs, a_cs );
}

void bli_cher2( uplo_t uplo, conj_t conj, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_conj;
	scomplex* y_conj;
	int       incx_conj;
	int       incy_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
		bli_toggle_conj( conj );
	}

	// Initialize with values assuming no conjugation of ( x * y' ) or
	// ( y * x' ).
	x_conj    = x;
	incx_conj = incx;
	y_conj    = y;
	incy_conj = incy;

	// We want to handle the case where ( x * y' ) and ( y * x' ) are
	// conjugated, but without explicitly conjugating the matrices. To do
	// so, we leverage the fact that computing the products conj( x * y' )
	// and conj( y * x' ) is equivalent to computing ( conj(x) * conj(y)' )
	// and ( conj(y) * conj(x)' ), respectively.
	if ( bli_is_conj( conj ) )
	{
		x_conj    = bli_callocv( m );
		incx_conj = 1;

		y_conj    = bli_callocv( m );
		incy_conj = 1;

		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );

		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    y,      incy,
                    y_conj, incy_conj );
	}

	bli_cher2_blas( uplo,
	                m,
	                alpha,
	                x_conj, incx_conj,
	                y_conj, incy_conj,
	                a,      lda );

	// Free the temporary conjugated x and y vectors.
	if ( bli_is_conj( conj ) )
	{
		bli_cfree( x_conj );
		bli_cfree( y_conj );
	}

	// Free the temporary contiguous matrix.
	bli_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_zher2( uplo_t uplo, conj_t conj, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_conj;
	dcomplex* y_conj;
	int       incx_conj;
	int       incy_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bli_is_row_storage( a_rs, a_cs ) )
	{
		bli_swap_ints( lda, inca );
		bli_toggle_uplo( uplo );
		bli_toggle_conj( conj );
	}

	// Initialize with values assuming no conjugation of ( x * y' ) or
	// ( y * x' ).
	x_conj    = x;
	incx_conj = incx;
	y_conj    = y;
	incy_conj = incy;

	// We want to handle the case where ( x * y' ) and ( y * x' ) are
	// conjugated, but without explicitly conjugating the matrices. To do
	// so, we leverage the fact that computing the products conj( x * y' )
	// and conj( y * x' ) is equivalent to computing ( conj(x) * conj(y)' )
	// and ( conj(y) * conj(x)' ), respectively.
	if ( bli_is_conj( conj ) )
	{
		x_conj    = bli_zallocv( m );
		incx_conj = 1;

		y_conj    = bli_zallocv( m );
		incy_conj = 1;

		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );

		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    y,      incy,
                    y_conj, incy_conj );
	}

	bli_zher2_blas( uplo,
	                m,
	                alpha,
	                x_conj, incx_conj,
	                y_conj, incy_conj,
	                a,      lda );

	// Free the temporary conjugated x and y vectors.
	if ( bli_is_conj( conj ) )
	{
		bli_zfree( x_conj );
		bli_zfree( y_conj );
	}

	// Free the temporary contiguous matrix.
	bli_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_cher2_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_cher2( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_cher2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bli_zher2_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zher2( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zher2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

