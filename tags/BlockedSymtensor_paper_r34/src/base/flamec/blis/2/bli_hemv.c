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

void bli_shemv( uplo_t uplo, conj_t conj, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	bli_ssymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bli_dhemv( uplo_t uplo, conj_t conj, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	bli_dsymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bli_chemv( uplo_t uplo, conj_t conj, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex  zero = bli_c0();
	scomplex  one  = bli_c1();
	scomplex* x_conj;
	scomplex* ax;
	int       lda, inca;
	int       incx_conj;
	int       incax;

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

	// We want to handle the case where A is conjugated, but without
	// explicitly or conjugating A. To do so, we leverage the fact that
	// computing the product conj(A) * x is equivalent to computing
	// conj( A * conj(x) ).
	if ( bli_is_conj( conj ) )
	{
		// We need a temporary vector so we can create a conjugated copy of x.
		x_conj    = bli_callocv( m );
		incx_conj = 1;

		bli_ccopyv( BLIS_CONJUGATE,
		            m,
		            x,      incx,
                    x_conj, incx_conj );

		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y (and x).
		ax    = bli_callocv( m );
		incax = 1;
		
		// Compute A * conj(x) where x is the temporary copy of x created above.
		bli_chemv_blas( uplo,
		                m,
                        &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax,     incax );

		// Scale y by beta.
		bli_cscalv( BLIS_NO_CONJUGATE,
                    m,
                    beta,
                    y, incy );

		// And finally, accumulate alpha * conj( A * conj(x) ) into y.
		bli_caxpyv( BLIS_CONJUGATE,
                    m,
		            alpha,
                    ax, incax,
                    y,  incy);

		// Free the temporary vectors for x and Ax.
		bli_cfree( x_conj );
		bli_cfree( ax );
	}
	else // noconj
	{
		bli_chemv_blas( uplo,
		                m,
		                alpha,
		                a, lda,
		                x, incx,
		                beta,
		                y, incy );
	}

	// Free the temporary contiguous matrix.
	bli_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_zhemv( uplo_t uplo, conj_t conj, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex  zero = bli_z0();
	dcomplex  one  = bli_z1();
	dcomplex* x_conj;
	dcomplex* ax;
	int       lda, inca;
	int       incx_conj;
	int       incax;

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

	// We want to handle the case where A is conjugated, but without
	// explicitly or conjugating A. To do so, we leverage the fact that
	// computing the product conj(A) * x is equivalent to computing
	// conj( A * conj(x) ).
	if ( bli_is_conj( conj ) )
	{
		// We need a temporary vector so we can create a conjugated copy of x.
		x_conj    = bli_zallocv( m );
		incx_conj = 1;

		bli_zcopyv( BLIS_CONJUGATE,
		            m,
		            x,      incx,
                    x_conj, incx_conj );

		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y (and x).
		ax    = bli_zallocv( m );
		incax = 1;
		
		// Compute A * conj(x) where x is the temporary copy of x created above.
		bli_zhemv_blas( uplo,
		                m,
		                &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax,     incax );

		// Scale y by beta.
		bli_zscalv( BLIS_NO_CONJUGATE,
                    m,
                    beta,
                    y, incy );

		// And finally, accumulate alpha * conj( A * conj(x) ) into y.
		bli_zaxpyv( BLIS_CONJUGATE,
                    m,
                    alpha,
                    ax, incax,
                    y,  incy);

		// Free the temporary vectors for x and Ax.
		bli_zfree( x_conj );
		bli_zfree( ax );
	}
	else // noconj
	{
		bli_zhemv_blas( uplo,
		                m,
		                alpha,
		                a, lda,
		                x, incx,
		                beta,
		                y, incy );
	}

	// Free the temporary contiguous matrix.
	bli_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_chemv_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_chemv( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_chemv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_zhemv_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zhemv( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zhemv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

