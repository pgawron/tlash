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

void bli_sgemv( trans_t transa, conj_t conjx, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) )
	{
		int n_elem;

		if ( bli_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bli_sscalv( BLIS_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

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
		bli_toggle_trans( transa );
	}

	bli_sgemv_blas( transa,
	                m,
	                n,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_dgemv( trans_t transa, conj_t conjx, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) )
	{
		int n_elem;

		if ( bli_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bli_dscalv( BLIS_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

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
		bli_toggle_trans( transa );
	}

	bli_dgemv_blas( transa,
	                m,
	                n,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_cgemv( trans_t transa, conj_t conjx, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex  zero = bli_c0();
	scomplex  one  = bli_c1();
	scomplex* x_conj;
	scomplex* ax;
	int       lda, inca;
	int       n_x;
	int       incx_conj;
	int       incax;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) )
	{
		int n_elem;

		if ( bli_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bli_cscalv( BLIS_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

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
		bli_toggle_trans( transa );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated, and
	// also for the cases where A is conjugated.
	if ( bli_is_conj( conjx ) || bli_is_conjnotrans( transa ) )
	{
		if ( bli_does_trans( transa ) ) n_x = m;
		else                            n_x = n;

		x_conj    = bli_callocv( n_x );
		incx_conj = 1;

		bli_ccopyv( conjx,
		            n_x,
		            x,      incx,
                    x_conj, incx_conj );
	}

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	if ( bli_is_conjnotrans( transa ) )
	{
		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y. We know we are not transposing, so y is length m.
		ax    = bli_callocv( m );
		incax = 1;
		
		// Start by conjugating the contents of the temporary copy of x.
		bli_cconjv( n,
		            x_conj, incx_conj );

		// Compute A * conj(x) where x is the temporary copy of x created above.
		bli_cgemv_blas( BLIS_NO_TRANSPOSE,
		                m,
		                n,
		                &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax, incax );

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

		// Free the temporary vector for Ax.
		bli_cfree( ax );
	}
	else // notrans, trans, or conjtrans
	{
		bli_cgemv_blas( transa,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                x_conj, incx_conj,
		                beta,
		                y, incy );
	}

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conjx ) || bli_is_conjnotrans( transa ) )
		bli_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_zgemv( trans_t transa, conj_t conjx, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex  zero = bli_z0();
	dcomplex  one  = bli_z1();
	dcomplex* x_conj;
	dcomplex* ax;
	int       lda, inca;
	int       n_x;
	int       incx_conj;
	int       incax;

	// Return early if possible.
	if ( bli_zero_dim2( m, n ) )
	{
		int n_elem;

		if ( bli_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bli_zscalv( BLIS_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

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
		bli_toggle_trans( transa );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated, and
	// also for the cases where A is conjugated.
	if ( bli_is_conj( conjx ) || bli_is_conjnotrans( transa ) )
	{
		if ( bli_does_trans( transa ) ) n_x = m;
		else                            n_x = n;

		x_conj    = bli_zallocv( n_x );
		incx_conj = 1;

		bli_zcopyv( conjx,
		            n_x,
		            x,      incx,
                    x_conj, incx_conj );
	}

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	if ( bli_is_conjnotrans( transa ) )
	{
		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y. We know we are not transposing, so y is length m.
		ax    = bli_zallocv( m );
		incax = 1;
		
		// Start by conjugating the contents of the temporary copy of x.
		bli_zconjv( n,
		            x_conj, incx_conj );

		// Compute A * conj(x) where x is the temporary copy of x created above.
		bli_zgemv_blas( BLIS_NO_TRANSPOSE,
		                m,
		                n,
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

		// Free the temporary vector for Ax.
		bli_zfree( ax );
	}
	else // notrans, trans, or conjtrans
	{
		bli_zgemv_blas( transa,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                x_conj, incx_conj,
		                beta,
		                y,      incy );
	}

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conjx ) || bli_is_conjnotrans( transa ) )
		bli_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_sgemv_blas( trans_t transa, int m, int n, float* alpha, float* a, int lda, float* x, int incx, float* beta, float* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_sgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_transa;

	bli_param_map_to_netlib_trans( transa, &blas_transa );

	F77_sgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_dgemv_blas( trans_t transa, int m, int n, double* alpha, double* a, int lda, double* x, int incx, double* beta, double* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_dgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_transa;

	bli_param_map_to_netlib_trans( transa, &blas_transa );

	F77_dgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_cgemv_blas( trans_t transa, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_cgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_transa;

	bli_param_map_to_netlib_trans( transa, &blas_transa );

	F77_cgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_zgemv_blas( trans_t transa, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bli_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_zgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_transa;

	bli_param_map_to_netlib_trans( transa, &blas_transa );

	F77_zgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

