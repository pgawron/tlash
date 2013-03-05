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

void bli_ssymv( uplo_t uplo, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_screate_contigmr( uplo,
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
	}

	bli_ssymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_dsymv( uplo_t uplo, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_dcreate_contigmr( uplo,
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
	}

	bli_dsymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_csymv( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
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
	}

	bli_csymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_zsymv( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
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
	}

	bli_zsymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bli_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_ssymv_blas( uplo_t uplo, int m, float* alpha, float* a, int lda, float* x, int incx, float* beta, float* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssymv( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssymv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_dsymv_blas( uplo_t uplo, int m, double* alpha, double* a, int lda, double* x, int incx, double* beta, double* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsymv( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsymv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bli_csymv_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* x_copy;
	scomplex* y_copy;
	int       n   = 1;
	int       ldx = m;
	int       ldy = m;

#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_side( BLIS_LEFT, &cblas_side );
	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	x_copy = bli_callocv( m );
	y_copy = bli_callocv( m );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	cblas_csymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a,      lda,
	             x_copy, ldx,
	             beta,
	             y_copy, ldy );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bli_cfree( x_copy );
	bli_cfree( y_copy );

#else
	char blas_side;
	char blas_uplo;

	bli_param_map_to_netlib_side( BLIS_LEFT, &blas_side );
	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	x_copy = bli_callocv( m );
	y_copy = bli_callocv( m );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	F77_csymm ( &blas_side,
	            &blas_uplo,
	            &m,
	            &n,
	            alpha,
	            a,      &lda,
	            x_copy, &ldx,
	            beta,
	            y_copy, &ldy );

	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bli_cfree( x_copy );
	bli_cfree( y_copy );
#endif
}

void bli_zsymv_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex*        x_copy;
	dcomplex*        y_copy;
	int              n   = 1;
	int              ldx = m;
	int              ldy = m;

#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_side( BLIS_LEFT, &cblas_side );
	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	x_copy = bli_zallocv( m );
	y_copy = bli_zallocv( m );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	cblas_zsymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a,      lda,
	             x_copy, ldx,
	             beta,
	             y_copy, ldy );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bli_zfree( x_copy );
	bli_zfree( y_copy );

#else
	char blas_side;
	char blas_uplo;

	bli_param_map_to_netlib_side( BLIS_LEFT, &blas_side );
	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	x_copy = bli_zallocv( m );
	y_copy = bli_zallocv( m );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	F77_zsymm ( &blas_side,
	            &blas_uplo,
	            &m,
	            &n,
	            alpha,
	            a,      &lda,
	            x_copy, &ldx,
	            beta,
	            y_copy, &ldy );

	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bli_zfree( x_copy );
	bli_zfree( y_copy );
#endif
}

