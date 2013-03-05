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

void bli_ssyr2( uplo_t uplo, int m, float* alpha, float* x, int incx, float* y, int incy, float* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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

	bli_ssyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bli_sfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_dsyr2( uplo_t uplo, int m, double* alpha, double* x, int incx, double* y, int incy, double* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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

	bli_dsyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bli_dfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_csyr2( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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

	bli_csyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bli_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_zsyr2( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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

	bli_zsyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bli_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_ssyr2_blas( uplo_t uplo, int m, float* alpha, float* x, int incx, float* y, int incy, float* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssyr2( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssyr2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bli_dsyr2_blas( uplo_t uplo, int m, double* alpha, double* x, int incx, double* y, int incy, double* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsyr2( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsyr2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bli_csyr2_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda )
{
	scomplex* x_copy;
	scomplex* y_copy;
	scomplex  beta;
	int       k   = 1;
	int       ldx = m;
	int       ldy = m;

#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( BLIS_NO_TRANSPOSE, &cblas_trans );

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

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_csyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              x_copy, ldx,
	              y_copy, ldy,
	              &beta,
	              a,      lda );

	bli_cfree( x_copy );
	bli_cfree( y_copy );
#else
	char blas_uplo;
	char blas_trans;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( BLIS_NO_TRANSPOSE, &blas_trans );

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

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_csyr2k ( &blas_uplo,
	             &blas_trans,
	             &m,
	             &k,
	             alpha,
	             x_copy, &ldx,
	             y_copy, &ldy,
	             &beta,
	             a,      &lda );

	bli_cfree( x_copy );
	bli_cfree( y_copy );
#endif
}

void bli_zsyr2_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda )
{
	dcomplex* x_copy;
	dcomplex* y_copy;
	dcomplex  beta;
	int       k   = 1;
	int       ldx = m;
	int       ldy = m;

#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( BLIS_NO_TRANSPOSE, &cblas_trans );

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

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_zsyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              x_copy, ldx,
	              y_copy, ldy,
	              &beta,
	              a,      lda );

	bli_zfree( x_copy );
	bli_zfree( y_copy );
#else
	char blas_uplo;
	char blas_trans;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( BLIS_NO_TRANSPOSE, &blas_trans );

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

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_zsyr2k ( &blas_uplo,
	             &blas_trans,
	             &m,
	             &k,
	             alpha,
	             x_copy, &ldx,
	             y_copy, &ldy,
	             &beta,
	             a,      &lda );

	bli_zfree( x_copy );
	bli_zfree( y_copy );
#endif
}

