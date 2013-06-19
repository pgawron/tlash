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

void bli_sher( uplo_t uplo, conj_t conj, int m, float* alpha, float* x, int incx, float* a, int a_rs, int a_cs )
{
	bli_ssyr( uplo,
	          m,
	          alpha,
	          x, incx,
	          a, a_rs, a_cs );
}

void bli_dher( uplo_t uplo, conj_t conj, int m, double* alpha, double* x, int incx, double* a, int a_rs, int a_cs )
{
	bli_dsyr( uplo,
	          m,
	          alpha,
	          x, incx,
	          a, a_rs, a_cs );
}

void bli_cher( uplo_t uplo, conj_t conj, int m, float* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_conj;
	int       incx_conj;
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

	// Initialize with values assuming no conjugation of ( x * x' ).
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the case where ( x * x' ) is conjugated, but
	// without explicitly conjugating the matrix. To do so, we leverage
	// the fact that computing the product conj( x * x' ) is equivalent
	// to computing ( conj(x) * conj(x)' ), since ( x * x' ) is Hermitian.
	if ( bli_is_conj( conj ) )
	{
		x_conj    = bli_callocv( m );
		incx_conj = 1;

		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bli_cher_blas( uplo,
	               m,
	               alpha,
	               x_conj, incx_conj,
	               a,      lda );

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conj ) )
		bli_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bli_zher( uplo_t uplo, conj_t conj, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_conj;
	int       incx_conj;
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

	// Initialize with values assuming no conjugation of ( x * x' ).
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the case where ( x * x' ) is conjugated, but
	// without explicitly conjugating the matrix. To do so, we leverage
	// the fact that computing the product conj( x * x' ) is equivalent
	// to computing ( conj(x) * conj(x)' ), since ( x * x' ) is Hermitian.
	if ( bli_is_conj( conj ) )
	{
		x_conj    = bli_zallocv( m );
		incx_conj = 1;

		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bli_zher_blas( uplo,
	               m,
	               alpha,
	               x_conj, incx_conj,
	               a,      lda );

	// Free the temporary conjugated x vector.
	if ( bli_is_conj( conj ) )
		bli_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bli_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_cher_blas( uplo_t uplo, int m, float* alpha, scomplex* x, int incx, scomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_cher( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_cher( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

void bli_zher_blas( uplo_t uplo, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int lda )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zher( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zher( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

