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

void bli_strmv( uplo_t uplo, trans_t trans, diag_t diag, int m, float* a, int a_rs, int a_cs, float* x, int incx )
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
		bli_toggle_trans( trans );
	}

	bli_strmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a, lda,
	                x, incx );

	// Free the temporary contiguous matrix.
	bli_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_dtrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, double* a, int a_rs, int a_cs, double* x, int incx )
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
		bli_toggle_trans( trans );
	}

	bli_dtrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a, lda,
	                x, incx );

	// Free the temporary contiguous matrix.
	bli_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_ctrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx )
{
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
		bli_toggle_trans( trans );
	}

	// Initialize with values assuming that trans is not conjnotrans.
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	// Note: strictly speaking, we don't need to create a copy of x since
	// the operation is simpler than, say, gemv. However, we create a copy
	// anyway since in practice it performs better due to increased spatial
	// locality.
	if ( bli_is_conjnotrans( trans ) )
	{
		x_conj    = bli_callocv( m );
		incx_conj = 1;

		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bli_ctrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a,      lda,
	                x_conj, incx_conj );

	// Save the contents of and then free the temporary conjugated x vector.
	if ( bli_is_conjnotrans( trans ) )
	{
		bli_ccopyv( BLIS_CONJUGATE,
                    m,
                    x_conj, incx_conj,
                    x,      incx );

		bli_cfree( x_conj );
	}

	// Free the temporary contiguous matrix.
	bli_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_ztrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx )
{
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
		bli_toggle_trans( trans );
	}

	// Initialize with values assuming that trans is not conjnotrans.
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	// Note: strictly speaking, we don't need to create a copy of x since
	// the operation is simpler than, say, gemv. However, we create a copy
	// anyway since in practice it performs better due to increased spatial
	// locality.
	if ( bli_is_conjnotrans( trans ) )
	{
		x_conj    = bli_zallocv( m );
		incx_conj = 1;

		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bli_ztrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a,      lda,
	                x_conj, incx_conj );

	// Save the contents of and then free the temporary conjugated x vector.
	if ( bli_is_conjnotrans( trans ) )
	{
		bli_zcopyv( BLIS_CONJUGATE,
                    m,
                    x_conj, incx_conj,
                    x,      incx );

		bli_zfree( x_conj );
	}

	// Free the temporary contiguous matrix.
	bli_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bli_strmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, float* a, int lda, float* x, int incx )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( trans, &cblas_trans );
	bli_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_strmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( trans, &blas_trans );
	bli_param_map_to_netlib_diag( diag, &blas_diag );

	F77_strmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bli_dtrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, double* a, int lda, double* x, int incx )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( trans, &cblas_trans );
	bli_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_dtrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( trans, &blas_trans );
	bli_param_map_to_netlib_diag( diag, &blas_diag );

	F77_dtrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bli_ctrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int lda, scomplex* x, int incx )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( trans, &cblas_trans );
	bli_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ctrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( trans, &blas_trans );
	bli_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ctrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bli_ztrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bli_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bli_param_map_to_netlib_trans( trans, &cblas_trans );
	bli_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ztrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bli_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bli_param_map_to_netlib_trans( trans, &blas_trans );
	bli_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ztrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

