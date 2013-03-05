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

void bli_strmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	float*    x_temp;
    int       incx_temp;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_screate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bli_sallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bli_scopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bli_strmv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bli_sscalv( BLIS_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bli_saxpyv( BLIS_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bli_sfree( x_temp );

	// Free the temporary contiguous matrix.
	bli_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_dtrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	double*   x_temp;
    int       incx_temp;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_dcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bli_dallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bli_dcopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bli_dtrmv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bli_dscalv( BLIS_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bli_daxpyv( BLIS_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bli_dfree( x_temp );

	// Free the temporary contiguous matrix.
	bli_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_ctrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_temp;
    int       incx_temp;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bli_callocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bli_ccopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bli_ctrmv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bli_cscalv( BLIS_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bli_caxpyv( BLIS_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bli_cfree( x_temp );

	// Free the temporary contiguous matrix.
	bli_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bli_ztrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_temp;
    int       incx_temp;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bli_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bli_zallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bli_zcopyv( BLIS_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bli_ztrmv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bli_zscalv( BLIS_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bli_zaxpyv( BLIS_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bli_zfree( x_temp );

	// Free the temporary contiguous matrix.
	bli_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

