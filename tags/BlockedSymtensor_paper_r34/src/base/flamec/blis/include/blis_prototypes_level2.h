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

// --- Level-2 BLAS-like prototypes --------------------------------------------

// --- gemv ---

void bli_sgemv( trans_t transa, conj_t conjx, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_dgemv( trans_t transa, conj_t conjx, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_cgemv( trans_t transa, conj_t conjx, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zgemv( trans_t transa, conj_t conjx, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bli_sgemv_blas( trans_t transa, int m, int n, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_dgemv_blas( trans_t transa, int m, int n, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_cgemv_blas( trans_t transa, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zgemv_blas( trans_t transa, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- ger ---

void bli_sger( conj_t conjx, conj_t conjy, int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bli_dger( conj_t conjx, conj_t conjy, int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bli_cger( conj_t conjx, conj_t conjy, int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bli_zger( conj_t conjx, conj_t conjy, int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bli_sger_blas(  int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bli_dger_blas(  int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bli_cgerc_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bli_cgeru_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bli_zgerc_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );
void bli_zgeru_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- hemv ---

void bli_shemv( uplo_t uplo, conj_t conj, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_dhemv( uplo_t uplo, conj_t conj, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_chemv( uplo_t uplo, conj_t conj, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zhemv( uplo_t uplo, conj_t conj, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bli_chemv_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zhemv_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- her ---

void bli_sher( uplo_t uplo, conj_t conj, int m, float*  alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bli_dher( uplo_t uplo, conj_t conj, int m, double* alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bli_cher( uplo_t uplo, conj_t conj, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bli_zher( uplo_t uplo, conj_t conj, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bli_cher_blas( uplo_t uplo, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int lda );
void bli_zher_blas( uplo_t uplo, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- her2 ---

void bli_sher2( uplo_t uplo, conj_t conj, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bli_dher2( uplo_t uplo, conj_t conj, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bli_cher2( uplo_t uplo, conj_t conj, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bli_zher2( uplo_t uplo, conj_t conj, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bli_cher2_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bli_zher2_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- symv ---

void bli_ssymv( uplo_t uplo, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_dsymv( uplo_t uplo, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_csymv( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zsymv( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bli_ssymv_blas( uplo_t uplo, int m, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bli_dsymv_blas( uplo_t uplo, int m, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bli_csymv_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_zsymv_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- syr ---

void bli_ssyr( uplo_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bli_dsyr( uplo_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bli_csyr( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bli_zsyr( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bli_ssyr_blas( uplo_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int lda );
void bli_dsyr_blas( uplo_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int lda );
void bli_csyr_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int lda );
void bli_zsyr_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- syr2 ---

void bli_ssyr2( uplo_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bli_dsyr2( uplo_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bli_csyr2( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bli_zsyr2( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bli_ssyr2_blas( uplo_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bli_dsyr2_blas( uplo_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bli_csyr2_blas( uplo_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bli_zsyr2_blas( uplo_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- trmv ---

void bli_strmv( uplo_t uplo, trans_t trans, diag_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bli_dtrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bli_ctrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bli_ztrmv( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bli_strmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, float*    a, int lda, float*    x, int incx );
void bli_dtrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, double*   a, int lda, double*   x, int incx );
void bli_ctrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bli_ztrmv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trsv ---

void bli_strsv( uplo_t uplo, trans_t trans, diag_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bli_dtrsv( uplo_t uplo, trans_t trans, diag_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bli_ctrsv( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bli_ztrsv( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bli_strsv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, float*    a, int lda, float*    x, int incx );
void bli_dtrsv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, double*   a, int lda, double*   x, int incx );
void bli_ctrsv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bli_ztrsv_blas( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trmvsx ---

void bli_strmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bli_dtrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bli_ctrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_ztrmvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- trsvsx ---

void bli_strsvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bli_dtrsvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bli_ctrsvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bli_ztrsvsx( uplo_t uplo, trans_t trans, diag_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

