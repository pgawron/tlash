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

// --- Level-3 BLAS-like prototypes --------------------------------------------

// --- gemm ---

void bli_sgemm( trans_t transa, trans_t transb, int m, int k, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dgemm( trans_t transa, trans_t transb, int m, int k, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_cgemm( trans_t transa, trans_t transb, int m, int k, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_zgemm( trans_t transa, trans_t transb, int m, int k, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bli_sgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bli_dgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bli_cgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bli_zgemm_blas( trans_t transa, trans_t transb, int m, int n, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- hemm ---

void bli_shemm( side_t side, uplo_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dhemm( side_t side, uplo_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_chemm( side_t side, uplo_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_zhemm( side_t side, uplo_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bli_chemm_blas( side_t side, uplo_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bli_zhemm_blas( side_t side, uplo_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- herk ---

void bli_sherk( uplo_t uplo, trans_t trans, int m, int k, float*  alpha, float*    a, int a_rs, int a_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bli_dherk( uplo_t uplo, trans_t trans, int m, int k, double* alpha, double*   a, int a_rs, int a_cs, double* beta, double*   c, int c_rs, int c_cs );
void bli_cherk( uplo_t uplo, trans_t trans, int m, int k, float*  alpha, scomplex* a, int a_rs, int a_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bli_zherk( uplo_t uplo, trans_t trans, int m, int k, double* alpha, dcomplex* a, int a_rs, int a_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bli_cherk_blas( uplo_t uplo, trans_t trans, int m, int k, float*  alpha, scomplex* a, int lda, float*  beta, scomplex* c, int ldc );
void bli_zherk_blas( uplo_t uplo, trans_t trans, int m, int k, double* alpha, dcomplex* a, int lda, double* beta, dcomplex* c, int ldc );

// --- her2k ---

void bli_sher2k( uplo_t uplo, trans_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bli_dher2k( uplo_t uplo, trans_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double* beta, double*   c, int c_rs, int c_cs );
void bli_cher2k( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bli_zher2k( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bli_cher2k_blas( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, float*  beta, scomplex* c, int ldc );
void bli_zher2k_blas( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, double* beta, dcomplex* c, int ldc );

// --- symm ---

void bli_ssymm( side_t side, uplo_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dsymm( side_t side, uplo_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_csymm( side_t side, uplo_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_zsymm( side_t side, uplo_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bli_ssymm_blas( side_t side, uplo_t uplo, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bli_dsymm_blas( side_t side, uplo_t uplo, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bli_csymm_blas( side_t side, uplo_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bli_zsymm_blas( side_t side, uplo_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- syrk ---

void bli_ssyrk( uplo_t uplo, trans_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dsyrk( uplo_t uplo, trans_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_csyrk( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_zsyrk( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bli_ssyrk_blas( uplo_t uplo, trans_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    beta, float*    c, int ldc );
void bli_dsyrk_blas( uplo_t uplo, trans_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   beta, double*   c, int ldc );
void bli_csyrk_blas( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* beta, scomplex* c, int ldc );
void bli_zsyrk_blas( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* beta, dcomplex* c, int ldc );

// --- syr2k ---

void bli_ssyr2k( uplo_t uplo, trans_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dsyr2k( uplo_t uplo, trans_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_csyr2k( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_zsyr2k( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bli_ssyr2k_blas( uplo_t uplo, trans_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bli_dsyr2k_blas( uplo_t uplo, trans_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bli_csyr2k_blas( uplo_t uplo, trans_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bli_zsyr2k_blas( uplo_t uplo, trans_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- trmm ---

void bli_strmm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dtrmm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_ctrmm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_ztrmm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bli_strmm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bli_dtrmm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bli_ctrmm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bli_ztrmm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trsm ---

void bli_strsm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bli_dtrsm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bli_ctrsm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bli_ztrsm( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bli_strsm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bli_dtrsm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bli_ctrsm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bli_ztrsm_blas( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trmmsx ---

void bli_strmmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dtrmmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_ctrmmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_ztrmmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

// --- trsmsx ---

void bli_strsmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bli_dtrsmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bli_ctrsmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bli_ztrsmsx( side_t side, uplo_t uplo, trans_t trans, diag_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

