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

void FLA_Cntl_init_flash( void );
void FLA_Cntl_finalize_flash( void );


// --- Base library prototypes -------------------------------------------------
void FLASH_Transpose_cntl_init( void );

void FLASH_Transpose_cntl_finalize( void );


// --- Level-1 BLAS prototypes -------------------------------------------------
void FLASH_Axpy_cntl_init( void );
void FLASH_Axpyt_cntl_init( void );
void FLASH_Copy_cntl_init( void );
void FLASH_Copyt_cntl_init( void );
void FLASH_Copyr_cntl_init( void );
void FLASH_Scal_cntl_init( void );
void FLASH_Scalr_cntl_init( void );

void FLASH_Axpy_cntl_finalize( void );
void FLASH_Axpyt_cntl_finalize( void );
void FLASH_Copy_cntl_finalize( void );
void FLASH_Copyt_cntl_finalize( void );
void FLASH_Copyr_cntl_finalize( void );
void FLASH_Scal_cntl_finalize( void );
void FLASH_Scalr_cntl_finalize( void );


// --- Level-2 BLAS prototypes -------------------------------------------------
void FLASH_Gemv_cntl_init( void );
void FLASH_Trsv_cntl_init( void );

void FLASH_Gemv_cntl_finalize( void );
void FLASH_Trsv_cntl_finalize( void );


// --- Level-3 BLAS prototypes -------------------------------------------------
void FLASH_Gemm_cntl_init( void );
void FLASH_Hemm_cntl_init( void );
void FLASH_Herk_cntl_init( void );
void FLASH_Her2k_cntl_init( void );
void FLASH_Symm_cntl_init( void );
void FLASH_Syrk_cntl_init( void );
void FLASH_Syr2k_cntl_init( void );
void FLASH_Trmm_cntl_init( void );
void FLASH_Trsm_cntl_init( void );

void FLASH_Gemm_cntl_finalize( void );
void FLASH_Hemm_cntl_finalize( void );
void FLASH_Herk_cntl_finalize( void );
void FLASH_Her2k_cntl_finalize( void );
void FLASH_Symm_cntl_finalize( void );
void FLASH_Syrk_cntl_finalize( void );
void FLASH_Syr2k_cntl_finalize( void );
void FLASH_Trmm_cntl_finalize( void );
void FLASH_Trsm_cntl_finalize( void );


// --- LAPACK-level prototypes -------------------------------------------------
/*
void FLASH_Apply_pivots_cntl_init( void );
void FLASH_Chol_cntl_init( void );
void FLASH_LU_nopiv_cntl_init( void );
void FLASH_LU_piv_cntl_init( void );
void FLASH_LU_incpiv_cntl_init( void );
void FLASH_Trinv_cntl_init( void );
void FLASH_Ttmm_cntl_init( void );
void FLASH_SPDinv_cntl_init( void );
void FLASH_Sylv_cntl_init( void );
void FLASH_Lyap_cntl_init( void );
void FLASH_QR_UT_cntl_init( void );
void FLASH_QR2_UT_cntl_init( void );
void FLASH_LQ_UT_cntl_init( void );
void FLASH_CAQR2_UT_cntl_init( void );
void FLASH_UDdate_UT_cntl_init( void );
void FLASH_QR_UT_inc_cntl_init( void );
void FLASH_CAQR_UT_inc_cntl_init( void );
void FLASH_UDdate_UT_inc_cntl_init( void );
void FLASH_Apply_Q_UT_cntl_init( void );
void FLASH_Apply_Q2_UT_cntl_init( void );
void FLASH_Apply_CAQ2_UT_cntl_init( void );
void FLASH_Apply_QUD_UT_cntl_init( void );
void FLASH_Apply_Q_UT_inc_cntl_init( void );
void FLASH_Apply_CAQ_UT_inc_cntl_init( void );
void FLASH_Apply_QUD_UT_inc_cntl_init( void );
void FLASH_Eig_gest_cntl_init( void );

void FLASH_Apply_pivots_cntl_finalize( void );
void FLASH_Chol_cntl_finalize( void );
void FLASH_LU_nopiv_cntl_finalize( void );
void FLASH_LU_piv_cntl_finalize( void );
void FLASH_LU_incpiv_cntl_finalize( void );
void FLASH_Trinv_cntl_finalize( void );
void FLASH_Ttmm_cntl_finalize( void );
void FLASH_SPDinv_cntl_finalize( void );
void FLASH_Sylv_cntl_finalize( void );
void FLASH_Lyap_cntl_finalize( void );
void FLASH_QR_UT_cntl_finalize( void );
void FLASH_QR2_UT_cntl_finalize( void );
void FLASH_LQ_UT_cntl_finalize( void );
void FLASH_CAQR2_UT_cntl_finalize( void );
void FLASH_UDdate_UT_cntl_finalize( void );
void FLASH_QR_UT_inc_cntl_finalize( void );
void FLASH_CAQR_UT_inc_cntl_finalize( void );
void FLASH_UDdate_UT_inc_cntl_finalize( void );
void FLASH_Apply_Q_UT_cntl_finalize( void );
void FLASH_Apply_Q2_UT_cntl_finalize( void );
void FLASH_Apply_CAQ2_UT_cntl_finalize( void );
void FLASH_Apply_QUD_UT_cntl_finalize( void );
void FLASH_Apply_Q_UT_inc_cntl_finalize( void );
void FLASH_Apply_CAQ_UT_inc_cntl_finalize( void );
void FLASH_Apply_QUD_UT_inc_cntl_finalize( void );
void FLASH_Eig_gest_cntl_finalize( void );
*/
