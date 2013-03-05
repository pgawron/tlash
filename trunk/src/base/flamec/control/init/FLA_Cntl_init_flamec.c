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

#include "FLAME.h"

void FLA_Cntl_init_flamec()
{
  FLA_Transpose_cntl_init();

  // Level-1 BLAS
  FLA_Axpy_cntl_init();
  FLA_Axpyt_cntl_init();
  FLA_Copy_cntl_init();
  FLA_Copyt_cntl_init();
  FLA_Copyr_cntl_init();
  FLA_Scal_cntl_init();
  FLA_Scalr_cntl_init();

  // Level-2 BLAS
  FLA_Gemv_cntl_init();
  FLA_Trsv_cntl_init();

  // Level-3 BLAS
  // Note gemm must be first since it is used by all other level-3 BLAS.
  FLA_Gemm_cntl_init();
  FLA_Hemm_cntl_init();
  FLA_Herk_cntl_init();
  FLA_Her2k_cntl_init();
  FLA_Symm_cntl_init();
  FLA_Syrk_cntl_init();
  FLA_Syr2k_cntl_init();
  FLA_Trmm_cntl_init();
  FLA_Trsm_cntl_init();

//Original libFLAME lines
/*
  // LAPACK-level
  // These require level-3 BLAS operations to be initialized.
  FLA_Apply_pivots_cntl_init();
  FLA_Chol_cntl_init();
  FLA_LU_nopiv_cntl_init();
  FLA_LU_piv_cntl_init();
  FLA_Trinv_cntl_init();
  FLA_Ttmm_cntl_init();
  FLA_Sylv_cntl_init();
  FLA_QR2_UT_cntl_init();
  FLA_CAQR2_UT_cntl_init();
  FLA_Apply_Q_UT_cntl_init();
  FLA_Apply_Q2_UT_cntl_init();
  FLA_Apply_CAQ2_UT_cntl_init();
  FLA_Apply_QUD_UT_cntl_init();
  FLA_Eig_gest_cntl_init();

  // Compound LAPACK operations
  // These require previous LAPACK operations to already be initialized.
  FLA_Lyap_cntl_init();
  FLA_SPDinv_cntl_init();
  FLA_QR_UT_cntl_init();
  FLA_LQ_UT_cntl_init();
  FLA_UDdate_UT_cntl_init();
  FLA_Hess_UT_cntl_init();
  FLA_Tridiag_UT_cntl_init();
  FLA_Bidiag_UT_cntl_init();
*/
}

void FLA_Cntl_finalize_flamec()
{
  FLA_Transpose_cntl_finalize();

  // Level-1 BLAS
  FLA_Axpy_cntl_finalize();
  FLA_Axpyt_cntl_finalize();
  FLA_Copy_cntl_finalize();
  FLA_Copyt_cntl_finalize();
  FLA_Copyr_cntl_finalize();
  FLA_Scal_cntl_finalize();
  FLA_Scalr_cntl_finalize();

  // Level-2 BLAS
  FLA_Gemv_cntl_finalize();
  FLA_Trsv_cntl_finalize();

  // Level-3 BLAS
  FLA_Gemm_cntl_finalize();
  FLA_Hemm_cntl_finalize();
  FLA_Herk_cntl_finalize();
  FLA_Her2k_cntl_finalize();
  FLA_Symm_cntl_finalize();
  FLA_Syrk_cntl_finalize();
  FLA_Syr2k_cntl_finalize();
  FLA_Trmm_cntl_finalize();
  FLA_Trsm_cntl_finalize();

//Original libFLAME lines
/*
  // LAPACK-level
  FLA_Apply_pivots_cntl_finalize();
  FLA_Chol_cntl_finalize();
  FLA_LU_nopiv_cntl_finalize();
  FLA_LU_piv_cntl_finalize();
  FLA_Trinv_cntl_finalize();
  FLA_Ttmm_cntl_finalize();
  FLA_Sylv_cntl_finalize();
  FLA_QR2_UT_cntl_finalize();
  FLA_CAQR2_UT_cntl_finalize();
  FLA_Apply_Q_UT_cntl_finalize();
  FLA_Apply_Q2_UT_cntl_finalize();
  FLA_Apply_CAQ2_UT_cntl_finalize();
  FLA_Apply_QUD_UT_cntl_finalize();
  FLA_Eig_gest_cntl_finalize();

  // Compound LAPACK operations
  FLA_Lyap_cntl_finalize();
  FLA_SPDinv_cntl_finalize();
  FLA_QR_UT_cntl_finalize();
  FLA_LQ_UT_cntl_finalize();
  FLA_UDdate_UT_cntl_finalize();
  FLA_Hess_UT_cntl_finalize();
  FLA_Tridiag_UT_cntl_finalize();
  FLA_Bidiag_UT_cntl_finalize();
*/
}

