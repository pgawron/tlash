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

// --- top-level wrapper prototypes --------------------------------------------

// Implemented:
FLA_Error FLASH_Chol( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_nopiv( FLA_Obj A );
FLA_Error FLASH_LU_nopiv_solve( FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_piv( FLA_Obj A, FLA_Obj p );
FLA_Error FLASH_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L );
FLA_Error FLASH_FS_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b );
FLA_Error FLASH_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLASH_Ttmm( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_SPDinv( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLASH_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLASH_Apply_Q2_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E );
FLA_Error FLASH_QR2_UT( FLA_Obj B, FLA_Obj D, FLA_Obj T );
FLA_Error FLASH_QR_UT( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_QR_UT_inc( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_inc_solve( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_Apply_Q_UT_inc( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B );
FLA_Error FLASH_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );
FLA_Error FLASH_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

// Not yet implemented:
FLA_Error FLASH_LQ_UT_inv( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_LQ2_UT( FLA_Obj B, FLA_Obj C, FLA_Obj T );
