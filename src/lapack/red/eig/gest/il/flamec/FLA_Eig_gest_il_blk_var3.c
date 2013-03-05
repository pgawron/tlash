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

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_Eig_gest_il_blk_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02,
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj BTL,   BTR,      B00, B01, B02,
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  FLA_Obj YTL,   YTR,      Y00, Y01, Y02, 
          YBL,   YBR,      Y10, Y11, Y12,
                           Y20, Y21, Y22;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x2( Y,    &YTL, &YTR,
                      &YBL, &YBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00, /**/ &B01, &B02,
                        /* ************* */   /* ******************** */
                                                &B10, /**/ &B11, &B12,
                           BBL, /**/ BBR,       &B20, /**/ &B21, &B22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( YTL, /**/ YTR,       &Y00, /**/ &Y01, &Y02,
                        /* ************* */   /* ******************** */
                                                &Y10, /**/ &Y11, &Y12,
                           YBL, /**/ YBR,       &Y20, /**/ &Y21, &Y22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    // A10 = A10 - 1/2 * Y10;
    FLA_Axpy_internal( FLA_MINUS_ONE_HALF, Y10, A10,
                       FLA_Cntl_sub_axpy1( cntl ) );

    // A11 = A11 - A10 * B10' - B10 * A10';
    FLA_Her2k_internal( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
                        FLA_MINUS_ONE, A10, B10, FLA_ONE, A11,
                        FLA_Cntl_sub_her2k( cntl ) );

    // A11 = inv( tril( B11 ) ) * A11 * inv( tril( B11 )' );
    FLA_Eig_gest_internal( FLA_INVERSE, FLA_LOWER_TRIANGULAR,
                           A11, Y11, B11,
                           FLA_Cntl_sub_eig_gest( cntl ) );

    // A21 = A21 - A20 * B10';
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_MINUS_ONE, A20, B10, FLA_ONE, A21,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // A21 = A21 * inv( tril( B11 )' );
    FLA_Trsm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B11, A21,
                       FLA_Cntl_sub_trsm1( cntl ) );

    // A10 = A10 - 1/2 * Y10;
    FLA_Axpy_internal( FLA_MINUS_ONE_HALF, Y10, A10,
                       FLA_Cntl_sub_axpy2( cntl ) );

    // A10 = inv( tril( B11 ) ) * A10;
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B11, A10,
                       FLA_Cntl_sub_trsm2( cntl ) );

    // Y20 = Y20 + B21 * A10;
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_ONE, B21, A10, FLA_ONE, Y20,
                       FLA_Cntl_sub_gemm2( cntl ) );

    // Y21 = B21 * A11;
    FLA_Hemm_internal( FLA_RIGHT, FLA_LOWER_TRIANGULAR,
                       FLA_ONE, A11, B21, FLA_ZERO, Y21,
                       FLA_Cntl_sub_hemm( cntl ) );

    // Y21 = Y21 + B20 * A10';
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, B20, A10, FLA_ONE, Y21,
                       FLA_Cntl_sub_gemm3( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00, B01, /**/ B02,
                                                     B10, B11, /**/ B12,
                            /* ************** */  /* ****************** */
                              &BBL, /**/ &BBR,       B20, B21, /**/ B22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &YTL, /**/ &YTR,       Y00, Y01, /**/ Y02,
                                                     Y10, Y11, /**/ Y12,
                            /* ************** */  /* ****************** */
                              &YBL, /**/ &YBR,       Y20, Y21, /**/ Y22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

#endif