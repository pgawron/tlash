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

FLA_Error FLA_Apply_Q_UT_rnfc_blk_var3( FLA_Obj A, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TWTL,  TWTR,     TW00, TW01, TW02, 
          TWBL,  TWBR,     TW10,  T11,  W12,
                           TW20, TW21, TW22;

  FLA_Obj WTL,   WTR,
          WBL,   WBR;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( TW,   &TWTL, &TWTR,
                      &TWBL, &TWBR,     0, 0, FLA_TL );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TWTL, /**/ TWTR,       &TW00, /**/ &TW01, &TW02,
                        /* *************** */   /* *********************** */
                                                  &TW10, /**/  &T11,  &W12,
                           TWBL, /**/ TWBR,       &TW20, /**/ &TW21, &TW22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( W,     &WTL, &WTR,
                         &WBL, &WBR,     b, FLA_Obj_length( B1 ), FLA_TL );

    // WTL = B1^T;

    FLA_Copyt_internal( FLA_TRANSPOSE, B1, WTL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // U11 = trilu( A11 );
    // U21 = A21;
    // Let WTL^T be conformal to B1.
    //
    // WTL^T = ( B1 * U11 + B2 * U21 ) * inv( triu(T11) );
    // WTL   = inv( triu(T11)^T ) * ( U11^T * B1^T + U21^T * B2^T );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm1( cntl ) );

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_TRANSPOSE, 
                       FLA_ONE, A21, B2, FLA_ONE, WTL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T11, WTL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // B2 = B2 - WTL^T * U21';
    // B1 = B1 - WTL^T * U11';
    //    = B1 - ( conj(U11) * WTL )^T;

    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_MINUS_ONE, WTL, A21, FLA_ONE, B2,
                       FLA_Cntl_sub_gemm2( cntl ) );

    FLA_Trmm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                       FLA_CONJ_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_MINUS_ONE, A11, WTL,
                       FLA_Cntl_sub_trmm2( cntl ) );

    FLA_Axpyt_internal( FLA_TRANSPOSE, FLA_ONE, WTL, B1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TWTL, /**/ &TWTR,       TW00, TW01, /**/ TW02,
                                                       TW10,  T11, /**/  W12,
                            /* **************** */  /* ********************* */
                              &TWBL, /**/ &TWBR,       TW20, TW21, /**/ TW22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

