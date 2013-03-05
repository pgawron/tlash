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

FLA_Error FLA_LQ_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TTL,   TTR,      T00, T01, T02, 
          TBL,   TBR,      T10, T11, W12,
                           T20, T21, T22;

  FLA_Obj AR1,
          AR2;

  dim_t   b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00, /**/ &T01, &T02,
                        /* ************* */   /* ******************** */
                                                &T10, /**/ &T11, &W12,
                           TBL, /**/ TBR,       &T20, /**/ &T21, &T22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    FLA_Merge_1x2( A11, A12,   &AR1 );

    // Perform an LQ factorization via the UT transform on AR1:
    //
    //   ( A11 A12 ) -> L11 QR1
    //
    // where:
    //  - QR1 is formed from UR1 (which is stored row-wise above the
    //    diagonal of AR1) and T11 (which is stored to the upper triangle
    //    of T11).
    //  - L11 is stored to the lower triangle of AR1.
  
    FLA_LQ_UT_internal( AR1, T11, 
                        FLA_Cntl_sub_lqut( cntl ) );


    if ( FLA_Obj_length( A21 ) > 0 )
    {
      FLA_Merge_1x2( A21, A22,   &AR2 );

      // Apply the Householder transforms associated with UR1 and T11 to 
      // AR2:
      //
      //   ( A21 A22 ) := ( A21 A22 ) Q1
      //
      // where QR1 is formed from UR1 and T11.

      FLA_Apply_Q_UT_internal( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                               AR1, T11, W12, AR2,
                               FLA_Cntl_sub_apqut( cntl ) );
    }


    // Update T
    //
    // T01 = conj(A01) * U11^T + conj(A02) * U12^T;
    //
    // Recall: U11 = triuu( A11 );
    //         U12 = A12;

    FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, A01, T01 );
    FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                       FLA_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, T01 );
    FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_TRANSPOSE, 
                       FLA_ONE, A02, A12, FLA_ONE, T01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00, T01, /**/ T02,
                                                     T10, T11, /**/ W12,
                            /* ************** */  /* ****************** */
                              &TBL, /**/ &TBR,       T20, T21, /**/ T22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

