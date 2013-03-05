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

FLA_Error FLA_Eig_gest_iu_blk_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02,
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj BTL,   BTR,      B00, B01, B02,
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  FLA_Obj YT,              Y01,
          YB,              Y11,
                           Y21;

  FLA_Obj Y01_l, Y01_r;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x1( Y,    &YT,
                      &YB,            0, FLA_TOP );

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

    FLA_Repart_2x1_to_3x1( YT,                  &Y01,
                        /* ** */              /* *** */
                                                &Y11,
                           YB,                  &Y21,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( Y01,    &Y01_l, &Y01_r,    b, FLA_LEFT );

    // Y01 = A00 * B01;
    FLA_Hemm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_ONE, A00, B01, FLA_ZERO, Y01_l,
                       FLA_Cntl_sub_hemm( cntl ) );

    // A01 = inv( triu( B00 )' ) * A01;
    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B00, A01,
                       FLA_Cntl_sub_trsm1( cntl ) );

    // A01 = A01 - 1/2 * Y01;
    FLA_Axpy_internal( FLA_MINUS_ONE_HALF, Y01_l, A01,
                       FLA_Cntl_sub_axpy1( cntl ) );

    // A11 = A11 - A01' * B01 - B01' * A01;
    FLA_Her2k_internal( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                        FLA_MINUS_ONE, A01, B01, FLA_ONE, A11,
                        FLA_Cntl_sub_her2k( cntl ) );

    // A11 = inv( triu( B11 )' ) * A11 * inv( triu( B11 ) );
    FLA_Eig_gest_internal( FLA_INVERSE, FLA_UPPER_TRIANGULAR,
                           A11, Y11, B11,
                           FLA_Cntl_sub_eig_gest( cntl ) );

    // A01 = A01 - 1/2 * Y01;
    FLA_Axpy_internal( FLA_MINUS_ONE_HALF, Y01_l, A01,
                       FLA_Cntl_sub_axpy2( cntl ) );

    // A01 = A01 * inv( triu( B11 ) );
    FLA_Trsm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, B11, A01,
                       FLA_Cntl_sub_trsm2( cntl ) );
    
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

    FLA_Cont_with_3x1_to_2x1( &YT,                   Y01,
                                                     Y11,
                            /* ** */              /* *** */
                              &YB,                   Y21,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

#endif
