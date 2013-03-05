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

FLA_Error FLA_Apply_Q2_UT_lnfc_blk_var1( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
                                                                           FLA_Obj E, fla_apq2ut_t* cntl )
{
  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj T1T,
          T2B;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj W1TL,  W1TR,
          W1BL,  W1BR;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_RIGHT );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_RIGHT );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width( DR ) < FLA_Obj_width( D ) ){

    b = min( b_alg, FLA_Obj_width( DL ) );

    // Since T was filled from left to right, and since we need to access them
    // in reverse order, we need to handle the case where the last block is
    // smaller than the other b x b blocks.
    if ( FLA_Obj_width( TR ) == 0 && FLA_Obj_width( T ) % b_alg > 0 )
      b = FLA_Obj_width( T ) % b_alg;

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, &D1, /**/ &D2,
                           b, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, &T1, /**/ &T2,
                           b, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &C1, 
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,    &T1T, 
                         &T2B,     b, FLA_TOP );

    FLA_Part_2x2( W1,    &W1TL, &W1TR,
                         &W1BL, &W1BR,     b, FLA_Obj_width( C1 ), FLA_TL );

    // W1TL = C1;

    FLA_Copyt_internal( FLA_NO_TRANSPOSE, C1, W1TL,
                        FLA_Cntl_sub_copyt( cntl ) );

    // W1TL = inv( triu( T1T ) ) * ( C1 + D1' * E );

    FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       FLA_ONE, D1, E, FLA_ONE, W1TL,
                       FLA_Cntl_sub_gemm1( cntl ) );

    FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, T1T, W1TL,
                       FLA_Cntl_sub_trsm( cntl ) );

    // C1 = C1 - W1TL;
    // E  = E  - D1 * W1TL;

    FLA_Axpyt_internal( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, W1TL, C1,
                        FLA_Cntl_sub_axpyt( cntl ) );

    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, D1, W1TL, FLA_ONE, E,
                       FLA_Cntl_sub_gemm2( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, /**/ D1, D2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, /**/ T1, T2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* ** */
                                                  C1, 
                              &CB,                C2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

