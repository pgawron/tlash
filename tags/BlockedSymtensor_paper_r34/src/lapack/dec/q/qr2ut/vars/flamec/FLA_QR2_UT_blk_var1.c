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

FLA_Error FLA_QR2_UT_blk_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl )
{
  FLA_Obj UTL,   UTR,      U00, U01, U02, 
          UBL,   UBR,      U10, U11, U12,
                           U20, U21, U22;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj TL,    TR,       T0,  T1,  W12;

  FLA_Obj W12T, W12B;

  FLA_Obj T1T, T2B;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( UBR ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( UBR ) );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // T1T = FLA_Top_part( T1, b );

    FLA_Part_2x1( T1,   &T1T, 
                        &T2B, b, FLA_TOP );

    // [ U11, ...
    //   D1, T1 ] = FLA_QR2_UT( U11
    //                          D1, T1T );

    FLA_QR2_UT_internal( U11,
                         D1, T1T, 
                         FLA_Cntl_sub_qr2ut( cntl ) );


    if ( FLA_Obj_width( U12 ) > 0 )
    {
      // W12T = FLA_Top_part( W12, b );

      FLA_Part_2x1( W12,  &W12T, 
                          &W12B, b, FLA_TOP );

      // W12T = inv( triu( T1T ) )' * ( U12 + D1' * D2 );

      FLA_Copy_internal( U12, W12T,
                         FLA_Cntl_sub_copy( cntl ) );

      FLA_Gemm_internal( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE,
                         FLA_ONE, D1, D2, FLA_ONE, W12T,
                         FLA_Cntl_sub_gemm1( cntl ) );

      FLA_Trsm_internal( FLA_LEFT, FLA_UPPER_TRIANGULAR, 
                         FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, 
                         FLA_ONE, T1T, W12T,
                         FLA_Cntl_sub_trsm( cntl ) );

      // U12 = U12 - W12T;
      // D2  = D2  - D1 * W12T;

      FLA_Axpy_internal( FLA_MINUS_ONE, W12T, U12,
                         FLA_Cntl_sub_axpy( cntl ) );

      FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
                         FLA_MINUS_ONE, D1, W12T, FLA_ONE, D2,
                         FLA_Cntl_sub_gemm2( cntl ) );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

