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

FLA_Error FLA_LU_nopiv_blk_var3( FLA_Obj A, fla_lu_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  dim_t b;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
          FLA_Obj_width( ATL ) < FLA_Obj_width( A )){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    // A01 = trilu( A00 ) \ A10
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A00, A01,
                       FLA_Cntl_sub_trsm1( cntl ) );

    // A11 = LU( A11 - A10 * A01 )
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A01, FLA_ONE, A11,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // A11 = LU_nopiv( A11 )
    FLA_LU_nopiv_internal( A11, FLA_Cntl_sub_lu( cntl ) );

    // A21 = A21 - A20 * A01
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A20, A01, FLA_ONE, A21,
                       FLA_Cntl_sub_gemm2( cntl ) );

    // A21 = A21 / triu( A11 )
    FLA_Trsm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, A11, A21,
                       FLA_Cntl_sub_trsm2( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
  }

  if ( FLA_Obj_width( ATR ) > 0 )
    /* ATR = trilu( ATL ) \ ATR */ 
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, ATL, ATR );

  return FLA_SUCCESS;
}

#endif
