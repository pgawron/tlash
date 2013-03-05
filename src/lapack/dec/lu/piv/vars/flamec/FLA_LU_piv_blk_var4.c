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

FLA_Error FLA_LU_piv_blk_var4( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj pT,              p0,
          pB,              p1,
                           p2;

  FLA_Obj AB0, AB1, AB2;

  dim_t b;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
          FLA_Obj_width( ATL ) < FLA_Obj_width( A )){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* ** */
                                              &p1, 
                           pB,                &p2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // A11 = A11 - A10 * A0
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A01, FLA_ONE, A11,
                       FLA_Cntl_sub_gemm1( cntl ) );

    // A21 = A21 - A20 * A01
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A20, A01, FLA_ONE, A21,
                       FLA_Cntl_sub_gemm3( cntl ) );

    // AB1 = / A11 \
    //       \ A21 /
    FLA_Merge_2x1( A11,
                   A21,      &AB1 );

    // AB1, p1 = LU_piv( AB1 )
    FLA_LU_piv_internal( AB1, p1, 
                         FLA_Cntl_sub_lu( cntl ) );

    // AB0 = / A10 \
    //       \ A20 /
    FLA_Merge_2x1( A10,
                   A20,      &AB0 );

    // AB2 = / A12 \
    //       \ A22 /
    FLA_Merge_2x1( A12,
                   A22,      &AB2 );

    // Apply pivots to remaining columns
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB0,
                               FLA_Cntl_sub_appiv1( cntl ) );
    FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB2,
                               FLA_Cntl_sub_appiv1( cntl ) );

    // A12 = A12 - A10 * A02
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                       FLA_MINUS_ONE, A10, A02, FLA_ONE, A12,
                       FLA_Cntl_sub_gemm2( cntl ) );

    // A12 = trilu( A11 ) \ A12 
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, A11, A12,
                       FLA_Cntl_sub_trsm1( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  p1, 
                            /* ** */           /* ** */
                              &pB,                p2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

#endif
