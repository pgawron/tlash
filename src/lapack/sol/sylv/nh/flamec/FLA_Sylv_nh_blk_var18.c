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

FLA_Error FLA_Sylv_nh_blk_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
  FLA_Obj BTL,   BTR,      B00, B01, B02, 
          BBL,   BBR,      B10, B11, B12,
                           B20, B21, B22;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  dim_t b;

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_BR );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_RIGHT );

  while ( FLA_Obj_length( BBR ) < FLA_Obj_length( B ) ){

    b = FLA_Determine_blocksize( CL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00, &B01, /**/ &B02,
                                                &B10, &B11, /**/ &B12,
                        /* ************* */   /* ******************** */
                           BBL, /**/ BBR,       &B20, &B21, /**/ &B22,
                           b, b, FLA_TL );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, &C1, /**/ &C2,
                           b, FLA_LEFT );

    // Loop Invariant:
    // CL = 
    // CR = 

    /*------------------------------------------------------------*/

    // C1 = sylv( A, B11', C1 );
    FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       isgn, A, B11, C1, scale,
                       FLA_Cntl_sub_sylv1( cntl ) );

    // C0 = C0 -/+ C1 * B01';
    FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, 
                       FLA_NEGATE( isgn ), C1, B01, FLA_ONE, C0,
                       FLA_Cntl_sub_gemm1( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00, /**/ B01, B02,
                            /* ************** */  /* ****************** */
                                                     B10, /**/ B11, B12,
                              &BBL, /**/ &BBR,       B20, /**/ B21, B22,
                              FLA_BR );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, /**/ C1, C2,
                              FLA_RIGHT );
  }

  return FLA_SUCCESS;
}
