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

FLA_Error FLA_Lyap_n_blk_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;
  FLA_Obj CTL,   CTR,      C00, C01, C02, 
          CBL,   CBR,      C10, C11, C12,
                           C20, C21, C22;
  dim_t   b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_BR );

  while ( FLA_Obj_length( CTL ) > 0 ){

    b = FLA_Determine_blocksize( CTL, FLA_TL, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                &A10, &A11, /**/ &A12,
                        /* ************* */   /* ******************** */
                           ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                           b, b, FLA_TL );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00, &C01, /**/ &C02,
                                                &C10, &C11, /**/ &C12,
                        /* ************* */   /* ******************** */
                           CBL, /**/ CBR,       &C20, &C21, /**/ &C22,
                           b, b, FLA_TL );

    /*------------------------------------------------------------*/

    // C12 = isgn * C12 - A12 * C22;
    // C12 = sylv( A11, A22', C12 );
    FLA_Hemm_internal( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                       FLA_MINUS_ONE, C22, A12, isgn, C12,
                       FLA_Cntl_sub_hemm( cntl ) );
    FLA_Sylv_internal( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                       FLA_ONE, A11, A22, C12, scale,
                       FLA_Cntl_sub_sylv( cntl ) );

    // C11 = isgn * C11 - A12 * C12' - C12 * A12';
    FLA_Her2k_internal( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE,
                        FLA_MINUS_ONE, A12, C12, isgn, C11,
                        FLA_Cntl_sub_her2k( cntl ) );

    // C11 = lyap_n( A11, C11 );
    FLA_Lyap_internal( FLA_NO_TRANSPOSE, FLA_ONE, A11, C11, scale,
                       FLA_Cntl_sub_lyap( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                            /* ************** */  /* ****************** */
                                                     A10, /**/ A11, A12,
                              &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                              FLA_BR );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00, /**/ C01, C02,
                            /* ************** */  /* ****************** */
                                                     C10, /**/ C11, C12,
                              &CBL, /**/ &CBR,       C20, /**/ C21, C22,
                              FLA_BR );
  }

  return FLA_SUCCESS;
}

