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

FLA_Error FLA_Apply_CAQ_UT_inc_lhfc_blk_var1( FLA_Obj A, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apcaqutinc_t* cntl )
{
  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TTL,   WTR,      T00, W01, W02, 
          TBL,   TBR,      T10, T11, W12,
                           T20, T21, T22;

  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;

  dim_t   b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( TW,   &TTL, &WTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  FLA_Part_2x1( B,    &BT,
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = FLA_Determine_blocksize( ABR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_2x2_to_3x3( TTL, /**/ WTR,       &T00, /**/ &W01, &W02,
                        /* ************* */   /* ******************** */
                                                &T10, /**/ &T11, &W12,
                           TBL, /**/ TBR,       &T20, /**/ &T21, &T22,
                           b, b, FLA_BR );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /*
        Apply Q^H to B1 and B2 from the left:
      
            / B1 \ = Q^H * / B1 \
            \ B2 /         \ B2 / 
      
        where Q is formed from A21 and T21.
    */

    FLA_Apply_CAQ2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                A21, T21, W, B1,
                                             B2,
                                FLA_Cntl_sub_apcaq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &WTR,       T00, W01, /**/ W02,
                                                     T10, T11, /**/ W12,
                            /* ************** */  /* ****************** */
                              &TBL, /**/ &TBR,       T20, T21, /**/ T22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

