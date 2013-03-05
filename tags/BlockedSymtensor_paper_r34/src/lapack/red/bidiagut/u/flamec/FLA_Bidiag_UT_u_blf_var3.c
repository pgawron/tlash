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

FLA_Error FLA_Bidiag_UT_u_blf_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Obj  ATL,   ATR,      A00, A01, A02, 
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  TUL,   TUR,      TU0, TU1, TU2; 
  FLA_Obj  TVL,   TVR,      TV0, TV1, TV2; 

  FLA_Obj  TU1_tl;
  FLA_Obj  TV1_tl;
  FLA_Obj  none, none2, none3;
  dim_t    b_alg, b;

  b_alg = FLA_Obj_length( TU );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,   0, 0, FLA_TL );
  FLA_Part_1x2( TU,   &TUL, &TUR,      0, FLA_LEFT ); 
  FLA_Part_1x2( TV,   &TVL, &TVR,      0, FLA_LEFT ); 

  while ( FLA_Obj_min_dim( ABR ) > 0 )
  {
    b = min( FLA_Obj_min_dim( ABR ), b_alg );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( TUL, /**/ TUR,       &TU0, /**/ &TU1, &TU2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( TVL, /**/ TVR,       &TV0, /**/ &TV1, &TV2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x2( TU1,     &TU1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    FLA_Part_2x2( TV1,     &TV1_tl, &none,   
                           &none2,  &none3,   b, b, FLA_TL ); 

    // [ ABR, T1 ] = FLA_Bidiag_UT_u_step_unb_var3( ABR, TU1, TV1, b );
    //FLA_Bidiag_UT_u_step_unb_var3( ABR, TU1_tl, TV1_tl );
    FLA_Bidiag_UT_u_step_ofu_var3( ABR, TU1_tl, TV1_tl );
    //FLA_Bidiag_UT_u_step_opt_var3( ABR, TU1_tl, TV1_tl );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &TUL, /**/ &TUR,       TU0, TU1, /**/ TU2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &TVL, /**/ &TVR,       TV0, TV1, /**/ TV2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

