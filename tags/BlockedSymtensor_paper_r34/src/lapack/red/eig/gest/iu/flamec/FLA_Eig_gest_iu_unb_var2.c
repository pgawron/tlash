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

FLA_Error FLA_Eig_gest_iu_unb_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  FLA_Obj yT,              y01,
          yB,              psi11,
                           y21;

  FLA_Obj y01_l, y01_r;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  FLA_Part_2x1( Y,    &yT,
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x2_to_3x3( BTL, /**/ BTR,       &B00,  /**/ &b01,    &B02,
                        /* ************* */   /* ************************* */
                                                &b10t, /**/ &beta11, &b12t,
                           BBL, /**/ BBR,       &B20,  /**/ &b21,    &B22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( yT,                  &y01,
                        /* ** */              /* ***** */
                                                &psi11,
                           yB,                  &y21,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( y01,    &y01_l, &y01_r,     1, FLA_LEFT );

    // y01 = 1/2 * A00 * b01;
    FLA_Hemvc_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
                        FLA_ONE_HALF, A00, b01, FLA_ZERO, y01_l );

    // a01 = a01 - y01;
    FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );

    // alpha11 = alpha11 - a01' * b01 - b01' * a01;
    FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );

    // alpha11 = inv(beta11) * alpha11 * inv(conj(beta11));
    //         = inv(beta11) * alpha11 * inv(beta11);
    FLA_Inv_scal_external( beta11, alpha11 );
    FLA_Inv_scal_external( beta11, alpha11 );

    // a12t   = a12t   - b01' * A02; 
    // a12t^T = a12t^T - A02^T * conj(b01); 
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
                        FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );

    // a12t = inv(conj(beta11)) * a12t;
    // a12t = inv(beta11) * a12t;
    FLA_Inv_scal_external( beta11, a12t );

    // a01 = a01 - y01;
    FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );

    // a01 = a01 * inv(beta11);
    FLA_Inv_scal_external( beta11, a01 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x3_to_2x2( &BTL, /**/ &BTR,       B00,  b01,    /**/ B02,
                                                     b10t, beta11, /**/ b12t,
                            /* ************** */  /* *********************** */
                              &BBL, /**/ &BBR,       B20,  b21,    /**/ B22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &yT,                   y01,
                                                     psi11,
                            /* ** */              /* ***** */
                              &yB,                   y21,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

