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

FLA_Error FLA_Eig_gest_iu_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  //FLA_Obj yL,    yR,       y10t, psi11,  y12t;

  //FLA_Obj y12t_t,
  //        y12t_b;

  FLA_Obj psi11, y12t,
          y21,   Y22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x2( B,    &BTL, &BTR,
                      &BBL, &BBR,     0, 0, FLA_TL );

  //FLA_Part_1x2( Y,    &yL,  &yR,      0, FLA_LEFT );

  FLA_Part_2x2( Y,    &psi11, &y12t,
                      &y21,   &Y22,     1, 1, FLA_TL );

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

    //FLA_Repart_1x2_to_1x3( yL,  /**/ yR,        &y10t, /**/ &psi11,  &y12t,
    //                       1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    //FLA_Part_2x1( y12t,   &y12t_t,
    //                      &y12t_b,    1, FLA_TOP );

    // a01 = a01 * inv(beta11);
    FLA_Inv_scal_external( beta11, a01 );

    // A02 = A02 - a01 * b12t;
    FLA_Ger_external( FLA_MINUS_ONE, a01, b12t, A02 );

    // alpha11 = inv(beta11) * alpha11 * inv(conj(beta11));
    //         = inv(beta11) * alpha11 * inv(beta11);
    FLA_Inv_scal_external( beta11, alpha11 );
    FLA_Inv_scal_external( beta11, alpha11 );

    //// y12t = alpha11 * b12t;
    //FLA_Copy_external( b12t, y12t_t );
    //FLA_Scal_external( alpha11, y12t_t );
    // psi11 = - 1/2 * alpha11;
    FLA_Copy_external( alpha11, psi11 );
    FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );

    // a12t = inv(conj(beta11)) * a12t;
    //      = inv(beta11) * a12t;
    FLA_Inv_scal_external( beta11, a12t );

    //// a12t = a12t - 1/2 * y12t;
    //FLA_Axpy_external( FLA_MINUS_ONE_HALF, y12t_t, a12t );
    // a12t = a12t - 1/2 * alpha11 * b12t;
    FLA_Axpy_external( psi11, b12t, a12t );

    // A22 = A22 - a12t' * b12t - b12t' * a12t;
    FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
                        FLA_MINUS_ONE, a12t, b12t, A22 );

    //// a12t = a12t - 1/2 * y12t;
    //FLA_Axpy_external( FLA_MINUS_ONE_HALF, y12t_t, a12t );
    // a12t = a12t - 1/2 * alpha11 * b12t;
    FLA_Axpy_external( psi11, b12t, a12t );

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

    //FLA_Cont_with_1x3_to_1x2( &yL,  /**/ &yR,        y10t, psi11, /**/  y12t,
    //                          FLA_LEFT );
  }

  return FLA_SUCCESS;
}

