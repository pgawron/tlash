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

FLA_Error FLA_Eig_gest_nl_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj BTL,   BTR,      B00,  b01,    B02, 
          BBL,   BBR,      b10t, beta11, b12t,
                           B20,  b21,    B22;

  //FLA_Obj yL,    yR,       y10t, psi11,  y12t;

  //FLA_Obj y10t_t,
  //        y10t_b;

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

    //FLA_Part_2x1( y10t,   &y10t_t,
    //                      &y10t_b,    1, FLA_TOP );

    //// y10t = alpha11 * b10t;
    //FLA_Copy_external( b10t, y10t_t );
    //FLA_Scal_external( alpha11, y10t_t );
    // psi11 = 1/2 * alpha11;
    FLA_Copy_external( alpha11, psi11 );
    FLA_Scal_external( FLA_ONE_HALF, psi11 );

    //// a10t = a10t + 1/2 * y10t;
    //FLA_Axpy_external( FLA_ONE_HALF, y10t_t, a10t );
    // a10t = a10t + 1/2 * alpha11 * b10t;
    FLA_Axpy_external( psi11, b10t, a10t );

    // A00 = A00 + a10t' * b10t + b10t' * a10t;
    FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
                        FLA_ONE, a10t, b10t, A00 );

    //// a10t = a10t + 1/2 * y10t;
    //FLA_Axpy_external( FLA_ONE_HALF, y10t_t, a10t );
    // a10t = a10t + 1/2 * alpha11 * b10t;
    FLA_Axpy_external( psi11, b10t, a10t );

    // a10t = conj(beta11) * a10t;
    //      = beta11 * a10t;
    FLA_Scal_external( beta11, a10t );

    // alpha11 = conj(beta11) * alpha11 * beta11;
    //         = beta11 * alpha11 * beta11;
    FLA_Scal_external( beta11, alpha11 );
    FLA_Scal_external( beta11, alpha11 );

    // A20 = A20 + a21 * b10t;
    FLA_Ger_external( FLA_ONE, a21, b10t, A20 );

    // a21 = a21 * beta11;
    FLA_Scal_external( beta11, a21 );

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

#endif