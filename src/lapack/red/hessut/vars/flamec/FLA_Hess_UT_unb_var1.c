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

FLA_Error FLA_Hess_UT_unb_var1( FLA_Obj A, FLA_Obj T )
{
  return FLA_Hess_UT_step_unb_var1( A, T );
}

FLA_Error FLA_Hess_UT_step_unb_var1( FLA_Obj A, FLA_Obj T )
{
  FLA_Obj  ATL,   ATR,      A00,  a01,     A02, 
           ABL,   ABR,      a10t, alpha11, a12t,
                            A20,  a21,     A22;
  FLA_Obj  AL,    AR,       A0,   a1,      A2;
  FLA_Obj  TTL,   TTR,      T00,  t01,   T02, 
           TBL,   TBR,      t10t, tau11, t12t,
                            T20,  t21,   T22;
           
  FLA_Obj  a21_t,
           a21_b;

  FLA_Obj  A22_t,
           A22_b;

  FLA_Obj  A2_l, A2_r;

  FLA_Obj  first_elem;

  dim_t        b_alg;
  FLA_Datatype datatype_A;


  b_alg      = FLA_Obj_length( T );
  datatype_A = FLA_Obj_datatype( A );

  FLA_Obj_create( datatype_A, 1, 1, 0, 0, &first_elem );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_1x2( A,    &AL,  &AR,         0, FLA_LEFT );
  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < b_alg )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************** */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_length( A22 ) > 0 )
    {
      FLA_Part_2x1( a21,    &a21_t,
                            &a21_b,        1, FLA_TOP );

      FLA_Part_2x1( A22,    &A22_t,
                            &A22_b,        1, FLA_TOP );

      FLA_Part_1x2( A2,     &A2_l, &A2_r,        1, FLA_LEFT );

      // [ u21, tau11, a21 ] = House( a21 );
      FLA_Househ2_UT( FLA_LEFT,
                      a21_t,
                      a21_b, tau11 );

      // Save first element of a21_t and set it to one so we can use a21 as
      // u21 in subsequent computations. We will restore a21_t later on.
      FLA_Copy( a21_t, first_elem );
      FLA_Set( FLA_ONE, a21_t );

      // A22 = ( I - inv( tau ) * u21 * u21' ) * A22;
      FLA_Apply_H2_UT( FLA_LEFT, tau11, a21_b, A22_t,
                                               A22_b );

      // A02  = A02  * ( I - inv( tau ) * u21 * u21' );
      // a12t = a12t * ( I - inv( tau ) * u21 * u21' );
      // A22  = A22  * ( I - inv( tau ) * u21 * u21' );
      FLA_Apply_H2_UT( FLA_RIGHT, tau11, a21_b, A2_l, A2_r );

      // t01 = U20' * u21;
      FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );

      // Restore first element of a21.
      FLA_Copy( first_elem, a21_t );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ************************ */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
  }

  FLA_Obj_free( &first_elem );

  return FLA_SUCCESS;
}

