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

FLA_Error FLA_Hess_UT_recover_tau_submatrix( FLA_Obj T, FLA_Obj t );

FLA_Error FLA_Hess_UT_recover_tau( FLA_Obj T, FLA_Obj t )
{
  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj tT,              t0,
          tB,              t1,
                           t2;

  dim_t b_alg, b;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Hess_UT_recover_tau_check( T, t );

  b_alg = FLA_Obj_length( T );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_2x1( t,    &tT, 
                      &tB,            0, FLA_TOP );

  while ( FLA_Obj_length( tT ) < FLA_Obj_length( t ) ){

    b = min( FLA_Obj_length( tB ), b_alg );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( tT,                &t0, 
                        /* ** */            /* ** */
                                              &t1, 
                           tB,                &t2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Hess_UT_recover_tau_submatrix( T1, t1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &tT,                t0, 
                                                  t1, 
                            /* ** */           /* ** */
                              &tB,                t2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}


FLA_Error FLA_Hess_UT_recover_tau_submatrix( FLA_Obj T, FLA_Obj t )
{
  FLA_Obj TTL,   TTR,      T00,  t01,   T02, 
          TBL,   TBR,      t10t, tau11, t12t,
                           T20,  t21,   T22;

  FLA_Obj tT,              t0,
          tB,              tau1,
                           t2;


  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  FLA_Part_2x1( t,    &tT, 
                      &tB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( TBR ) > 0 ){

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &t12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( tT,                  &t0, 
                        /* ** */              /* ** */
                                                &tau1, 
                           tB,                  &t2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // tau1 = tau11;
    FLA_Copy_external( tau11, tau1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ t12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &tT,                   t0, 
                                                     tau1, 
                            /* ** */              /* ** */
                              &tB,                   t2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}

