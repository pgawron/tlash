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

FLA_Error FLA_UDdate_UT_unb_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T )
{
  FLA_Obj RTL,   RTR,      R00,  r01,   R02, 
          RBL,   RBR,      r10t, rho11, r12t,
                           R20,  r21,   R22;

  FLA_Obj CL,    CR,       C0,  c1,  C2;

  FLA_Obj DL,    DR,       D0,  d1,  D2;

  FLA_Obj TTL,   TTR,      T00,  t01,   T02,
          TBL,   TBR,      t10t, tau11, w12t,
                           T20,  t21,   T22;

  FLA_Part_2x2( R,    &RTL, &RTR,
                      &RBL, &RBR,     0, 0, FLA_TL );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_2x2( T,    &TTL, &TTR,
                      &TBL, &TBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( RBR ) > 0 ){

    FLA_Repart_2x2_to_3x3( RTL, /**/ RTR,       &R00,  /**/ &r01,   &R02,
                        /* ************* */   /* ************************** */
                                                &r10t, /**/ &rho11, &r12t,
                           RBL, /**/ RBR,       &R20,  /**/ &r21,   &R22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &c1, &C2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &d1, &D2,
                           1, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00,  /**/ &t01,   &T02,
                        /* ************* */   /* ************************ */
                                                &t10t, /**/ &tau11, &w12t,
                           TBL, /**/ TBR,       &T20,  /**/ &t21,   &T22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/

    // Compute tau11, u1, and v2 from rho11, c1, and d1 such that tau11, u1,
    // and v1 determine an up/downdating UT Householder transform H such that
    // applying H from the left to the column vector consisting of rho11, c1,
    // and d1 annihilates the entries in c1 and d1 (and updates rho11).
    FLA_Househ3UD_UT( rho11,
                      c1,
                      d1, tau11 );

    // / r12t \       / r12t \ 
    // |  C2  | =  H' |  C2  | 
    // \  D2  /       \  D2  / 
    //
    // where H is formed from tau11, u1 (stored in c1) and v1 (stored in d1).
    FLA_Apply_HUD_UT( FLA_LEFT,
                      tau11, w12t,
                             r12t,
                      c1,    C2,
                      d1,    D2 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &RTL, /**/ &RTR,       R00,  r01,   /**/ R02,
                                                     r10t, rho11, /**/ r12t,
                            /* ************** */  /* ************************ */
                              &RBL, /**/ &RBR,       R20,  r21,   /**/ R22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, c1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, d1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00,  t01,   /**/ T02,
                                                     t10t, tau11, /**/ w12t,
                            /* ************** */  /* ********************** */
                              &TBL, /**/ &TBR,       T20,  t21,   /**/ T22,
                              FLA_TL );
  }

  // T = I + C' * C - D' * D;
  // T = striu( T ) + 0.5*diag( T );
  // NOTE: The only reason this 'herk' method of computing T works is because
  // up-and-downdating is used to up/downdate a system that is being solved
  // either by QR factorization, or the method of normal equations (Cholesky
  // factorization on A' * A), and in either case, R will have a real diagonal.

  FLA_Set_to_identity( T );

  FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                     FLA_ONE, C, FLA_ONE, T );

  FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
                     FLA_MINUS_ONE, D, FLA_ONE, T );

  FLA_Scale_diag( FLA_NO_CONJUGATE, FLA_ONE_HALF, T );

  return FLA_SUCCESS;
}

