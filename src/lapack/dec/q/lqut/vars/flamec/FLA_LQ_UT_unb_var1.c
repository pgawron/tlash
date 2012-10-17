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

FLA_Error FLA_LQ_UT_unb_var1( FLA_Obj A, FLA_Obj t )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj tLt,   tRt,      t0t,  tau1,  t2t;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( t,    &tLt,  &tRt,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_1x2_to_1x3( tLt,  /**/ tRt,      &t0t, /**/ &tau1, &t2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // Compute tau11 and u12t from alpha11 and a12t such that tau11 and u12t
    // determine a Householder transform H such that applying H from the
    // right to the row vector consisting of alpha11 and a12t annihilates
    // the entries in a12t (and updates alpha11).
    FLA_Househ2_UT( FLA_RIGHT, alpha11, a12t,
                    tau1 );

    // ( a21 A22 ) = ( a21 A22 ) H
    //
    // where H is formed from tau11 and u12t.
    FLA_Apply_H2_UT( FLA_RIGHT, tau1, a12t, a21, A22 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &tLt,  /**/ &tRt,      t0t, tau1, /**/ t2t,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

