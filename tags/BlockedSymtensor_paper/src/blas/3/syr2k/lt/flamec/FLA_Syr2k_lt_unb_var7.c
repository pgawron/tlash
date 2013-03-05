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

FLA_Error FLA_Syr2k_lt_unb_var7( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj BL,    BR,       B0,  b1,  B2;

  FLA_Obj CTL,   CTR,      C00,  c01,     C02, 
          CBL,   CBR,      c10t, gamma11, c12t,
                           C20,  c21,     C22;


  FLA_Scalr_external( FLA_LOWER_TRIANGULAR, beta, C );

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_RIGHT );

  FLA_Part_2x2( C,    &CTL, &CTR,
                      &CBL, &CBR,     0, 0, FLA_BR );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){


    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &a1, /**/ &A2,
                           1, FLA_LEFT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, &b1, /**/ &B2,
                           1, FLA_LEFT );

    FLA_Repart_2x2_to_3x3( CTL, /**/ CTR,       &C00,  &c01,     /**/ &C02,
                                                &c10t, &gamma11, /**/ &c12t,
                        /* ************* */   /* ************************** */
                           CBL, /**/ CBR,       &C20,  &c21,     /**/ &C22,
                           1, 1, FLA_TL );

    /*------------------------------------------------------------*/

    /* c10t = c10t + b1' * A0 */
    FLA_Gemv_external( FLA_TRANSPOSE, alpha, A0, b1, FLA_ONE, c10t );

    /* c21 = c21 + A2' * b1 */
    FLA_Gemv_external( FLA_TRANSPOSE, alpha, A2, b1, FLA_ONE, c21 );    

    /* gamma11 = gamma11 + a1' * b1 + b1' * a1 */
    FLA_Dot2s_external( alpha, a1, b1, FLA_ONE, gamma11 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ a1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, /**/ b1, B2,
                              FLA_RIGHT );

    FLA_Cont_with_3x3_to_2x2( &CTL, /**/ &CTR,       C00,  /**/ c01,     C02,
                            /* ************** */  /* ************************ */
                                                     c10t, /**/ gamma11, c12t,
                              &CBL, /**/ &CBR,       C20,  /**/ c21,     C22,
                              FLA_BR );

  }

  return FLA_SUCCESS;
}

#endif
