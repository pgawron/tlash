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

FLA_Error FLA_LU_piv_unb_var3b( FLA_Obj A, FLA_Obj p )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj pT,              p0,
          pB,              pi1,
                           p2;

  FLA_Obj AB0, aB1, AB2;


  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( p,    &pT, 
                      &pB,            0, FLA_TOP );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
          FLA_Obj_width( ATL ) < FLA_Obj_width( A )){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( pT,                &p0, 
                        /* ** */            /* *** */
                                              &pi1, 
                           pB,                &p2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // a01 = trilu( A00 ) \ a01
    FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );

    // alpha11 = alpha11 - a10t * a01
    FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );

    // a21 = a21 - A20 * a01
    FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );

    // aB1 = / alpha11 \
    //       \  a21    /
    FLA_Merge_2x1( alpha11,
                   a21,      &aB1 );

    // Determine pivot index
    FLA_Amax_external( aB1, pi1 );

    // Apply pivots to current column
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );

    // AB0 = / a10t \
    //       \ A20  /
    FLA_Merge_2x1( a10t,
                   A20,      &AB0 );

    // AB2 = / a12t \
    //       \ A22  /
    FLA_Merge_2x1( a12t,
                   A22,      &AB2 );

    // Apply pivots to remaining columns
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
    FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );

    // a21 = a21 / alpha11
    FLA_Inv_scal_external( alpha11, a21 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &pT,                p0, 
                                                  pi1, 
                            /* ** */           /* *** */
                              &pB,                p2,     FLA_TOP );

  }

  if ( FLA_Obj_width( ATR ) > 0 )
  {
    /* ATR = trilu( ATL ) \ ATR */ 
    FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
                       FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                       FLA_ONE, ATL, ATR );
  }

  return FLA_SUCCESS;
}

#endif
