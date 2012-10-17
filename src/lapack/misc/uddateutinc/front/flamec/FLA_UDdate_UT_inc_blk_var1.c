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

FLA_Error FLA_UDdate_UT_inc_blk_var1( FLA_Obj R,
                                      FLA_Obj C,
                                      FLA_Obj D, FLA_Obj T, FLA_Obj W, fla_uddateutinc_t* cntl )
{
  FLA_Obj RTL,   RTR,      R00, R01, R02, 
          RBL,   RBR,      R10, R11, R12,
                           R20, R21, R22;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj DL,    DR,       D0,  D1,  D2;

  FLA_Obj TL,    TR,       T0,  T1,  T2;

  FLA_Obj WTL,   WTR,      W00, W01, W02, 
          WBL,   WBR,      W10, W11, W12,
                           W20, W21, W22;

  dim_t b;

  FLA_Part_2x2( R,    &RTL, &RTR,
                      &RBL, &RBR,     0, 0, FLA_TL );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  FLA_Part_2x2( W,    &WTL, &WTR,
                      &WBL, &WBR,     0, 0, FLA_TL );

  while ( FLA_Obj_min_dim( RBR ) > 0 ){

    b = FLA_Determine_blocksize( RBR, FLA_BR, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x2_to_3x3( RTL, /**/ RTR,       &R00, /**/ &R01, &R02,
                        /* ************* */   /* ******************** */
                                                &R10, /**/ &R11, &R12,
                           RBL, /**/ RBR,       &R20, /**/ &R21, &R22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &T2,
                           b, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( WTL, /**/ WTR,       &W00, /**/ &W01, &W02,
                        /* ************* */   /* ******************** */
                                                &W10, /**/ &W11, &W12,
                           WBL, /**/ WBR,       &W20, /**/ &W21, &W22,
                           b, b, FLA_BR );

    /*------------------------------------------------------------*/

    /*
       Perform an up/downdate of the upper triangular factor R11 via 
       up/downdating UT Householder transformations:

         [ R11, ...
           C1,  ...
           D1, T1 ] = FLA_UDdate_UT( R11, ...
                                     C1, ...
                                     D1, T1 );

       by updating R11 in such a way that removes the contributions of
       the rows in D1 while simultaneously adding new contributions to
       the factorization from the rows of C1. Note that C1 and D1 are
       also updated in the process.
    */

    FLA_UDdate_UT_internal( R11,
                            C1,
                            D1, T1,
                            FLA_Cntl_sub_uddateut( cntl ) );



    if ( FLA_Obj_width( R12 ) > 0 )
    {
      /*
           Apply Q' to R12, C2, and D2 from the left:
  
             / R12 \          / R12 \
             | C2  |  =  Q' * | C2  |
             \ D2  /          \ D2  /
  
           where Q is formed from C1, D1, and T1.
      */

      FLA_Apply_QUD_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                 T1, W12,
                                     R12,
                                 C1, C2,
                                 D1, D2, FLA_Cntl_sub_apqudut( cntl ) );
    }


    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &RTL, /**/ &RTR,       R00, R01, /**/ R02,
                                                     R10, R11, /**/ R12,
                            /* ************** */  /* ****************** */
                              &RBL, /**/ &RBR,       R20, R21, /**/ R22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ T2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &WTL, /**/ &WTR,       W00, W01, /**/ W02,
                                                     W10, W11, /**/ W12,
                            /* ************** */  /* ****************** */
                              &WBL, /**/ &WBR,       W20, W21, /**/ W22,
                              FLA_TL );
  }

  return FLA_SUCCESS;
}

