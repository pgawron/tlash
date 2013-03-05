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

FLA_Error FLA_Apply_CAQ2_UT_lhfc_blk_var3( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, 
                                                                            FLA_Obj E, fla_apcaq2ut_t* cntl )
{
  FLA_Obj WL,    WR,       W0,  W1,  W2;

  FLA_Obj CL,    CR,       C0,  C1,  C2;

  FLA_Obj EL,    ER,       E0,  E1,  E2;

  dim_t b;

  FLA_Part_1x2( W,    &WL,  &WR,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( E,    &EL,  &ER,      0, FLA_LEFT );

  while ( FLA_Obj_width( CL ) < FLA_Obj_width( C ) ){

    b = FLA_Determine_blocksize( CR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( WL,  /**/ WR,        &W0, /**/ &W1, &W2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( EL,  /**/ ER,        &E0, /**/ &E1, &E2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    //  / C1 \ =  Q' / C1 \
    //  \ E1 /       \ E1 /
    //
    //  where Q is formed from D and T.

    FLA_Apply_CAQ2_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                                D, T, W1, C1,
                                          E1, FLA_Cntl_sub_apcaq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &WL,  /**/ &WR,        W0, W1, /**/ W2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &EL,  /**/ &ER,        E0, E1, /**/ E2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

