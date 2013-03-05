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

FLA_Error FLA_Apply_Q2_UT_lnfc_blk_var2( FLA_Obj D, FLA_Obj T, FLA_Obj W1, FLA_Obj C, 
                                                                           FLA_Obj E, fla_apq2ut_t* cntl )
{
  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  FLA_Obj ET,              E0,
          EB,              E1,
                           E2;

  dim_t b;

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_BOTTOM );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_BOTTOM );

  FLA_Part_2x1( E,    &ET, 
                      &EB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( DB ) < FLA_Obj_length( D ) ){

    b = FLA_Determine_blocksize( DT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                                              &D1, 
                        /* ** */            /* ** */
                           DB,                &D2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                                              &T1, 
                        /* ** */            /* ** */
                           TB,                &T2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( ET,                &E0, 
                                              &E1, 
                        /* ** */            /* ** */
                           EB,                &E2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    //  / C  \ =  Q / C  \
    //  \ E1 /      \ E1 /
    //
    //  where Q is formed from D1 and T1.

    FLA_Apply_Q2_UT_internal( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                              D1, T1, W1, C,
                                          E1, FLA_Cntl_sub_apq2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                            /* ** */           /* ** */
                                                  D1, 
                              &DB,                D2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                            /* ** */           /* ** */
                                                  T1, 
                              &TB,                T2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &ET,                E0, 
                            /* ** */           /* ** */
                                                  E1, 
                              &EB,                E2,     FLA_BOTTOM );
  }

  return FLA_SUCCESS;
}

