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

FLA_Error FLA_QR2_UT_blk_var2( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl )
{
  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  dim_t b;

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_TOP );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_TOP );

  while ( FLA_Obj_length( DT ) < FLA_Obj_length( D ) ){

    b = FLA_Determine_blocksize( DB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                        /* ** */            /* ****** */
                                              &D1, 
                           DB,                &D2,        b, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                        /* ** */            /* ****** */
                                              &T1, 
                           TB,                &T2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    // [ U, ...
    //   D1, T ] = FLA_QR2_UT( U
    //                         D1, T1 );

    FLA_QR2_UT_internal( U,
                         D1, T1, 
                         FLA_Cntl_sub_qr2ut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                                                  D1, 
                            /* ** */           /* ****** */
                              &DB,                D2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                                                  T1, 
                            /* ** */           /* ****** */
                              &TB,                T2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

