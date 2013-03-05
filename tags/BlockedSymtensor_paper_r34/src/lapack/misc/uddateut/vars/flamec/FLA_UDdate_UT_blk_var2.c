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

FLA_Error FLA_UDdate_UT_blk_var2( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl )
{
  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  dim_t   b_C, b_D, b_T;

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_TOP );

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_TOP );

  while ( FLA_Obj_length( CT ) < FLA_Obj_length( C ) &&
          FLA_Obj_length( DT ) < FLA_Obj_length( D ) ){

    b_C = FLA_Determine_blocksize( CB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );
    b_D = FLA_Determine_blocksize( DB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );
    b_T = FLA_Determine_blocksize( TB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */          /* ****** */
                                              &C1, 
                           CB,                &C2,        b_C, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                        /* ** */          /* ****** */
                                              &D1, 
                           DB,                &D2,        b_D, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                        /* ** */          /* ****** */
                                              &T1, 
                           TB,                &T2,        b_T, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /*
       Perform an up/downdate of the upper triangular Cholesky factor R via
       "UD" UT Householder transformations:

         [ R, ...
           C1, ...
           D1, T1 ] = FLA_UDdate_UT( R, ...
                                     C1, ...
                                     D1, T1 );

       by updating R in such a way that removes the contributions of the rows
       in D1 while simultaneously adding new contributions to the factorization
       from the rows of C1. Note that C1 and D1 are also updated in the process.
       Also note that either C1 or D1 may become empty at any iteration.
    */

    FLA_UDdate_UT_internal( R,
                            C1,
                            D1, T1,
                            FLA_Cntl_sub_uddateut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */          /* ****** */
                              &CB,                C2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                                                  D1, 
                            /* ** */          /* ****** */
                              &DB,                D2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                                                  T1, 
                            /* ** */          /* ****** */
                              &TB,                T2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

