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

FLA_Error FLA_Apply_QUD_UT_lhfc_blk_var3( FLA_Obj T, FLA_Obj W,
                                                     FLA_Obj R,
                                          FLA_Obj U, FLA_Obj C,
                                          FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl )
{
  FLA_Obj TT,              T0,
          TB,              T1,
                           T2;

  FLA_Obj UT,              U0,
          UB,              U1,
                           U2;

  FLA_Obj VT,              V0,
          VB,              V1,
                           V2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj DT,              D0,
          DB,              D1,
                           D2;

  dim_t b_T, b_UC, b_VD;

  FLA_Part_2x1( T,    &TT, 
                      &TB,            0, FLA_TOP );

  FLA_Part_2x1( U,    &UT, 
                      &UB,            0, FLA_TOP );

  FLA_Part_2x1( V,    &VT, 
                      &VB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );

  FLA_Part_2x1( D,    &DT, 
                      &DB,            0, FLA_TOP );

  while ( FLA_Obj_length( TT ) < FLA_Obj_length( T ) ){

    b_T  = FLA_Determine_blocksize( TB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );
    b_UC = FLA_Determine_blocksize( UB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );
    b_VD = FLA_Determine_blocksize( VB, FLA_BOTTOM, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( TT,                &T0, 
                        /* ** */            /* ** */
                                              &T1, 
                           TB,                &T2,        b_T, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( UT,                &U0, 
                        /* ** */            /* ** */
                                              &U1, 
                           UB,                &U2,        b_UC, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( VT,                &V0, 
                        /* ** */            /* ** */
                                              &V1, 
                           VB,                &V2,        b_VD, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b_UC, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( DT,                &D0, 
                        /* ** */            /* ** */
                                              &D1, 
                           DB,                &D2,        b_VD, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    //  / R  \        / R  \
    //  | C1 |  =  Q' | C1 |
    //  \ D1 /        \ D1 /
    //
    //  where Q is formed from U1, V1, and T1.

    FLA_Apply_QUD_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                               T1, W,
                                   R,
                               U1, C1,
                               V1, D1, FLA_Cntl_sub_apqudut( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &TT,                T0, 
                                                  T1, 
                            /* ** */           /* ** */
                              &TB,                T2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &UT,                U0, 
                                                  U1, 
                            /* ** */           /* ** */
                              &UB,                U2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &VT,                V0, 
                                                  V1, 
                            /* ** */           /* ** */
                              &VB,                V2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &DT,                D0, 
                                                  D1, 
                            /* ** */           /* ** */
                              &DB,                D2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

