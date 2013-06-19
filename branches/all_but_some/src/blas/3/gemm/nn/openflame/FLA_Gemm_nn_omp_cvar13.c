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

FLA_Error FLA_Gemm_nn_omp_var13( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  FLA_Obj BL,    BR,       B0,  B1,  B2;

  FLA_Obj C1L,   C1R,      C10,  C11,  C12;

  int b_m, b_n;
 

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_TOP );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_TOP );


  #pragma intel omp parallel taskq
  {
  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) )
  {
    b_m = FLA_Determine_blocksize( A, AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                        /* ** */            /* ** */
                                              &A1, 
                           AB,                &A2,        b_m, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                        /* ** */            /* ** */
                                              &C1, 
                           CB,                &C2,        b_m, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

    FLA_Part_1x2( C1,    &C1L,  &C1R,      0, FLA_LEFT );

    while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

      b_n = FLA_Determine_blocksize( B, BL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

      FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                             b_n, FLA_RIGHT );

      FLA_Repart_1x2_to_1x3( C1L,  /**/ C1R,        &C10, /**/ &C11, &C12,
                             b_n, FLA_RIGHT );

      /*------------------------------------------------------------*/

      #pragma intel omp task captureprivate( A1, B1, C11 )
      {
      /*    C1 = alpha * A * B1 + C1; */
      FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                         alpha, A1, B1, FLA_ONE, C11 );
      }

      /*------------------------------------------------------------*/

      FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                                FLA_LEFT );

      FLA_Cont_with_1x3_to_1x2( &C1L,  /**/ &C1R,        C10, C11, /**/ C12,
                                FLA_LEFT );

    }


    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                                                  A1, 
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                                                  C1, 
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );

  }
  }

  return FLA_SUCCESS;
}

