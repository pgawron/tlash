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

FLA_Error FLA_Gemm_tn_blk_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  FLA_Obj AL,    AR,       A0,  A1,  A2;

  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;

  dim_t b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_RIGHT );

  FLA_Part_2x1( C,    &CT, 
                      &CB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width( AR ) < FLA_Obj_width( A ) ){

    b = FLA_Determine_blocksize( AL, FLA_LEFT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, &A1, /**/ &A2,
                           b, FLA_LEFT );

    FLA_Repart_2x1_to_3x1( CT,                &C0, 
                                              &C1, 
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /* C1 = alpha * A1' * B + C1; */
    FLA_Gemm_internal( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
                       alpha, A1, B, beta, C1,
                       FLA_Cntl_sub_gemm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, /**/ A1, A2,
                              FLA_RIGHT );

    FLA_Cont_with_3x1_to_2x1( &CT,                C0, 
                            /* ** */           /* ** */
                                                  C1, 
                              &CB,                C2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
