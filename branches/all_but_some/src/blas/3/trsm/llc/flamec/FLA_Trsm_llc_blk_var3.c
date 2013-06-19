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

FLA_Error FLA_Trsm_llc_blk_var3( FLA_Diag diagA, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  dim_t b;

  FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

  while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) ){

    b = FLA_Determine_blocksize( BR, FLA_RIGHT, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    /* B1 = tril( A ) * B1; */
    FLA_Trsm_internal( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diagA,
                       alpha, A, B1,
                       FLA_Cntl_sub_trsm( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}

