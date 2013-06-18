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

FLA_Error FLA_Gemv_n_blk_var2( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;

  FLA_Obj yT,              y0,
          yB,              y1,
                           y2;

  dim_t b;

  FLA_Part_2x1( A,    &AT, 
                      &AB,            0, FLA_BOTTOM );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( AB ) < FLA_Obj_length( A ) ){

    b = FLA_Determine_blocksize( AT, FLA_TOP, FLA_Cntl_blocksize( cntl ) );

    FLA_Repart_2x1_to_3x1( AT,                &A0, 
                                              &A1, 
                        /* ** */            /* ** */
                           AB,                &A2,        b, FLA_TOP );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                                              &y1, 
                        /* ** */            /* ** */
                           yB,                &y2,        b, FLA_TOP );

    /*------------------------------------------------------------*/

    /*   y1 = alpha * A1 * x + y1;   */
    FLA_Gemv_internal( FLA_NO_TRANSPOSE, 
                       alpha, A1, x, beta, y1, 
                       FLA_Cntl_sub_gemv( cntl ) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0, 
                            /* ** */           /* ** */
                                                  A1, 
                              &AB,                A2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                            /* ** */           /* ** */
                                                  y1, 
                              &yB,                y2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}

#endif
