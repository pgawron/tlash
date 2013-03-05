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

FLA_Error FLASH_LQ_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X )
{
  FLA_Obj W;
  FLA_Obj AL, AR;
  FLA_Obj XT, XB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LQ_UT_solve_check( A, T, B, X );

  FLASH_Apply_Q_UT_create_workspace( T, X, &W );

  FLA_Part_1x2( A,   &AL, &AR,    FLA_Obj_length( A ), FLA_LEFT );
  FLA_Part_2x1( X,   &XT,
                     &XB,    FLA_Obj_length( B ), FLA_TOP );

  FLASH_Copy( B, XT );

  FLASH_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE,
              FLA_NONUNIT_DIAG, FLA_ONE, AL, XT );

  FLASH_Set( FLA_ZERO, XB );

  FLASH_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
                    A, T, W, X );

  FLASH_Obj_free( &W );

  return FLA_SUCCESS;
}

