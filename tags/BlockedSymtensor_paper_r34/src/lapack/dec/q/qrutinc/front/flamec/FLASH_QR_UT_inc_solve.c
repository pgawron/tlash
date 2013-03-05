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

FLA_Error FLASH_QR_UT_inc_solve( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X )
{
  FLA_Obj W, Y;
  FLA_Obj AT, AB;
  FLA_Obj YT, YB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_inc_solve_check( A, TW, B, X );

  FLASH_Apply_Q_UT_inc_create_workspace( TW, B, &W );

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &Y );

  FLASH_Apply_Q_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                        A, TW, W, Y );

  // Create a temporary hierarchical view of only the top n-by-n part of A in
  // case m > n so that AT captures the upper triangular factor R. We do the
  // same for Y to ensure conformality.
  FLASH_Part_create_2x1( A,   &AT,    
                              &AB,    FLASH_Obj_scalar_width( A ), FLA_TOP );
  FLASH_Part_create_2x1( Y,   &YT,    
                              &YB,    FLASH_Obj_scalar_width( A ), FLA_TOP );

  FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, AT, YT );

  FLASH_Copy( YT, X );

  // Free the temporary hierarchical views.
  FLASH_Part_free_2x1( &AT,
                       &AB );
  FLASH_Part_free_2x1( &YT,
                       &YB );

  FLASH_Obj_free( &Y );
  FLASH_Obj_free( &W );

  return FLA_SUCCESS;
}

