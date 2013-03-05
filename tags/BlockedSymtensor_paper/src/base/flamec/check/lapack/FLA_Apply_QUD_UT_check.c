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

FLA_Error FLA_Apply_QUD_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_storev( storev );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, W );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, U );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, V );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, D );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( U ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_width_equals( W, FLA_Obj_width( R ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, U, R, C );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, V, R, D );
    FLA_Check_error_code( e_val );
  }
  else
  {
  }

  return FLA_SUCCESS;
}

