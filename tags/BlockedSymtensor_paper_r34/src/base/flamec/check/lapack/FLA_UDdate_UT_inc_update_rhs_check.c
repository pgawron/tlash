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

FLA_Error FLA_UDdate_UT_inc_update_rhs_check( FLA_Obj T, FLA_Obj bR, FLA_Obj C, FLA_Obj bC, FLA_Obj D, FLA_Obj bD )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bC );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bD );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( T, max( FLA_Obj_length( C ),
                                                  FLA_Obj_length( D ) ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( C, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( D, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, C, bR, bC );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, bR, bD );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}
