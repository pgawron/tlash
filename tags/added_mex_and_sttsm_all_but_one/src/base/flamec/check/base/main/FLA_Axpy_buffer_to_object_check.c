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

FLA_Error FLA_Axpy_buffer_to_object_check( FLA_Trans trans, FLA_Obj alpha, dim_t m, dim_t n, void* A_buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_real_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( B, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A_buffer );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_dims( trans, m, n, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_strides( m, n, rs, cs );
  FLA_Check_error_code( e_val );

  if ( trans == FLA_NO_TRANSPOSE )
  {
    e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, B );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_submatrix_dims_and_offset( n, m, i, j, B );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

