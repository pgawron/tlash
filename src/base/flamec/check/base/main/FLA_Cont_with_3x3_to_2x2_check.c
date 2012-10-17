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

FLA_Error FLA_Cont_with_3x3_to_2x2_check( FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                       FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                          FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                       FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( ATL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABL );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ATR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( ABR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A00 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A10 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A20 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A01 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A02 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

