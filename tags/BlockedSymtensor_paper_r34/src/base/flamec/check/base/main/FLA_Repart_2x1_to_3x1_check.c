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

FLA_Error FLA_Repart_2x1_to_3x1_check( FLA_Obj AT,  FLA_Obj *A0,
                                                    FLA_Obj *A1,
                                       FLA_Obj AB,  FLA_Obj *A2,
                                       dim_t   mb,  FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A0 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_topbottom_side( side );
  FLA_Check_error_code( e_val );

  if      ( side == FLA_TOP )
  {
    e_val = FLA_Check_attempted_repart_2x1( AT, mb );
    FLA_Check_error_code( e_val );
  }
  else if ( side == FLA_BOTTOM )
  {
    e_val = FLA_Check_attempted_repart_2x1( AB, mb );
    FLA_Check_error_code( e_val );
  }

  // Needed: check for adjacency, similar to those in FLA_Merge_*().

  return FLA_SUCCESS;
}

