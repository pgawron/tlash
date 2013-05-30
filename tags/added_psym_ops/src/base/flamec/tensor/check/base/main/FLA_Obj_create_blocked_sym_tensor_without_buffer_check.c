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

FLA_Error FLA_Obj_create_blocked_sym_tensor_without_buffer_check( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t b, FLA_Obj *obj )
{
  dim_t i;
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  for(i = 0; i < order; i++){
	if(size[i] != size[0])
		printf("FLASHY ERROR: sym tensor mode sizes must be equal\n");
  }
  if((size[0] % b) != 0)
	printf("FLASHY ERROR: sym tensor blk size must equally divide size\n");
  return FLA_SUCCESS;
}

