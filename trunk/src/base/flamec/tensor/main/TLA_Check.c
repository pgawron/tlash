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

FLA_Error FLA_Check_adjacent_objects_1xmode2( FLA_Obj AT,
                                          FLA_Obj AB, dim_t mode )
{
  dim_t i;
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_dimsize( AT, mode ) != FLA_Obj_dimsize( AB, mode ) )
    e_val = FLA_ADJACENT_OBJECT_DIM_MISMATCH;

  for(i = 0; i < AT.order; i++){
    if(i == mode){
      if ( AB.offset[i] != AT.offset[i] + FLA_Obj_dimsize( AT, mode ) )
        e_val = FLA_OBJECTS_NOT_VERTICALLY_ADJ;
    }
    else{
      if ( AB.offset[i] != AT.offset[i] )
        e_val = FLA_OBJECTS_NOT_VERTICALLY_ALIGNED;
    }
  }
  return e_val;
}

FLA_Error FLA_Check_attempted_repart_1xmode2( FLA_Obj A_side, dim_t mode, dim_t b )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( b > FLA_Obj_dimsize( A_side, mode ) )
    e_val = FLA_ATTEMPTED_OVER_REPART_2X1;

  return e_val;
}
