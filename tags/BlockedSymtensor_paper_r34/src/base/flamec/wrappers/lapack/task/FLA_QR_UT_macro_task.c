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

extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_macro_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
   FLA_Error r_val;
   FLA_Obj   A_flat;

   if ( FLA_Obj_length( A ) > 1 )
   {
      FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
  
      r_val = FLA_QR_UT_internal( A_flat, T, fla_qrut_cntl_leaf );
  
      FLASH_Copy_flat_to_hier( A_flat, 0, 0, A );
  
      FLA_Obj_free( &A_flat );
   }
   else
   {
      r_val = FLA_QR_UT_task( *FLASH_OBJ_PTR_AT( A ), T, cntl );
   }

   return r_val;
}
