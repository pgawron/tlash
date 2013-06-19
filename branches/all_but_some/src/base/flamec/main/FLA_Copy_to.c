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

FLA_Error FLA_Copy_buffer_to_object( FLA_Trans trans, dim_t m, dim_t n, void* A_buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj B )
{
  FLA_Obj  A;
  FLA_Obj  BTL, BTR, 
           BBL, Bij;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Copy_buffer_to_object_check( trans, m, n, A_buffer, rs, cs, i, j, B );

  FLA_Part_2x2( B,  &BTL, &BTR,
                    &BBL, &Bij,     i, j, FLA_TL );

  FLA_Obj_create_without_buffer( FLA_Obj_datatype( B ), m, n, &A );
  FLA_Obj_attach_buffer( A_buffer, rs, cs, &A );

  FLA_Copyt_external( trans, A, Bij );

  FLA_Obj_free_without_buffer( &A );

  return FLA_SUCCESS;
}



FLA_Error FLA_Copy_object_to_buffer( FLA_Trans trans, dim_t i, dim_t j, FLA_Obj A, dim_t m, dim_t n, void* B_buffer, dim_t rs, dim_t cs )
{
  FLA_Obj  B;
  FLA_Obj  ATL, ATR, 
           ABL, Aij;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Copy_object_to_buffer_check( trans, i, j, A, m, n, B_buffer, rs, cs );

  FLA_Part_2x2( A,  &ATL, &ATR,
                    &ABL, &Aij,     i, j, FLA_TL );

  FLA_Obj_create_without_buffer( FLA_Obj_datatype( A ), m, n, &B );
  FLA_Obj_attach_buffer( B_buffer, rs, cs, &B );

  FLA_Copyt_external( trans, Aij, B );

  FLA_Obj_free_without_buffer( &B );

  return FLA_SUCCESS;
}

