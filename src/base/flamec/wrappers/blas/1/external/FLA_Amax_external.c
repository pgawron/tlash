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

FLA_Error FLA_Amax_external( FLA_Obj x, FLA_Obj index )
{
  FLA_Datatype datatype;
  int          num_elem;
  int          inc_x;
  int         *buff_index;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Amax_check( x, index );

  buff_index = ( int * ) FLA_INT_PTR( index );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    *buff_index = 0;
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  num_elem = FLA_Obj_vector_dim( x );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bli_samax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }
  
  case FLA_DOUBLE:
  {
    double* buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bli_damax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }
  
  case FLA_COMPLEX:
  {
    scomplex* buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bli_camax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bli_zamax( num_elem,
               buff_x, inc_x,
               buff_index );

    break;
  }

  }

  return FLA_SUCCESS;
}

