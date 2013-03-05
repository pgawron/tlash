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

FLA_Error FLA_Invert( FLA_Conj conj, FLA_Obj x )
{
  FLA_Datatype datatype;
  int          n_elem;
  int          inc_x;
  conj_t       blis_conj;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING ) 
    FLA_Invert_check( conj, x );

  if ( FLA_Obj_has_zero_dim( x ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( x );

  n_elem   = FLA_Obj_vector_dim( x );
  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bli_sinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bli_dinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bli_cinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  { 
    dcomplex *buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bli_zinvertv( blis_conj,
                  n_elem,
                  buff_x, inc_x );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

