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

FLA_Error FLA_Nrm2_external( FLA_Obj x, FLA_Obj norm_x )
{
  FLA_Datatype datatype;
  int          num_elem;
  int          inc_x;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Nrm2_check( x, norm_x );

  if ( FLA_Obj_has_zero_dim( x ) )
  {
    FLA_Set( FLA_ZERO, norm_x );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( x );

  inc_x    = FLA_Obj_vector_inc( x );
  num_elem = FLA_Obj_vector_dim( x );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_x      = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_norm_x = ( float * ) FLA_FLOAT_PTR( norm_x );

    bli_snrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_x      = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_norm_x = ( double * ) FLA_DOUBLE_PTR( norm_x );

    bli_dnrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_x      = ( scomplex * ) FLA_COMPLEX_PTR( x );
    float    *buff_norm_x = ( float    * ) FLA_COMPLEX_PTR( norm_x );

    bli_cnrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_x      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    double   *buff_norm_x = ( double   * ) FLA_DOUBLE_COMPLEX_PTR( norm_x );

    bli_znrm2( num_elem,
               buff_x, inc_x,
               buff_norm_x );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

