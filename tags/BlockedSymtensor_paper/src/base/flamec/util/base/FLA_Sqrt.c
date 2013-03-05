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

FLA_Error FLA_Sqrt( FLA_Obj alpha )
{
  FLA_Datatype datatype;
  int          r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Sqrt_check( alpha );

  datatype = FLA_Obj_datatype( alpha );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    if ( *buff_alpha < 0.0F )
      r_val = FLA_FAILURE;
    else
      *buff_alpha = ( float ) sqrt( *buff_alpha );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    if ( *buff_alpha < 0.0 )
      r_val = FLA_FAILURE;
    else
      *buff_alpha = ( double ) sqrt( *buff_alpha );
    
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    if ( buff_alpha->real < 0.0F )
      r_val = FLA_FAILURE;
    else
      buff_alpha->real = ( float ) sqrt( buff_alpha->real );
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    if ( buff_alpha->real < 0.0 )
      r_val = FLA_FAILURE;
    else
      buff_alpha->real = ( double ) sqrt( buff_alpha->real );
    
    break;
  }

  }

  return r_val;
}

