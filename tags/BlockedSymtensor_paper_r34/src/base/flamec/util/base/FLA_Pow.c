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

FLA_Error FLA_Pow( FLA_Obj base, FLA_Obj exp, FLA_Obj btoe )
{
  FLA_Datatype datatype;
  int          r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Pow_check( base, exp, btoe );

  datatype = FLA_Obj_datatype( base );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_base = ( float * ) FLA_FLOAT_PTR( base );
    float *buff_exp  = ( float * ) FLA_FLOAT_PTR( exp );
    float *buff_btoe = ( float * ) FLA_FLOAT_PTR( btoe );

    *buff_btoe = ( float ) pow( *buff_base, *buff_exp );
    
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_base = ( double * ) FLA_DOUBLE_PTR( base );
    double *buff_exp  = ( double * ) FLA_DOUBLE_PTR( exp );
    double *buff_btoe = ( double * ) FLA_DOUBLE_PTR( btoe );

    *buff_btoe = ( double ) pow( *buff_base, *buff_exp );
    
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_base = ( scomplex * ) FLA_COMPLEX_PTR( base );
    scomplex *buff_exp  = ( scomplex * ) FLA_COMPLEX_PTR( exp );
    scomplex *buff_btoe = ( scomplex * ) FLA_COMPLEX_PTR( btoe );

    buff_btoe->real = ( float ) pow( buff_base->real, buff_exp->real );
    buff_btoe->imag = 0.0;
    
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_base = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( base );
    dcomplex *buff_exp  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( exp );
    dcomplex *buff_btoe = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( btoe );

    buff_btoe->real = ( double ) pow( buff_base->real, buff_exp->real );
    buff_btoe->imag = 0.0;
    
    break;
  }

  }

  return r_val;
}

