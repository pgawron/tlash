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

FLA_Error FLA_Add_to_diag( void* diag_value, FLA_Obj A )
{
  FLA_Datatype datatype;
  dim_t        j, min_m_n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Add_to_diag_check( diag_value, A );

  datatype = FLA_Obj_datatype( A );
  min_m_n  = FLA_Obj_min_dim( A );
  rs       = FLA_Obj_row_stride( A );
  cs       = FLA_Obj_col_stride( A );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *value_ptr = ( float * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
      buff_A[ j*cs + j*rs ] += *value_ptr;

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *value_ptr = ( double * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
      buff_A[ j*cs + j*rs ] += *value_ptr;

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *value_ptr = ( scomplex * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
    {
      buff_A[ j*cs + j*rs ].real += value_ptr->real;
      buff_A[ j*cs + j*rs ].imag += value_ptr->imag;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *value_ptr = ( dcomplex * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
    {
      buff_A[ j*cs + j*rs ].real += value_ptr->real;
      buff_A[ j*cs + j*rs ].imag += value_ptr->imag;
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

