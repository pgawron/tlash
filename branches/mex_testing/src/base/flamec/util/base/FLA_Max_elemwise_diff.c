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

double FLA_Max_elemwise_diff( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype;
  dim_t        i, j;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;
  dim_t        rs_B, cs_B;
  double       diff;
  double       d_max = 0.0;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Max_elemwise_diff_check( A, B );

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_a = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_b = ( float * ) FLA_FLOAT_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ] - buff_b[ j*cs_B + i*rs_B ] );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_a = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_b = ( double * ) FLA_DOUBLE_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ] - buff_b[ j*cs_B + i*rs_B ] );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_a = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_b = ( scomplex * ) FLA_COMPLEX_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].real - buff_b[ j*cs_B + i*rs_B ].real );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);

        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].imag - buff_b[ j*cs_B + i*rs_B ].imag );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_a = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_b = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].real - buff_b[ j*cs_B + i*rs_B ].real );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);

        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].imag - buff_b[ j*cs_B + i*rs_B ].imag );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  }

  
  return d_max;
}

