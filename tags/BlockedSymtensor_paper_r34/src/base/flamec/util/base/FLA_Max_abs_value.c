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

FLA_Error FLA_Max_abs_value( FLA_Obj A, FLA_Obj maxabs )
{
  FLA_Datatype datatype;
  FLA_Datatype dt_maxabs;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Max_abs_value_check( A, maxabs );

  datatype  = FLA_Obj_datatype( A );
  dt_maxabs = FLA_Obj_datatype( maxabs );

  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );
  rs_A      = FLA_Obj_row_stride( A );
  cs_A      = FLA_Obj_col_stride( A );
 
 
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float*    buff_A      = ( float    * ) FLA_FLOAT_PTR( A );
    float*    buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

    bli_smaxabsm( m_A,
                  n_A,
                  buff_A, rs_A, cs_A,
                  buff_maxabs );

    break;
  }

  case FLA_DOUBLE:
  {
    double*   buff_A      = ( double   * ) FLA_DOUBLE_PTR( A );
    double*   buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

    bli_dmaxabsm( m_A,
                  n_A,
                  buff_A, rs_A, cs_A,
                  buff_maxabs );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( dt_maxabs == FLA_FLOAT )
    {
      scomplex* buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float*    buff_maxabs = ( float    * ) FLA_FLOAT_PTR( maxabs );

      bli_cmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    buff_maxabs );
    }
    else
    {
      scomplex* buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buff_maxabs = ( scomplex * ) FLA_COMPLEX_PTR( maxabs );

      bli_cmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    &(buff_maxabs->real) );

      buff_maxabs->imag = 0.0;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( dt_maxabs == FLA_DOUBLE )
    {
      dcomplex* buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_maxabs = ( double   * ) FLA_DOUBLE_PTR( maxabs );

      bli_zmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    buff_maxabs );
    }
    else
    {
      dcomplex* buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_maxabs = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( maxabs );

      bli_zmaxabsm( m_A,
                    n_A,
                    buff_A, rs_A, cs_A,
                    &(buff_maxabs->real) );

      buff_maxabs->imag = 0.0;
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

