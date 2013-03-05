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

FLA_Error FLA_Bidiag_UT_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_extract_diagonals_check( A, d, e );

  if ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) )
    r_val = FLA_Bidiag_UT_u_extract_diagonals( A, d, e );
  else
    r_val = FLA_Bidiag_UT_l_extract_diagonals( A, d, e );

  return r_val;
}



FLA_Error FLA_Bidiag_UT_u_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          n_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  inc_e    = FLA_Obj_vector_inc( e );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A = FLA_FLOAT_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = FLA_FLOAT_PTR( e );

      for ( i = 0; i < n_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = *a12t_l;
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A = FLA_DOUBLE_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = FLA_DOUBLE_PTR( e );

      for ( i = 0; i < n_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = *a12t_l;
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = FLA_FLOAT_PTR( e );

      for ( i = 0; i < n_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = a12t_l->real;
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = FLA_DOUBLE_PTR( e );

      for ( i = 0; i < n_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = a12t_l->real;
      }

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_l_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  inc_e    = FLA_Obj_vector_inc( e );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A = FLA_FLOAT_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = FLA_FLOAT_PTR( e );

      for ( i = 0; i < m_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = *a21_t;
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A = FLA_DOUBLE_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = FLA_DOUBLE_PTR( e );

      for ( i = 0; i < m_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = *a21_t;
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = FLA_FLOAT_PTR( e );

      for ( i = 0; i < m_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = a21_t->real;
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = FLA_DOUBLE_PTR( e );

      for ( i = 0; i < m_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = a21_t->real;
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

