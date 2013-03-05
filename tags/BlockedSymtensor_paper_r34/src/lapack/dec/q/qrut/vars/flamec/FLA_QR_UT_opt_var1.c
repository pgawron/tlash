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

FLA_Error FLA_QR_UT_opt_var1( FLA_Obj A, FLA_Obj t )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_t;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_t    = FLA_Obj_vector_inc( t );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_t = FLA_FLOAT_PTR( t );

      FLA_QR_UT_ops_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_t = FLA_DOUBLE_PTR( t );

      FLA_QR_UT_opd_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_t = FLA_COMPLEX_PTR( t );

      FLA_QR_UT_opc_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_t = FLA_DOUBLE_COMPLEX_PTR( t );

      FLA_QR_UT_opz_var1( m_A,
                          n_A,
                          buff_A, rs_A, cs_A,
                          buff_t, inc_t );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_ops_var1( int m_A,
                              int n_A,
                              float* buff_A, int rs_A, int cs_A, 
                              float* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float* tau1     = buff_t + (i  )*inc_t;

    int    m_ahead  = m_A - i - 1;
    int    n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau1 );
    FLA_Househ2_UT_l_ops( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    //                                       A22 );
    FLA_Apply_H2_UT_l_ops_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opd_var1( int m_A,
                              int n_A,
                              double* buff_A, int rs_A, int cs_A, 
                              double* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double* tau1     = buff_t + (i  )*inc_t;

    int     m_ahead  = m_A - i - 1;
    int     n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau1 );
    FLA_Househ2_UT_l_opd( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    //                                       A22 );
    FLA_Apply_H2_UT_l_opd_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opc_var1( int m_A,
                              int n_A,
                              scomplex* buff_A, int rs_A, int cs_A, 
                              scomplex* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* tau1     = buff_t + (i  )*inc_t;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau1 );
    FLA_Househ2_UT_l_opc( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    //                                       A22 );
    FLA_Apply_H2_UT_l_opc_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_QR_UT_opz_var1( int m_A,
                              int n_A,
                              dcomplex* buff_A, int rs_A, int cs_A, 
                              dcomplex* buff_t, int inc_t )
{
  int min_m_n = min( m_A, n_A );
  int i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* tau1     = buff_t + (i  )*inc_t;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau1 );
    FLA_Househ2_UT_l_opz( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau1 );

    // FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    //                                       A22 );
    FLA_Apply_H2_UT_l_opz_var1( m_ahead,
                                n_ahead,
                                tau1,
                                a21, rs_A,
                                a12t, cs_A,
                                A22, rs_A, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

