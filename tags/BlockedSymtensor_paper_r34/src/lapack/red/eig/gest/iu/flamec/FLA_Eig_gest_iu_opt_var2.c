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

FLA_Error FLA_Eig_gest_iu_opt_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_AB;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          inc_y;
  FLA_Obj      yL, yR;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  FLA_Part_1x2( Y,    &yL, &yR,     1, FLA_LEFT );

  inc_y    = FLA_Obj_vector_inc( yL );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_y = FLA_FLOAT_PTR( yL );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_iu_ops_var2( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_y = FLA_DOUBLE_PTR( yL );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_iu_opd_var2( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_y = FLA_COMPLEX_PTR( yL );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opc_var2( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( yL );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opz_var2( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_ops_var2( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_1h  = FLA_FLOAT_PTR( FLA_ONE_HALF );
  float*    buff_0   = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float*    y01      = buff_y + (0  )*inc_y;

    float*    b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE_HALF, A00, b01, FLA_ZERO, y01_l );
    bli_shemv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_behind,
               buff_1h,
               A00, rs_A, cs_A,
               b01, rs_B,
               buff_0,
               y01, inc_y );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );
    bli_sdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                b01, rs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_sinvscals( beta11, alpha11 );
    bli_sinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );
    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               b01, rs_B,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Inv_scal_external( beta11, a01 );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opd_var2( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_1h  = FLA_DOUBLE_PTR( FLA_ONE_HALF );
  double*   buff_0   = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double*   y01      = buff_y + (0  )*inc_y;

    double*   b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE_HALF, A00, b01, FLA_ZERO, y01_l );
    bli_dhemv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_behind,
               buff_1h,
               A00, rs_A, cs_A,
               b01, rs_B,
               buff_0,
               y01, inc_y );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );
    bli_ddot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                b01, rs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_dinvscals( beta11, alpha11 );
    bli_dinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );
    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               b01, rs_B,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Inv_scal_external( beta11, a01 );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opc_var2( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_1h  = FLA_COMPLEX_PTR( FLA_ONE_HALF );
  scomplex* buff_0   = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    scomplex* y01      = buff_y + (0  )*inc_y;

    scomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE_HALF, A00, b01, FLA_ZERO, y01_l );
    bli_chemv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_behind,
               buff_1h,
               A00, rs_A, cs_A,
               b01, rs_B,
               buff_0,
               y01, inc_y );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );
    bli_cdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                b01, rs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_cinvscals( beta11, alpha11 );
    bli_cinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );
    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               b01, rs_B,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Inv_scal_external( beta11, a01 );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opz_var2( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_1h  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  dcomplex* buff_0   = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    dcomplex* y01      = buff_y + (0  )*inc_y;

    dcomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Hemvc_external( FLA_UPPER_TRIANGULAR, FLA_NO_CONJUGATE,
    //                     FLA_ONE_HALF, A00, b01, FLA_ZERO, y01_l );
    bli_zhemv( BLIS_UPPER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_behind,
               buff_1h,
               A00, rs_A, cs_A,
               b01, rs_B,
               buff_0,
               y01, inc_y );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, b01, FLA_ONE, alpha11 );
    bli_zdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                b01, rs_B,
                buff_1,
                alpha11 );

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_zinvscals( beta11, alpha11 );
    bli_zinvscals( beta11, alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, A02, b01, FLA_ONE, a12t );
    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               b01, rs_B,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( FLA_MINUS_ONE, y01_l, a01 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1,
                y01, inc_y,
                a01, rs_A );

    // FLA_Inv_scal_external( beta11, a01 );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

