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

FLA_Error FLA_Eig_gest_iu_opt_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_AB;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          rs_Y, cs_Y;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Y = FLA_FLOAT_PTR( Y );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_iu_ops_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Y = FLA_DOUBLE_PTR( Y );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_iu_opd_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Y = FLA_COMPLEX_PTR( Y );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opc_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Y = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_iu_opz_var3( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_Y, rs_Y, cs_Y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_ops_var3( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_Y, int rs_Y, int cs_Y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float*    buff_m1h = FLA_FLOAT_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float*    y01      = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    float*    Y02      = buff_Y + (i+1)*cs_Y + (0  )*rs_Y;
    float*    y12t     = buff_Y + (i+1)*cs_Y + (i  )*rs_Y;

    float*    b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    float*    B02      = buff_B + (i+1)*cs_B + (0  )*rs_B;
    float*    b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
                a01, rs_A );
  
    // FLA_Inv_scal_external( beta11, a01 );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    // FLA_Ger_external( FLA_ONE, a01, b12t, Y02 );
    bli_sger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              a01,  rs_A,
              b12t, cs_B,
              Y02,  rs_Y, cs_Y );

    // FLA_Copy_external( b12t, y12t );
    // FLA_Scal_external( alpha11, y12t );
    bli_scopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                b12t, cs_B,
                y12t, cs_Y );
    bli_sscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y12t, cs_Y );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B02, a01, FLA_ONE, y12t );
    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_1,
               B02, rs_B, cs_B,
               a01, rs_A,
               buff_1,
               y12t, cs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opd_var3( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_Y, int rs_Y, int cs_Y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  double*   buff_m1h = FLA_DOUBLE_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double*   y01      = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    double*   Y02      = buff_Y + (i+1)*cs_Y + (0  )*rs_Y;
    double*   y12t     = buff_Y + (i+1)*cs_Y + (i  )*rs_Y;

    double*   b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    double*   B02      = buff_B + (i+1)*cs_B + (0  )*rs_B;
    double*   b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
                a01, rs_A );
  
    // FLA_Inv_scal_external( beta11, a01 );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    // FLA_Ger_external( FLA_ONE, a01, b12t, Y02 );
    bli_dger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              a01,  rs_A,
              b12t, cs_B,
              Y02,  rs_Y, cs_Y );

    // FLA_Copy_external( b12t, y12t );
    // FLA_Scal_external( alpha11, y12t );
    bli_dcopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                b12t, cs_B,
                y12t, cs_Y );
    bli_dscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y12t, cs_Y );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B02, a01, FLA_ONE, y12t );
    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_1,
               B02, rs_B, cs_B,
               a01, rs_A,
               buff_1,
               y12t, cs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opc_var3( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_Y, int rs_Y, int cs_Y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex* buff_m1h = FLA_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    scomplex* y01      = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    scomplex* Y02      = buff_Y + (i+1)*cs_Y + (0  )*rs_Y;
    scomplex* y12t     = buff_Y + (i+1)*cs_Y + (i  )*rs_Y;

    scomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    scomplex* B02      = buff_B + (i+1)*cs_B + (0  )*rs_B;
    scomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
                a01, rs_A );
  
    // FLA_Inv_scal_external( beta11, a01 );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    // FLA_Ger_external( FLA_ONE, a01, b12t, Y02 );
    bli_cger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              a01,  rs_A,
              b12t, cs_B,
              Y02,  rs_Y, cs_Y );

    // FLA_Copy_external( b12t, y12t );
    // FLA_Scal_external( alpha11, y12t );
    bli_ccopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                b12t, cs_B,
                y12t, cs_Y );
    bli_cscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y12t, cs_Y );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B02, a01, FLA_ONE, y12t );
    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_1,
               B02, rs_B, cs_B,
               a01, rs_A,
               buff_1,
               y12t, cs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opz_var3( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  dcomplex* buff_m1h = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    dcomplex* y01      = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    dcomplex* Y02      = buff_Y + (i+1)*cs_Y + (0  )*rs_Y;
    dcomplex* y12t     = buff_Y + (i+1)*cs_Y + (i  )*rs_Y;

    dcomplex* b01      = buff_B + (i  )*cs_B + (0  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    dcomplex* B02      = buff_B + (i+1)*cs_B + (0  )*rs_B;
    dcomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
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

    // FLA_Axpy_external( FLA_MINUS_ONE_HALF, y01, a01 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                buff_m1h,
                y01, rs_Y,
                a01, rs_A );
  
    // FLA_Inv_scal_external( beta11, a01 );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   m_behind,
                   beta11,
                   a01, rs_A );

    // FLA_Ger_external( FLA_ONE, a01, b12t, Y02 );
    bli_zger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_behind,
              m_ahead,
              buff_1,
              a01,  rs_A,
              b12t, cs_B,
              Y02,  rs_Y, cs_Y );

    // FLA_Copy_external( b12t, y12t );
    // FLA_Scal_external( alpha11, y12t );
    bli_zcopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                b12t, cs_B,
                y12t, cs_Y );
    bli_zscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                alpha11,
                y12t, cs_Y );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE,
    //                     FLA_ONE, B02, a01, FLA_ONE, y12t );
    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_1,
               B02, rs_B, cs_B,
               a01, rs_A,
               buff_1,
               y12t, cs_Y );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

