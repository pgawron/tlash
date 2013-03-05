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

FLA_Error FLA_Eig_gest_nl_opt_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
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

      FLA_Eig_gest_nl_ops_var1( m_AB,
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

      FLA_Eig_gest_nl_opd_var1( m_AB,
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

      FLA_Eig_gest_nl_opc_var1( m_AB,
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

      FLA_Eig_gest_nl_opz_var1( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_ops_var1( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_0   = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_1h  = FLA_FLOAT_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    y21      = buff_y + (i+1)*inc_y;

    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    float*    b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;
    float*    B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
    //                    FLA_ONE, A22, b21, FLA_ZERO, y21_l );
    bli_shemv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_ahead,
               buff_1,
               A22, rs_A, cs_A,
               b21, rs_B,
               buff_0,
               y21, inc_y );

    // FLA_Scal_external( beta11, a21 );
    bli_sscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_sscals( beta11, alpha11 );
    bli_sscals( beta11, alpha11 );
  
    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_ONE, a21, b21, FLA_ONE, alpha11 );
    bli_sdot2s( BLIS_CONJUGATE,
                m_ahead,
                buff_1,
                a21, rs_A,
                b21, rs_B,
                buff_1,
                alpha11 );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a21 );
    bli_strmv( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJ_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22, rs_B, cs_B,
               a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opd_var1( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_0   = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_1h  = FLA_DOUBLE_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   y21      = buff_y + (i+1)*inc_y;

    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    double*   b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;
    double*   B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
    //                    FLA_ONE, A22, b21, FLA_ZERO, y21_l );
    bli_dhemv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_ahead,
               buff_1,
               A22, rs_A, cs_A,
               b21, rs_B,
               buff_0,
               y21, inc_y );

    // FLA_Scal_external( beta11, a21 );
    bli_dscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_dscals( beta11, alpha11 );
    bli_dscals( beta11, alpha11 );
  
    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_ONE, a21, b21, FLA_ONE, alpha11 );
    bli_ddot2s( BLIS_CONJUGATE,
                m_ahead,
                buff_1,
                a21, rs_A,
                b21, rs_B,
                buff_1,
                alpha11 );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a21 );
    bli_dtrmv( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJ_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22, rs_B, cs_B,
               a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opc_var1( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_0   = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_1h  = FLA_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* y21      = buff_y + (i+1)*inc_y;

    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    scomplex* b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;
    scomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
    //                    FLA_ONE, A22, b21, FLA_ZERO, y21_l );
    bli_chemv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_ahead,
               buff_1,
               A22, rs_A, cs_A,
               b21, rs_B,
               buff_0,
               y21, inc_y );

    // FLA_Scal_external( beta11, a21 );
    bli_cscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_cscals( beta11, alpha11 );
    bli_cscals( beta11, alpha11 );
  
    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_ONE, a21, b21, FLA_ONE, alpha11 );
    bli_cdot2s( BLIS_CONJUGATE,
                m_ahead,
                buff_1,
                a21, rs_A,
                b21, rs_B,
                buff_1,
                alpha11 );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a21 );
    bli_ctrmv( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJ_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22, rs_B, cs_B,
               a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opz_var1( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_0   = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_1h  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* y21      = buff_y + (i+1)*inc_y;

    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    dcomplex* b21      = buff_B + (i  )*cs_B + (i+1)*rs_B;
    dcomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Hemv_external( FLA_LOWER_TRIANGULAR,
    //                    FLA_ONE, A22, b21, FLA_ZERO, y21_l );
    bli_zhemv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_CONJUGATE,
               m_ahead,
               buff_1,
               A22, rs_A, cs_A,
               b21, rs_B,
               buff_0,
               y21, inc_y );

    // FLA_Scal_external( beta11, a21 );
    bli_zscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_zscals( beta11, alpha11 );
    bli_zscals( beta11, alpha11 );
  
    // FLA_Dot2cs_external( FLA_CONJUGATE, FLA_ONE, a21, b21, FLA_ONE, alpha11 );
    bli_zdot2s( BLIS_CONJUGATE,
                m_ahead,
                buff_1,
                a21, rs_A,
                b21, rs_B,
                buff_1,
                alpha11 );

    // FLA_Axpy_external( FLA_ONE_HALF, y21_l, a21 );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                buff_1h,
                y21, inc_y,
                a21, rs_A );

    // FLA_Trmv_external( FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a21 );
    bli_ztrmv( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJ_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22, rs_B, cs_B,
               a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

