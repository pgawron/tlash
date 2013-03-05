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

FLA_Error FLA_Lyap_h_opt_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj C )
{
  FLA_Datatype datatype;
  int          m_AC;
  int          rs_A, cs_A;
  int          rs_W, cs_W;
  int          rs_C, cs_C;
  FLA_Obj      W;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &W );

  datatype = FLA_Obj_datatype( A );

  m_AC     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_W     = FLA_Obj_row_stride( W );
  cs_W     = FLA_Obj_col_stride( W );

  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );
 
  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A   = FLA_FLOAT_PTR( A );
      float* buff_W   = FLA_FLOAT_PTR( W );
      float* buff_C   = FLA_FLOAT_PTR( C );
      float* buff_sgn = FLA_FLOAT_PTR( isgn );

      FLA_Lyap_h_ops_var3( m_AC,
                           buff_sgn,
                           buff_A, rs_A, cs_A,
                           buff_W, rs_W, cs_W,
                           buff_C, rs_C, cs_C );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A   = FLA_DOUBLE_PTR( A );
      double* buff_W   = FLA_DOUBLE_PTR( W );
      double* buff_C   = FLA_DOUBLE_PTR( C );
      double* buff_sgn = FLA_DOUBLE_PTR( isgn );

      FLA_Lyap_h_opd_var3( m_AC,
                           buff_sgn,
                           buff_A, rs_A, cs_A,
                           buff_W, rs_W, cs_W,
                           buff_C, rs_C, cs_C );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A   = FLA_COMPLEX_PTR( A );
      scomplex* buff_W   = FLA_COMPLEX_PTR( W );
      scomplex* buff_C   = FLA_COMPLEX_PTR( C );
      scomplex* buff_sgn = FLA_COMPLEX_PTR( isgn );

      FLA_Lyap_h_opc_var3( m_AC,
                           buff_sgn,
                           buff_A, rs_A, cs_A,
                           buff_W, rs_W, cs_W,
                           buff_C, rs_C, cs_C );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A   = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_W   = FLA_DOUBLE_COMPLEX_PTR( W );
      dcomplex* buff_C   = FLA_DOUBLE_COMPLEX_PTR( C );
      dcomplex* buff_sgn = FLA_DOUBLE_COMPLEX_PTR( isgn );

      FLA_Lyap_h_opz_var3( m_AC,
                           buff_sgn,
                           buff_A, rs_A, cs_A,
                           buff_W, rs_W, cs_W,
                           buff_C, rs_C, cs_C );

      break;
    }
  }

  FLA_Obj_free( &W );

  return FLA_SUCCESS;
}



FLA_Error FLA_Lyap_h_ops_var3( int m_AC,
                               float* buff_sgn,
                               float* buff_A, int rs_A, int cs_A, 
                               float* buff_W, int rs_W, int cs_W, 
                               float* buff_C, int rs_C, int cs_C )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  bli_sscalm( BLIS_NO_CONJUGATE,
              m_AC,
              m_AC,
              buff_sgn,
              buff_C, rs_C, cs_C );

  for ( i = 0; i < m_AC; ++i )
  {
    float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
	float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
	float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    c01      = buff_C + (i  )*cs_C + (0  )*rs_C;
    float*    gamma11  = buff_C + (i  )*cs_C + (i  )*rs_C;
	float*    C02      = buff_C + (i+1)*cs_C + (0  )*rs_C;
    float*    c12t     = buff_C + (i+1)*cs_C + (i  )*rs_C;

	float*    W22      = buff_W + (i+1)*cs_W + (i+1)*rs_W;

    float     omega;

    int       m_behind = i;
    int       m_ahead  = m_AC - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Dot2cs( FLA_CONJUGATE, FLA_MINUS_ONE, a01, c01, FLA_ONE, gamma11 );
    bli_sdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                c01, rs_C,
                buff_1,
                gamma11 );

    // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, omega );
    // FLA_Mult_add( FLA_ONE, alpha11, omega );
    // FLA_Inv_scal( omega, gamma11 );
    bli_scopyconj( alpha11, &omega );
    bli_sadd3( alpha11, &omega, &omega );
    bli_sinvscals( &omega, gamma11 );

    // FLA_Axpys( FLA_MINUS_ONE, gamma11, a12t, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, c01, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, C02, a01, FLA_ONE, c12t );
    bli_saxpysv( m_ahead,
                 buff_m1,
                 gamma11,
                 a12t, cs_A,
                 buff_1,
                 c12t, cs_C );

    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               c01,  rs_C,
               buff_1,
               c12t, cs_C );

    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               C02,  rs_C, cs_C,
               a01,  rs_A,
               buff_1,
               c12t, cs_C );

    // FLA_Copyrt( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, A22, W22 );
    // FLA_Shift_diag( FLA_CONJUGATE, alpha11, W22 );
    // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, W22, c12t );
    bli_scopymrt( BLIS_UPPER_TRIANGULAR,
                  BLIS_NO_TRANSPOSE,
                  m_ahead,
                  m_ahead,
                  A22, rs_A, cs_A,
                  W22, rs_W, cs_W );

    bli_sshiftdiag( BLIS_CONJUGATE,
                    0,
                    m_ahead,
                    m_ahead,
                    alpha11,
                    W22, rs_W, cs_W );

    bli_strsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               W22,  rs_W, cs_W,
               c12t, cs_C );

    /*------------------------------------------------------------*/
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Lyap_h_opd_var3( int m_AC,
                               double* buff_sgn,
                               double* buff_A, int rs_A, int cs_A, 
                               double* buff_W, int rs_W, int cs_W, 
                               double* buff_C, int rs_C, int cs_C )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  bli_dscalm( BLIS_NO_CONJUGATE,
              m_AC,
              m_AC,
              buff_sgn,
              buff_C, rs_C, cs_C );

  for ( i = 0; i < m_AC; ++i )
  {
    double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
	double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
	double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   c01      = buff_C + (i  )*cs_C + (0  )*rs_C;
    double*   gamma11  = buff_C + (i  )*cs_C + (i  )*rs_C;
	double*   C02      = buff_C + (i+1)*cs_C + (0  )*rs_C;
    double*   c12t     = buff_C + (i+1)*cs_C + (i  )*rs_C;

	double*   W22      = buff_W + (i+1)*cs_W + (i+1)*rs_W;

    double    omega;

    int       m_behind = i;
    int       m_ahead  = m_AC - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Dot2cs( FLA_CONJUGATE, FLA_MINUS_ONE, a01, c01, FLA_ONE, gamma11 );
    bli_ddot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                c01, rs_C,
                buff_1,
                gamma11 );

    // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, omega );
    // FLA_Mult_add( FLA_ONE, alpha11, omega );
    // FLA_Inv_scal( omega, gamma11 );
    bli_dcopyconj( alpha11, &omega );
    bli_dadd3( alpha11, &omega, &omega );
    bli_dinvscals( &omega, gamma11 );

    // FLA_Axpys( FLA_MINUS_ONE, gamma11, a12t, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, c01, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, C02, a01, FLA_ONE, c12t );
    bli_daxpysv( m_ahead,
                 buff_m1,
                 gamma11,
                 a12t, cs_A,
                 buff_1,
                 c12t, cs_C );

    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               c01,  rs_C,
               buff_1,
               c12t, cs_C );

    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               C02,  rs_C, cs_C,
               a01,  rs_A,
               buff_1,
               c12t, cs_C );

    // FLA_Copyrt( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, A22, W22 );
    // FLA_Shift_diag( FLA_CONJUGATE, alpha11, W22 );
    // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, W22, c12t );
    bli_dcopymrt( BLIS_UPPER_TRIANGULAR,
                  BLIS_NO_TRANSPOSE,
                  m_ahead,
                  m_ahead,
                  A22, rs_A, cs_A,
                  W22, rs_W, cs_W );

    bli_dshiftdiag( BLIS_CONJUGATE,
                    0,
                    m_ahead,
                    m_ahead,
                    alpha11,
                    W22, rs_W, cs_W );

    bli_dtrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               W22,  rs_W, cs_W,
               c12t, cs_C );

    /*------------------------------------------------------------*/
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Lyap_h_opc_var3( int m_AC,
                               scomplex* buff_sgn,
                               scomplex* buff_A, int rs_A, int cs_A, 
                               scomplex* buff_W, int rs_W, int cs_W, 
                               scomplex* buff_C, int rs_C, int cs_C )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  bli_cscalm( BLIS_NO_CONJUGATE,
              m_AC,
              m_AC,
              buff_sgn,
              buff_C, rs_C, cs_C );

  for ( i = 0; i < m_AC; ++i )
  {
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
	scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
	scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* c01      = buff_C + (i  )*cs_C + (0  )*rs_C;
    scomplex* gamma11  = buff_C + (i  )*cs_C + (i  )*rs_C;
	scomplex* C02      = buff_C + (i+1)*cs_C + (0  )*rs_C;
    scomplex* c12t     = buff_C + (i+1)*cs_C + (i  )*rs_C;

	scomplex* W22      = buff_W + (i+1)*cs_W + (i+1)*rs_W;

    scomplex  omega;

    int       m_behind = i;
    int       m_ahead  = m_AC - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Dot2cs( FLA_CONJUGATE, FLA_MINUS_ONE, a01, c01, FLA_ONE, gamma11 );
    bli_cdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                c01, rs_C,
                buff_1,
                gamma11 );

    // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, omega );
    // FLA_Mult_add( FLA_ONE, alpha11, omega );
    // FLA_Inv_scal( omega, gamma11 );
    bli_ccopyconj( alpha11, &omega );
    bli_cadd3( alpha11, &omega, &omega );
    bli_cinvscals( &omega, gamma11 );

    // FLA_Axpys( FLA_MINUS_ONE, gamma11, a12t, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, c01, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, C02, a01, FLA_ONE, c12t );
    bli_caxpysv( m_ahead,
                 buff_m1,
                 gamma11,
                 a12t, cs_A,
                 buff_1,
                 c12t, cs_C );

    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               c01,  rs_C,
               buff_1,
               c12t, cs_C );

    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               C02,  rs_C, cs_C,
               a01,  rs_A,
               buff_1,
               c12t, cs_C );

    // FLA_Copyrt( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, A22, W22 );
    // FLA_Shift_diag( FLA_CONJUGATE, alpha11, W22 );
    // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, W22, c12t );
    bli_ccopymrt( BLIS_UPPER_TRIANGULAR,
                  BLIS_NO_TRANSPOSE,
                  m_ahead,
                  m_ahead,
                  A22, rs_A, cs_A,
                  W22, rs_W, cs_W );

    bli_cshiftdiag( BLIS_CONJUGATE,
                    0,
                    m_ahead,
                    m_ahead,
                    alpha11,
                    W22, rs_W, cs_W );

    bli_ctrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               W22,  rs_W, cs_W,
               c12t, cs_C );

    /*------------------------------------------------------------*/
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Lyap_h_opz_var3( int m_AC,
                               dcomplex* buff_sgn,
                               dcomplex* buff_A, int rs_A, int cs_A, 
                               dcomplex* buff_W, int rs_W, int cs_W, 
                               dcomplex* buff_C, int rs_C, int cs_C )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       i;

  bli_zscalm( BLIS_NO_CONJUGATE,
              m_AC,
              m_AC,
              buff_sgn,
              buff_C, rs_C, cs_C );

  for ( i = 0; i < m_AC; ++i )
  {
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
	dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
	dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* c01      = buff_C + (i  )*cs_C + (0  )*rs_C;
    dcomplex* gamma11  = buff_C + (i  )*cs_C + (i  )*rs_C;
	dcomplex* C02      = buff_C + (i+1)*cs_C + (0  )*rs_C;
    dcomplex* c12t     = buff_C + (i+1)*cs_C + (i  )*rs_C;

	dcomplex* W22      = buff_W + (i+1)*cs_W + (i+1)*rs_W;

    dcomplex  omega;

    int       m_behind = i;
    int       m_ahead  = m_AC - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Dot2cs( FLA_CONJUGATE, FLA_MINUS_ONE, a01, c01, FLA_ONE, gamma11 );
    bli_zdot2s( BLIS_CONJUGATE,
                m_behind,
                buff_m1,
                a01, rs_A,
                c01, rs_C,
                buff_1,
                gamma11 );

    // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, omega );
    // FLA_Mult_add( FLA_ONE, alpha11, omega );
    // FLA_Inv_scal( omega, gamma11 );
    bli_zcopyconj( alpha11, &omega );
    bli_zadd3( alpha11, &omega, &omega );
    bli_zinvscals( &omega, gamma11 );

    // FLA_Axpys( FLA_MINUS_ONE, gamma11, a12t, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, c01, FLA_ONE, c12t );
    // FLA_Gemvc( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, C02, a01, FLA_ONE, c12t );
    bli_zaxpysv( m_ahead,
                 buff_m1,
                 gamma11,
                 a12t, cs_A,
                 buff_1,
                 c12t, cs_C );

    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               c01,  rs_C,
               buff_1,
               c12t, cs_C );

    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_CONJUGATE,
               m_behind,
               m_ahead,
               buff_m1,
               C02,  rs_C, cs_C,
               a01,  rs_A,
               buff_1,
               c12t, cs_C );

    // FLA_Copyrt( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, A22, W22 );
    // FLA_Shift_diag( FLA_CONJUGATE, alpha11, W22 );
    // FLA_Trsv( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, W22, c12t );
    bli_zcopymrt( BLIS_UPPER_TRIANGULAR,
                  BLIS_NO_TRANSPOSE,
                  m_ahead,
                  m_ahead,
                  A22, rs_A, cs_A,
                  W22, rs_W, cs_W );

    bli_zshiftdiag( BLIS_CONJUGATE,
                    0,
                    m_ahead,
                    m_ahead,
                    alpha11,
                    W22, rs_W, cs_W );

    bli_ztrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               W22,  rs_W, cs_W,
               c12t, cs_C );

    /*------------------------------------------------------------*/
  }

  return FLA_SUCCESS;
}

