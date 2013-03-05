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

FLA_Error FLA_Tridiag_UT_l_ofu_var3( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;
  FLA_Obj   Z;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Z );

  r_val = FLA_Tridiag_UT_l_step_ofu_var3( A, Z, T );

  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Tridiag_UT_l_step_ofu_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, m_T;
  int          rs_A, cs_A;
  int          rs_Z, cs_Z;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  m_T      = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Z = FLA_FLOAT_PTR( Z );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Tridiag_UT_l_step_ofs_var3( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Z = FLA_DOUBLE_PTR( Z );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Tridiag_UT_l_step_ofd_var3( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Z = FLA_COMPLEX_PTR( Z );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Tridiag_UT_l_step_ofc_var3( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Z = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Tridiag_UT_l_step_ofz_var3( m_A,
                                      m_T,
                                      buff_A, rs_A, cs_A,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_ofs_var3( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_Z, int rs_Z, int cs_Z,
                                          float* buff_T, int rs_T, int cs_T )
{
  float*    buff_2  = FLA_FLOAT_PTR( FLA_TWO );
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     first_elem, last_elem;
  float     beta;
  float     inv_tau11;
  float     minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Set( FLA_ZERO, Z );
  bli_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    float*    Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    float*    z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    float*    ABL      = a10t;
    float*    ZBL      = z10t;

    float*    a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, z10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
    bli_sgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a2,   rs_A );
    bli_sgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_ops( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bli_sdiv3( buff_1, tau11, &inv_tau11 );
      bli_sneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bli_ssymv( BLIS_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f01 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f01, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d01, FLA_ONE, z21 );
      // FLA_Copy( d01, t01 );
      FLA_Fused_UZhu_ZUhu_ops_var1( m_ahead,
                                    n_behind,
                                    buff_m1,
                                    A20, rs_A, cs_A,
                                    Z20, rs_Z, cs_Z,
                                    t01, rs_T,
                                    a21, rs_A,
                                    z21, rs_Z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bli_sdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bli_sinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bli_sscals( &minus_inv_tau11, &beta );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bli_sscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_ofd_var3( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_Z, int rs_Z, int cs_Z,
                                          double* buff_T, int rs_T, int cs_T )
{
  double*   buff_2  = FLA_DOUBLE_PTR( FLA_TWO );
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    first_elem, last_elem;
  double    beta;
  double    inv_tau11;
  double    minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Set( FLA_ZERO, Z );
  bli_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    double*   Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    double*   z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    double*   ABL      = a10t;
    double*   ZBL      = z10t;

    double*   a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, z10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
    bli_dgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a2,   rs_A );
    bli_dgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opd( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bli_ddiv3( buff_1, tau11, &inv_tau11 );
      bli_dneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bli_dsymv( BLIS_LOWER_TRIANGULAR,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f01 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f01, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d01, FLA_ONE, z21 );
      // FLA_Copy( d01, t01 );
      FLA_Fused_UZhu_ZUhu_opd_var1( m_ahead,
                                    n_behind,
                                    buff_m1,
                                    A20, rs_A, cs_A,
                                    Z20, rs_Z, cs_Z,
                                    t01, rs_T,
                                    a21, rs_A,
                                    z21, rs_Z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bli_ddot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bli_dinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bli_dscals( &minus_inv_tau11, &beta );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bli_dscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_ofc_var3( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_Z, int rs_Z, int cs_Z,
                                          scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_2  = FLA_COMPLEX_PTR( FLA_TWO );
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  first_elem, last_elem;
  scomplex  beta;
  scomplex  inv_tau11;
  scomplex  minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Set( FLA_ZERO, Z );
  bli_csetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    scomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    scomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    scomplex* ABL      = a10t;
    scomplex* ZBL      = z10t;

    scomplex* a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, z10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
    bli_cgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a2,   rs_A );
    bli_cgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opc( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bli_cdiv3( buff_1, tau11, &inv_tau11 );
      bli_cneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bli_chemv( BLIS_LOWER_TRIANGULAR,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f01 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f01, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d01, FLA_ONE, z21 );
      // FLA_Copy( d01, t01 );
      FLA_Fused_UZhu_ZUhu_opc_var1( m_ahead,
                                    n_behind,
                                    buff_m1,
                                    A20, rs_A, cs_A,
                                    Z20, rs_Z, cs_Z,
                                    t01, rs_T,
                                    a21, rs_A,
                                    z21, rs_Z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bli_cdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bli_cinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bli_cscals( &minus_inv_tau11, &beta );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bli_cscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Tridiag_UT_l_step_ofz_var3( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_Z, int rs_Z, int cs_Z,
                                          dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_2  = FLA_DOUBLE_COMPLEX_PTR( FLA_TWO );
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  first_elem, last_elem;
  dcomplex  beta;
  dcomplex  inv_tau11;
  dcomplex  minus_inv_tau11;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Set( FLA_ZERO, Z );
  bli_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* z10t     = buff_Z + (0  )*cs_Z + (i  )*rs_Z;
    dcomplex* Z20      = buff_Z + (0  )*cs_Z + (i+1)*rs_Z;
    dcomplex* z21      = buff_Z + (i  )*cs_Z + (i+1)*rs_Z;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* a10t_r   = a10t   + (i-1)*cs_A + (0  )*rs_A;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    dcomplex* ABL      = a10t;
    dcomplex* ZBL      = z10t;

    dcomplex* a2       = alpha11;

    int       m_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a10t_r, last_elem );
      // FLA_Set( FLA_ONE, a10t_r );
      last_elem = *a10t_r;
      *a10t_r = *buff_1;
    }

    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, z10t, FLA_ONE, a2 );
    // FLA_Gemvc( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a10t, FLA_ONE, a2 );
    bli_zgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a2,   rs_A );
    bli_zgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL,  rs_Z, cs_Z,
               a10t, cs_A,
               buff_1,
               a2,   rs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a10t_r );
      *a10t_r = last_elem;
    }

    if ( m_ahead > 0 )
    {
      // FLA_Househ2_UT( FLA_LEFT,
      //                 a21_t,
      //                 a21_b, tau11 );
      FLA_Househ2_UT_l_opz( m_ahead - 1,
                            a21_t,
                            a21_b, rs_A,
                            tau11 );

      // FLA_Set( FLA_ONE, inv_tau11 );
      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, inv_tau11 );
      // FLA_Copy( inv_tau11, minus_inv_tau11 );
      // FLA_Scal( FLA_MINUS_ONE, minus_inv_tau11 );
      bli_zdiv3( buff_1, tau11, &inv_tau11 );
      bli_zneg2( &inv_tau11, &minus_inv_tau11 );

      // FLA_Copy( a21_t, first_elem );
      // FLA_Set( FLA_ONE, a21_t );
      first_elem = *a21_t;
      *a21_t = *buff_1;

      // FLA_Hemv( FLA_LOWER_TRIANGULAR, FLA_ONE, A22, a21, FLA_ZERO, z21 );
      bli_zhemv( BLIS_LOWER_TRIANGULAR,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 z21, rs_Z );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, d01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, a21, FLA_ZERO, f01 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f01, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, d01, FLA_ONE, z21 );
      // FLA_Copy( d01, t01 );
      FLA_Fused_UZhu_ZUhu_opz_var1( m_ahead,
                                    n_behind,
                                    buff_m1,
                                    A20, rs_A, cs_A,
                                    Z20, rs_Z, cs_Z,
                                    t01, rs_T,
                                    a21, rs_A,
                                    z21, rs_Z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z21, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      bli_zdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z21, rs_Z,
                &beta );
      bli_zinvscals( buff_2, &beta );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z21 );
      // FLA_Scal( inv_tau11, z21 );
      bli_zscals( &minus_inv_tau11, &beta );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z21, rs_Z );
      bli_zscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z21, rs_Z );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

