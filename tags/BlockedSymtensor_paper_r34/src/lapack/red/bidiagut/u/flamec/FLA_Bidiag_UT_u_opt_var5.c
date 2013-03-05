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

FLA_Error FLA_Bidiag_UT_u_opt_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error    r_val;
  FLA_Obj      Y, Z;
  FLA_Datatype datatype_A;
  dim_t        m_A, n_A;

  datatype_A = FLA_Obj_datatype( A );
  m_A        = FLA_Obj_length( A );
  n_A        = FLA_Obj_width( A );
  
  FLA_Obj_create( datatype_A, n_A, n_A, 0, 0, &Y );
  FLA_Obj_create( datatype_A, m_A, n_A, 0, 0, &Z );

  r_val = FLA_Bidiag_UT_u_step_opt_var5( A, Y, Z, TU, TV );

  FLA_Obj_free( &Y );
  FLA_Obj_free( &Z );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_u_step_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S )
{
  FLA_Datatype datatype;
  int          m_A, n_A, m_TS;
  int          rs_A, cs_A;
  int          rs_Y, cs_Y;
  int          rs_Z, cs_Z;
  int          rs_T, cs_T;
  int          rs_S, cs_S;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  m_TS     = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  
  rs_S     = FLA_Obj_row_stride( S );
  cs_S     = FLA_Obj_col_stride( S );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_Y = FLA_FLOAT_PTR( Y );
      float* buff_Z = FLA_FLOAT_PTR( Z );
      float* buff_T = FLA_FLOAT_PTR( T );
      float* buff_S = FLA_FLOAT_PTR( S );

      FLA_Bidiag_UT_u_step_ops_var5( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_Y = FLA_DOUBLE_PTR( Y );
      double* buff_Z = FLA_DOUBLE_PTR( Z );
      double* buff_T = FLA_DOUBLE_PTR( T );
      double* buff_S = FLA_DOUBLE_PTR( S );

      FLA_Bidiag_UT_u_step_opd_var5( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_Y = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z = FLA_COMPLEX_PTR( Z );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );
      scomplex* buff_S = FLA_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opc_var5( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_Y = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );
      dcomplex* buff_S = FLA_DOUBLE_COMPLEX_PTR( S );

      FLA_Bidiag_UT_u_step_opz_var5( m_A,
                                     n_A,
                                     m_TS,
                                     buff_A, rs_A, cs_A,
                                     buff_Y, rs_Y, cs_Y,
                                     buff_Z, rs_Z, cs_Z,
                                     buff_T, rs_T, cs_T,
                                     buff_S, rs_S, cs_S );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_ops_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_Y, int rs_Y, int cs_Y, 
                                         float* buff_Z, int rs_Z, int cs_Z, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     beta;
  float     last_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  float*    buff_u  = ( float*    ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_v  = ( float*    ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_d  = ( float*    ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_e  = ( float*    ) FLA_malloc( n_A * sizeof( *buff_A ) );
  float*    buff_f  = ( float*    ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_g  = ( float*    ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u   = 1;
  int       inc_v   = 1;
  int       inc_d   = 1;
  int       inc_e   = 1;
  int       inc_f   = 1;
  int       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bli_ssetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bli_ssetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    float*    a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    float*    a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    float*    A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    float*    A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    float*    y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    float*    Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    float*    y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    float*    z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    float*    Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    float*    z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    float*    t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    float*    tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    float*    s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    float*    sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    float*    u21       = buff_u  + (i+1)*inc_u;

    float*    v21       = buff_v  + (i+1)*inc_v;

    float*    d0        = buff_d  + (0  )*inc_d;

    float*    e0        = buff_e  + (0  )*inc_e;

    float*    f0        = buff_f  + (0  )*inc_f;

    float*    g0        = buff_g  + (0  )*inc_g;

    float*    v21_t     = v21     + (0  )*inc_v;
    float*    v21_b     = v21     + (1  )*inc_v;

    float*    a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    float*    a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    float*    a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    float*    ABL       = a10t;
    float*    ZBL       = z10t;

    float*    a2        = alpha11;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bli_sgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bli_sgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bli_sgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bli_sgemv( BLIS_CONJ_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21 );
    FLA_Househ2_UT_l_ops( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bli_scopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                a21, rs_A,
                u21, inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21, FLA_ONE, y21 );
      bli_scopyv( BLIS_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  rs_Y );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21, FLA_ZERO, e0 );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_0,
                 d0,  inc_d );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 u21, inc_u,
                 buff_0,
                 e0,  inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bli_scopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bli_sgemv( BLIS_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bli_sinvscalv( BLIS_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bli_saxpyv( BLIS_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  rs_Y,
                  a12t, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_ops( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bli_scopyv( BLIS_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bli_sdot( BLIS_CONJUGATE,
                n_ahead,
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bli_sscals( buff_m1, &beta );

      // FLA_Copy( u21, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, v21, beta, z21 );
      bli_scopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  u21, inc_u,
                  z21, rs_Z );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, rs_Z );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_ONE, Y20, v21, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, g0 );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 m_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 v21, inc_v,
                 buff_0,
                 f0,  inc_f );
      bli_sgemv( BLIS_CONJ_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 g0,  inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, z21 );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 z21, rs_Z );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bli_sinvscalv( BLIS_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );

      // FLA_Copy( g0, s01 );
      bli_scopyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
    }
	else // if ( n_ahead == 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
      bli_scopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_u );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opd_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_Y, int rs_Y, int cs_Y, 
                                         double* buff_Z, int rs_Z, int cs_Z, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    beta;
  double    last_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  double*   buff_u  = ( double*   ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_v  = ( double*   ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_d  = ( double*   ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_e  = ( double*   ) FLA_malloc( n_A * sizeof( *buff_A ) );
  double*   buff_f  = ( double*   ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_g  = ( double*   ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u   = 1;
  int       inc_v   = 1;
  int       inc_d   = 1;
  int       inc_e   = 1;
  int       inc_f   = 1;
  int       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bli_dsetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bli_dsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    double*   a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    double*   a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    double*   A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    double*   A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    double*   y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    double*   Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    double*   y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    double*   z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    double*   Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    double*   z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    double*   t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    double*   tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    double*   s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    double*   sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    double*   u21       = buff_u  + (i+1)*inc_u;

    double*   v21       = buff_v  + (i+1)*inc_v;

    double*   d0        = buff_d  + (0  )*inc_d;

    double*   e0        = buff_e  + (0  )*inc_e;

    double*   f0        = buff_f  + (0  )*inc_f;

    double*   g0        = buff_g  + (0  )*inc_g;

    double*   v21_t     = v21     + (0  )*inc_v;
    double*   v21_b     = v21     + (1  )*inc_v;

    double*   a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    double*   a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    double*   a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    double*   ABL       = a10t;
    double*   ZBL       = z10t;

    double*   a2        = alpha11;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bli_dgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bli_dgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bli_dgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bli_dgemv( BLIS_CONJ_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21 );
    FLA_Househ2_UT_l_opd( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bli_dcopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                a21, rs_A,
                u21, inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21, FLA_ONE, y21 );
      bli_dcopyv( BLIS_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  rs_Y );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21, FLA_ZERO, e0 );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_0,
                 d0,  inc_d );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 u21, inc_u,
                 buff_0,
                 e0,  inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bli_dcopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bli_dgemv( BLIS_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bli_dinvscalv( BLIS_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bli_daxpyv( BLIS_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  rs_Y,
                  a12t, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opd( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bli_ddot( BLIS_CONJUGATE,
                n_ahead,
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bli_dscals( buff_m1, &beta );

      // FLA_Copy( u21, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, v21, beta, z21 );
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  u21, inc_u,
                  z21, rs_Z );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, rs_Z );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_ONE, Y20, v21, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, g0 );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 m_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 v21, inc_v,
                 buff_0,
                 f0,  inc_f );
      bli_dgemv( BLIS_CONJ_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 g0,  inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, z21 );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 z21, rs_Z );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bli_dinvscalv( BLIS_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );

      // FLA_Copy( g0, s01 );
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
    }
	else // if ( n_ahead == 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
      bli_dcopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_u );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opc_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_Y, int rs_Y, int cs_Y, 
                                         scomplex* buff_Z, int rs_Z, int cs_Z, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  beta;
  scomplex  last_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  scomplex* buff_u  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_v  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_d  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_e  = ( scomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  scomplex* buff_f  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_g  = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u   = 1;
  int       inc_v   = 1;
  int       inc_d   = 1;
  int       inc_e   = 1;
  int       inc_f   = 1;
  int       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bli_csetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bli_csetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    scomplex* Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    scomplex* y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    scomplex* z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    scomplex* Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    scomplex* z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    scomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    scomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    scomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    scomplex* u21       = buff_u  + (i+1)*inc_u;

    scomplex* v21       = buff_v  + (i+1)*inc_v;

    scomplex* d0        = buff_d  + (0  )*inc_d;

    scomplex* e0        = buff_e  + (0  )*inc_e;

    scomplex* f0        = buff_f  + (0  )*inc_f;

    scomplex* g0        = buff_g  + (0  )*inc_g;

    scomplex* v21_t     = v21     + (0  )*inc_v;
    scomplex* v21_b     = v21     + (1  )*inc_v;

    scomplex* a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    scomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    scomplex* ABL       = a10t;
    scomplex* ZBL       = z10t;

    scomplex* a2        = alpha11;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bli_cgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bli_cgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bli_cgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bli_cgemv( BLIS_CONJ_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21 );
    FLA_Househ2_UT_l_opc( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bli_ccopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                a21, rs_A,
                u21, inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21, FLA_ONE, y21 );
      bli_ccopyv( BLIS_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  rs_Y );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21, FLA_ZERO, e0 );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_0,
                 d0,  inc_d );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 u21, inc_u,
                 buff_0,
                 e0,  inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bli_ccopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bli_cgemv( BLIS_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bli_cinvscalv( BLIS_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bli_caxpyv( BLIS_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  rs_Y,
                  a12t, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opc( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bli_cdot( BLIS_CONJUGATE,
                n_ahead,
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bli_cscals( buff_m1, &beta );

      // FLA_Copy( u21, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, v21, beta, z21 );
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  u21, inc_u,
                  z21, rs_Z );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, rs_Z );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_ONE, Y20, v21, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, g0 );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 m_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 v21, inc_v,
                 buff_0,
                 f0,  inc_f );
      bli_cgemv( BLIS_CONJ_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 g0,  inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, z21 );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 z21, rs_Z );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bli_cinvscalv( BLIS_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );

      // FLA_Copy( g0, s01 );
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
    }
	else // if ( n_ahead == 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
      bli_ccopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_u );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_step_opz_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                         dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  beta;
  dcomplex  last_elem;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_TS;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &d );
  // FLA_Obj_create( datatype_A, n_A, 1, 0, 0, &e );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &f );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &g );
  dcomplex* buff_u  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_v  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_d  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_e  = ( dcomplex* ) FLA_malloc( n_A * sizeof( *buff_A ) );
  dcomplex* buff_f  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_g  = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u   = 1;
  int       inc_v   = 1;
  int       inc_d   = 1;
  int       inc_e   = 1;
  int       inc_f   = 1;
  int       inc_g   = 1;

  // FLA_Set( FLA_ZERO, Y );
  // FLA_Set( FLA_ZERO, Z );
  bli_zsetm( n_A,
             b_alg,
             buff_0,
             buff_Y, rs_Y, cs_Y );
  bli_zsetm( m_A,
             b_alg,
             buff_0,
             buff_Z, rs_Z, cs_Z );

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* a10t      = buff_A  + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A  + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a01       = buff_A  + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A  + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A  + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02       = buff_A  + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A  + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22       = buff_A  + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* y10t      = buff_Y  + (0  )*cs_Y + (i  )*rs_Y;
    dcomplex* Y20       = buff_Y  + (0  )*cs_Y + (i+1)*rs_Y;
    dcomplex* y21       = buff_Y  + (i  )*cs_Y + (i+1)*rs_Y;

    dcomplex* z10t      = buff_Z  + (0  )*cs_Z + (i  )*rs_Z;
    dcomplex* Z20       = buff_Z  + (0  )*cs_Z + (i+1)*rs_Z;
    dcomplex* z21       = buff_Z  + (i  )*cs_Z + (i+1)*rs_Z;

    dcomplex* t01       = buff_T  + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11     = buff_T  + (i  )*cs_T + (i  )*rs_T;

    dcomplex* s01       = buff_S  + (i  )*cs_S + (0  )*rs_S;
    dcomplex* sigma11   = buff_S  + (i  )*cs_S + (i  )*rs_S;

    dcomplex* u21       = buff_u  + (i+1)*inc_u;

    dcomplex* v21       = buff_v  + (i+1)*inc_v;

    dcomplex* d0        = buff_d  + (0  )*inc_d;

    dcomplex* e0        = buff_e  + (0  )*inc_e;

    dcomplex* f0        = buff_f  + (0  )*inc_f;

    dcomplex* g0        = buff_g  + (0  )*inc_g;

    dcomplex* v21_t     = v21     + (0  )*inc_v;
    dcomplex* v21_b     = v21     + (1  )*inc_v;

    dcomplex* a01_b     = a01     + (0  )*cs_A + (i-1)*rs_A;

    dcomplex* a12t_l    = a12t    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a12t_r    = a12t    + (1  )*cs_A + (0  )*rs_A;

    dcomplex* ABL       = a10t;
    dcomplex* ZBL       = z10t;

    dcomplex* a2        = alpha11;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       m_behind  = i;
    int       n_behind  = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( a01_b, last_elem );
      // FLA_Set( FLA_ONE, a01_b );
      last_elem = *a01_b;
      *a01_b = *buff_1;
    }

    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ABL, y10t, FLA_ONE, a2 );
    // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, ZBL, a01,  FLA_ONE, a2 );
    bli_zgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ABL,  rs_A, cs_A,
               y10t, cs_Y,
               buff_1,
               a2,   rs_A );
    bli_zgemv( BLIS_NO_TRANSPOSE,
               BLIS_CONJUGATE,
               m_ahead + 1,
               n_behind,
               buff_m1,
               ZBL, rs_Z, cs_Z,
               a01, rs_A,
               buff_1,
               a2,  rs_A );

    // FLA_Gemv( FLA_CONJ_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, a10t, FLA_ONE, a12t );
    // FLA_Gemv( FLA_CONJ_TRANSPOSE,    FLA_MINUS_ONE, A02, z10t, FLA_ONE, a12t );
    bli_zgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               n_ahead,
               n_behind,
               buff_m1,
               Y20,  rs_Y, cs_Y,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );
    bli_zgemv( BLIS_CONJ_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_m1,
               A02,  rs_A, cs_A,
               z10t, cs_Z,
               buff_1,
               a12t, cs_A );

    if ( m_behind > 0 )
    {
      // FLA_Copy( last_elem, a01_b );
      *a01_b = last_elem;
    }

    // FLA_Househ2_UT( FLA_LEFT,
    //                 alpha11,
    //                 a21, tau11 );
    // FLA_Copy( a21, u21 );
    FLA_Househ2_UT_l_opz( m_ahead,
                          alpha11,
                          a21, rs_A,
                          tau11 );
    bli_zcopyv( BLIS_NO_CONJUGATE,
                m_ahead,
                a21, rs_A,
                u21, inc_u );

    if ( n_ahead > 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a12t, y21 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, u21, FLA_ONE, y21 );
      bli_zcopyv( BLIS_CONJUGATE,
                  n_ahead,
                  a12t, cs_A,
                  y21,  rs_Y );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 y21, rs_Y );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ZERO, d0 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, Z20, u21, FLA_ZERO, e0 );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_0,
                 d0,  inc_d );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 Z20, rs_Z, cs_Z,
                 u21, inc_u,
                 buff_0,
                 e0,  inc_e );

      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Axpy( FLA_ONE, d0, t01 );
      bli_zcopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  buff_1,
                  d0,  inc_d,
                  t01, rs_T );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Y20, d0, FLA_ONE, y21 );
      // FLA_Gemv( FLA_TRANSPOSE,    FLA_MINUS_ONE, A02, e0, FLA_ONE, y21 );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 n_behind,
                 buff_m1,
                 Y20, rs_Y, cs_Y,
                 d0,  inc_d,
                 buff_1,
                 y21, rs_Y );
      bli_zgemv( BLIS_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_m1,
                 A02, rs_A, cs_A,
                 e0,  inc_e,
                 buff_1,
                 y21, rs_Y );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, tau11, y21 );
      bli_zinvscalv( BLIS_NO_CONJUGATE,
                     n_ahead,
                     tau11,
                     y21, rs_Y );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, FLA_MINUS_ONE, y21, a12t );
      bli_zaxpyv( BLIS_CONJUGATE,
                  n_ahead,
                  buff_m1,
                  y21,  rs_Y,
                  a12t, cs_A );

      // FLA_Househ2_UT( FLA_RIGHT, a12t_l, a12t_r, sigma11 );
      FLA_Househ2_UT_r_opz( n_ahead - 1,
                            a12t_l,
                            a12t_r, cs_A,
                            sigma11 );

      // FLA_Set( FLA_ONE, v21_t );
      // FLA_Copyt( FLA_TRANSPOSE, a12t_r, v21_b );
      *v21_t = *buff_1;
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  n_ahead - 1,
                  a12t_r, cs_A,
                  v21_b,  inc_v );

      // FLA_Dotc( FLA_CONJUGATE, y21, v21, beta );
      // FLA_Scal( FLA_MINUS_ONE, beta );
      bli_zdot( BLIS_CONJUGATE,
                n_ahead,
                y21, rs_Y,
                v21, inc_v,
                &beta );
      bli_zscals( buff_m1, &beta );

      // FLA_Copy( u21, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, v21, beta, z21 );
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  u21, inc_u,
                  z21, rs_Z );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 v21, inc_v,
                 &beta,
                 z21, rs_Z );

      // FLA_Gemvc( FLA_CONJ_TRANSPOSE,    FLA_ONE, Y20, v21, FLA_ZERO, f0 );
      // FLA_Gemvc( FLA_CONJ_NO_TRANSPOSE, FLA_ONE, A02, v21, FLA_ZERO, g0 );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 n_ahead,
                 m_behind,
                 buff_1,
                 Y20, rs_Y, cs_Y,
                 v21, inc_v,
                 buff_0,
                 f0,  inc_f );
      bli_zgemv( BLIS_CONJ_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 v21, inc_v,
                 buff_0,
                 g0,  inc_g );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, f0, FLA_ONE, z21 );
      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, Z20, g0, FLA_ONE, z21 );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 A20, rs_A, cs_A,
                 f0,  inc_f,
                 buff_1,
                 z21, rs_Z );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_m1,
                 Z20, rs_Z, cs_Z,
                 g0,  inc_g,
                 buff_1,
                 z21, rs_Z );

      // FLA_Inv_scalc( FLA_NO_CONJUGATE, sigma11, z21 );
      bli_zinvscalv( BLIS_NO_CONJUGATE,
                     m_ahead,
                     sigma11,
                     z21, rs_Z );

      // FLA_Copy( g0, s01 );
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  n_behind,
                  g0,  inc_g,
                  s01, rs_S );
    }
	else // if ( n_ahead == 0 )
    {
      // FLA_Copyt( FLA_CONJ_TRANSPOSE, a10t, t01 );
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, u21, FLA_ONE, t01 );
      bli_zcopyv( BLIS_CONJUGATE,
                  n_behind,
                  a10t, cs_A,
                  t01,  rs_T );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 u21, inc_u,
                 buff_1,
                 t01, rs_T );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &d );
  // FLA_Obj_free( &e );
  // FLA_Obj_free( &f );
  // FLA_Obj_free( &g );
  FLA_free( buff_u );
  FLA_free( buff_v );
  FLA_free( buff_d );
  FLA_free( buff_e );
  FLA_free( buff_f );
  FLA_free( buff_g );

  return FLA_SUCCESS;
}

