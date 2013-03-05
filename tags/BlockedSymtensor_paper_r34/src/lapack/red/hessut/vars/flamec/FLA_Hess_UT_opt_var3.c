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

FLA_Error FLA_Hess_UT_opt_var3( FLA_Obj A, FLA_Obj T )
{
  return FLA_Hess_UT_step_opt_var3( A, T );
}

FLA_Error FLA_Hess_UT_step_opt_var3( FLA_Obj A, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, m_T;
  int          rs_A, cs_A;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  m_T      = FLA_Obj_length( T );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Hess_UT_step_ops_var3( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Hess_UT_step_opd_var3( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opc_var3( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Hess_UT_step_opz_var3( m_A,
                                 m_T,
                                 buff_A, rs_A, cs_A,
                                 buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_ops_var3( int m_A,
                                     int m_T,
                                     float* buff_A, int rs_A, int cs_A, 
                                     float* buff_T, int rs_T, int cs_T )
{
  float*    buff_2  = FLA_FLOAT_PTR( FLA_TWO );
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );

  float     first_elem;
  float     dot_product;
  float     beta, conj_beta;
  float     inv_tau11;
  float     minus_inv_tau11;
  float     minus_upsilon1, minus_conj_upsilon1;
  float     minus_psi1, minus_conj_psi1;
  float     minus_zeta1;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  float*    buff_u = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_y = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_z = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_v = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  float*    buff_w = ( float* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_y  = 1;
  int       inc_z  = 1;
  int       inc_v  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem           = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    float*    upsilon1 = buff_u + (i  )*inc_u;
    float*    u2       = buff_u + (i+1)*inc_u;

    float*    y0       = buff_y + (0  )*inc_y;
    float*    psi1     = buff_y + (i  )*inc_y;
    float*    y2       = buff_y + (i+1)*inc_y;

    float*    zeta1    = buff_z + (i  )*inc_z;
    float*    z2       = buff_z + (i+1)*inc_z;

    float*    v2       = buff_v + (i+1)*inc_v;

    float*    w2       = buff_w + (i+1)*inc_w;

    float*    a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    float*    a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon1, minus_upsilon1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, minus_conj_upsilon1 );
      bli_smult3( buff_m1, upsilon1, &minus_upsilon1 );
      bli_scopyconj( &minus_upsilon1, &minus_conj_upsilon1 );

      // FLA_Copy( psi1, minus_psi1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_psi1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_psi1, minus_conj_psi1 );
      bli_smult3( buff_m1, psi1, &minus_psi1 );
      bli_scopyconj( &minus_psi1, &minus_conj_psi1 );

      // FLA_Copy( zeta1, minus_zeta1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta1 );
      bli_smult3( buff_m1, zeta1, &minus_zeta1 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, psi1,     alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta1,    upsilon1, alpha11 );
      bli_saxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_upsilon1,
                  psi1,    1,
                  alpha11, 1 );
      bli_saxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_zeta1,
                  upsilon1, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon1, y2, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta1,    u2, a12t );
      bli_saxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_upsilon1,
                  y2,   inc_y,
                  a12t, cs_A );
      bli_saxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_zeta1,
                  u2,   inc_u,
                  a12t, cs_A );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi1,     u2, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon1, z2, a21 );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi1,
                  u2,  inc_u,
                  a21, rs_A );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon1,
                  z2,  inc_z,
                  a21, rs_A );
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_sger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_sger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, v2 );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 v2, inc_v );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, w2 );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w2, inc_w );

      // FLA_Copy( a21, u2 );
      // FLA_Copy( v2, y2 );
      // FLA_Copy( w2, z2 );
      bli_scopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u2,  inc_u );
      bli_scopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  v2, inc_v,
                  y2, inc_y );
      bli_scopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  w2, inc_w,
                  z2, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bli_sdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bli_sinvscals( buff_2, &beta );
      bli_scopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bli_sscals( &minus_inv_tau11, &conj_beta );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bli_sscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bli_sscals( &minus_inv_tau11, &beta );
      bli_saxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bli_sscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bli_sdot( BLIS_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bli_sscals( &minus_inv_tau11, &dot_product );
      bli_saxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bli_sgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bli_sger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bli_sgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_sger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_sger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_y );
  FLA_free( buff_z );
  FLA_free( buff_v );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opd_var3( int m_A,
                                     int m_T,
                                     double* buff_A, int rs_A, int cs_A, 
                                     double* buff_T, int rs_T, int cs_T )
{
  double*   buff_2  = FLA_DOUBLE_PTR( FLA_TWO );
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );

  double    first_elem;
  double    dot_product;
  double    beta, conj_beta;
  double    inv_tau11;
  double    minus_inv_tau11;
  double    minus_upsilon1, minus_conj_upsilon1;
  double    minus_psi1, minus_conj_psi1;
  double    minus_zeta1;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  double*   buff_u = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_y = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_z = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_v = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  double*   buff_w = ( double* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_y  = 1;
  int       inc_z  = 1;
  int       inc_v  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem           = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    double*   upsilon1 = buff_u + (i  )*inc_u;
    double*   u2       = buff_u + (i+1)*inc_u;

    double*   y0       = buff_y + (0  )*inc_y;
    double*   psi1     = buff_y + (i  )*inc_y;
    double*   y2       = buff_y + (i+1)*inc_y;

    double*   zeta1    = buff_z + (i  )*inc_z;
    double*   z2       = buff_z + (i+1)*inc_z;

    double*   v2       = buff_v + (i+1)*inc_v;

    double*   w2       = buff_w + (i+1)*inc_w;

    double*   a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    double*   a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon1, minus_upsilon1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, minus_conj_upsilon1 );
      bli_dmult3( buff_m1, upsilon1, &minus_upsilon1 );
      bli_dcopyconj( &minus_upsilon1, &minus_conj_upsilon1 );

      // FLA_Copy( psi1, minus_psi1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_psi1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_psi1, minus_conj_psi1 );
      bli_dmult3( buff_m1, psi1, &minus_psi1 );
      bli_dcopyconj( &minus_psi1, &minus_conj_psi1 );

      // FLA_Copy( zeta1, minus_zeta1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta1 );
      bli_dmult3( buff_m1, zeta1, &minus_zeta1 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, psi1,     alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta1,    upsilon1, alpha11 );
      bli_daxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_upsilon1,
                  psi1,    1,
                  alpha11, 1 );
      bli_daxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_zeta1,
                  upsilon1, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon1, y2, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta1,    u2, a12t );
      bli_daxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_upsilon1,
                  y2,   inc_y,
                  a12t, cs_A );
      bli_daxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_zeta1,
                  u2,   inc_u,
                  a12t, cs_A );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi1,     u2, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon1, z2, a21 );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi1,
                  u2,  inc_u,
                  a21, rs_A );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon1,
                  z2,  inc_z,
                  a21, rs_A );
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_dger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_dger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, v2 );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 v2, inc_v );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, w2 );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w2, inc_w );

      // FLA_Copy( a21, u2 );
      // FLA_Copy( v2, y2 );
      // FLA_Copy( w2, z2 );
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u2,  inc_u );
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  v2, inc_v,
                  y2, inc_y );
      bli_dcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  w2, inc_w,
                  z2, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bli_ddot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bli_dinvscals( buff_2, &beta );
      bli_dcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bli_dscals( &minus_inv_tau11, &conj_beta );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bli_dscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bli_dscals( &minus_inv_tau11, &beta );
      bli_daxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bli_dscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bli_ddot( BLIS_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bli_dscals( &minus_inv_tau11, &dot_product );
      bli_daxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bli_dgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bli_dger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bli_dgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_dger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_dger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_y );
  FLA_free( buff_z );
  FLA_free( buff_v );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opc_var3( int m_A,
                                     int m_T,
                                     scomplex* buff_A, int rs_A, int cs_A, 
                                     scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_2  = FLA_COMPLEX_PTR( FLA_TWO );
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );

  scomplex  first_elem;
  scomplex  dot_product;
  scomplex  beta, conj_beta;
  scomplex  inv_tau11;
  scomplex  minus_inv_tau11;
  scomplex  minus_upsilon1, minus_conj_upsilon1;
  scomplex  minus_psi1, minus_conj_psi1;
  scomplex  minus_zeta1;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  scomplex* buff_u = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_y = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_z = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_v = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  scomplex* buff_w = ( scomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_y  = 1;
  int       inc_z  = 1;
  int       inc_v  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem           = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    scomplex* upsilon1 = buff_u + (i  )*inc_u;
    scomplex* u2       = buff_u + (i+1)*inc_u;

    scomplex* y0       = buff_y + (0  )*inc_y;
    scomplex* psi1     = buff_y + (i  )*inc_y;
    scomplex* y2       = buff_y + (i+1)*inc_y;

    scomplex* zeta1    = buff_z + (i  )*inc_z;
    scomplex* z2       = buff_z + (i+1)*inc_z;

    scomplex* v2       = buff_v + (i+1)*inc_v;

    scomplex* w2       = buff_w + (i+1)*inc_w;

    scomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    scomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon1, minus_upsilon1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, minus_conj_upsilon1 );
      bli_cmult3( buff_m1, upsilon1, &minus_upsilon1 );
      bli_ccopyconj( &minus_upsilon1, &minus_conj_upsilon1 );

      // FLA_Copy( psi1, minus_psi1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_psi1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_psi1, minus_conj_psi1 );
      bli_cmult3( buff_m1, psi1, &minus_psi1 );
      bli_ccopyconj( &minus_psi1, &minus_conj_psi1 );

      // FLA_Copy( zeta1, minus_zeta1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta1 );
      bli_cmult3( buff_m1, zeta1, &minus_zeta1 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, psi1,     alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta1,    upsilon1, alpha11 );
      bli_caxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_upsilon1,
                  psi1,    1,
                  alpha11, 1 );
      bli_caxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_zeta1,
                  upsilon1, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon1, y2, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta1,    u2, a12t );
      bli_caxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_upsilon1,
                  y2,   inc_y,
                  a12t, cs_A );
      bli_caxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_zeta1,
                  u2,   inc_u,
                  a12t, cs_A );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi1,     u2, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon1, z2, a21 );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi1,
                  u2,  inc_u,
                  a21, rs_A );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon1,
                  z2,  inc_z,
                  a21, rs_A );
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_cger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_cger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, v2 );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 v2, inc_v );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, w2 );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w2, inc_w );

      // FLA_Copy( a21, u2 );
      // FLA_Copy( v2, y2 );
      // FLA_Copy( w2, z2 );
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u2,  inc_u );
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  v2, inc_v,
                  y2, inc_y );
      bli_ccopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  w2, inc_w,
                  z2, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bli_cdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bli_cinvscals( buff_2, &beta );
      bli_ccopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bli_cscals( &minus_inv_tau11, &conj_beta );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bli_cscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bli_cscals( &minus_inv_tau11, &beta );
      bli_caxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bli_cscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bli_cdot( BLIS_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bli_cscals( &minus_inv_tau11, &dot_product );
      bli_caxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bli_cgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bli_cger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bli_cgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_cger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_cger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_y );
  FLA_free( buff_z );
  FLA_free( buff_v );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}



FLA_Error FLA_Hess_UT_step_opz_var3( int m_A,
                                     int m_T,
                                     dcomplex* buff_A, int rs_A, int cs_A, 
                                     dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_2  = FLA_DOUBLE_COMPLEX_PTR( FLA_TWO );
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_0  = FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );

  dcomplex  first_elem;
  dcomplex  dot_product;
  dcomplex  beta, conj_beta;
  dcomplex  inv_tau11;
  dcomplex  minus_inv_tau11;
  dcomplex  minus_upsilon1, minus_conj_upsilon1;
  dcomplex  minus_psi1, minus_conj_psi1;
  dcomplex  minus_zeta1;
  int       i;

  // b_alg = FLA_Obj_length( T );
  int       b_alg = m_T;

  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &u );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &y );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &z );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &v );
  // FLA_Obj_create( datatype_A, m_A, 1, 0, 0, &w );
  dcomplex* buff_u = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_y = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_z = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_v = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  dcomplex* buff_w = ( dcomplex* ) FLA_malloc( m_A * sizeof( *buff_A ) );
  int       inc_u  = 1;
  int       inc_y  = 1;
  int       inc_z  = 1;
  int       inc_v  = 1;
  int       inc_w  = 1;

  // Initialize some variables (only to prevent compiler warnings).
  first_elem           = *buff_0;
  minus_inv_tau11 = *buff_0;

  for ( i = 0; i < b_alg; ++i )
  {
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;
    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

    dcomplex* upsilon1 = buff_u + (i  )*inc_u;
    dcomplex* u2       = buff_u + (i+1)*inc_u;

    dcomplex* y0       = buff_y + (0  )*inc_y;
    dcomplex* psi1     = buff_y + (i  )*inc_y;
    dcomplex* y2       = buff_y + (i+1)*inc_y;

    dcomplex* zeta1    = buff_z + (i  )*inc_z;
    dcomplex* z2       = buff_z + (i+1)*inc_z;

    dcomplex* v2       = buff_v + (i+1)*inc_v;

    dcomplex* w2       = buff_w + (i+1)*inc_w;

    dcomplex* a21_t    = a21    + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a21_b    = a21    + (0  )*cs_A + (1  )*rs_A;

    int       m_ahead  = m_A - i - 1;
    int       n_ahead  = m_A - i - 1;
    int       m_behind = i;
    int       n_behind = i;

    /*------------------------------------------------------------*/

    if ( m_behind > 0 )
    {
      // FLA_Copy( upsilon1, minus_upsilon1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_upsilon1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, minus_conj_upsilon1 );
      bli_zmult3( buff_m1, upsilon1, &minus_upsilon1 );
      bli_zcopyconj( &minus_upsilon1, &minus_conj_upsilon1 );

      // FLA_Copy( psi1, minus_psi1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_psi1 );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, minus_psi1, minus_conj_psi1 );
      bli_zmult3( buff_m1, psi1, &minus_psi1 );
      bli_zcopyconj( &minus_psi1, &minus_conj_psi1 );

      // FLA_Copy( zeta1, minus_zeta1 );
      // FLA_Scal( FLA_MINUS_ONE, minus_zeta1 );
      bli_zmult3( buff_m1, zeta1, &minus_zeta1 );

      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_upsilon1, psi1,     alpha11 );
      // FLA_Axpyt( FLA_CONJ_NO_TRANSPOSE, minus_zeta1,    upsilon1, alpha11 );
      bli_zaxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_upsilon1,
                  psi1,    1,
                  alpha11, 1 );
      bli_zaxpyv( BLIS_CONJUGATE,
                  1,
                  &minus_zeta1,
                  upsilon1, 1,
                  alpha11,  1 );

      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_upsilon1, y2, a12t );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, minus_zeta1,    u2, a12t );
      bli_zaxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_upsilon1,
                  y2,   inc_y,
                  a12t, cs_A );
      bli_zaxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &minus_zeta1,
                  u2,   inc_u,
                  a12t, cs_A );

      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_psi1,     u2, a21 );
      // FLA_Axpyt( FLA_NO_TRANSPOSE, minus_conj_upsilon1, z2, a21 );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_psi1,
                  u2,  inc_u,
                  a21, rs_A );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &minus_conj_upsilon1,
                  z2,  inc_z,
                  a21, rs_A );
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
    }

    if ( m_behind > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_zger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_zger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    if ( m_ahead > 0 )
    {
      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, v2 );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 v2, inc_v );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A22, a21, FLA_ZERO, w2 );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_ahead,
                 buff_1,
                 A22, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 w2, inc_w );

      // FLA_Copy( a21, u2 );
      // FLA_Copy( v2, y2 );
      // FLA_Copy( w2, z2 );
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  a21, rs_A,
                  u2,  inc_u );
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  v2, inc_v,
                  y2, inc_y );
      bli_zcopyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  w2, inc_w,
                  z2, inc_z );

      // FLA_Dotc( FLA_CONJUGATE, a21, z2, beta );
      // FLA_Inv_scal( FLA_TWO, beta );
      // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, beta, conj_beta );
      bli_zdot( BLIS_CONJUGATE,
                m_ahead,
                a21, rs_A,
                z2,  inc_z,
                &beta );
      bli_zinvscals( buff_2, &beta );
      bli_zcopyconj( &beta, &conj_beta );

      // FLA_Scal( minus_inv_tau11, conj_beta );
      // FLA_Axpy( conj_beta, a21, y2 );
      // FLA_Scal( inv_tau11, y2 );
      bli_zscals( &minus_inv_tau11, &conj_beta );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &conj_beta,
                  a21, rs_A,
                  y2, inc_y );
      bli_zscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  y2, inc_y );

      // FLA_Scal( minus_inv_tau11, beta );
      // FLA_Axpy( beta, a21, z2 );
      // FLA_Scal( inv_tau11, z2 );
      bli_zscals( &minus_inv_tau11, &beta );
      bli_zaxpyv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &beta,
                  a21, rs_A,
                  z2, inc_z );
      bli_zscalv( BLIS_NO_CONJUGATE,
                  m_ahead,
                  &inv_tau11,
                  z2, inc_z );

      // FLA_Dot( a12t, a21, dot_product );
      // FLA_Scal( minus_inv_tau11, dot_product );
      // FLA_Axpyt( FLA_CONJ_TRANSPOSE, dot_product, a21, a12t );
      bli_zdot( BLIS_NO_CONJUGATE,
                m_ahead,
                a12t, cs_A,
                a21,  rs_A,
                &dot_product );
      bli_zscals( &minus_inv_tau11, &dot_product );
      bli_zaxpyv( BLIS_CONJUGATE,
                  m_ahead,
                  &dot_product,
                  a21,  rs_A,
                  a12t, cs_A );

      // FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A02, a21, FLA_ZERO, y0 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, minus_inv_tau11, y0, a21, A02 );
      bli_zgemv( BLIS_NO_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_behind,
                 n_ahead,
                 buff_1,
                 A02, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 y0,  inc_y );
      bli_zger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_behind,
                n_ahead,
                &minus_inv_tau11,
                y0,  inc_y,
                a21, rs_A,
                A02, rs_A, cs_A );

      // FLA_Gemv( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, FLA_ZERO, t01 );
      bli_zgemv( BLIS_CONJ_TRANSPOSE,
                 BLIS_NO_CONJUGATE,
                 m_ahead,
                 n_behind,
                 buff_1,
                 A20, rs_A, cs_A,
                 a21, rs_A,
                 buff_0,
                 t01, rs_T );

      // FLA_Copy( first_elem, a21_t );
      *a21_t = first_elem;
    }

    if ( m_behind + 1 == b_alg && m_ahead > 0 )
    {
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, u2, y2, A22 );
      // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, z2, u2, A22 );
      bli_zger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                u2,  inc_u,
                y2,  inc_y,
                A22, rs_A, cs_A );
      bli_zger( BLIS_NO_CONJUGATE,
                BLIS_CONJUGATE,
                m_ahead,
                n_ahead,
                buff_m1,
                z2,  inc_z,
                u2,  inc_u,
                A22, rs_A, cs_A );
    }

    /*------------------------------------------------------------*/

  }

  // FLA_Obj_free( &u );
  // FLA_Obj_free( &y );
  // FLA_Obj_free( &z );
  // FLA_Obj_free( &v );
  // FLA_Obj_free( &w );
  FLA_free( buff_u );
  FLA_free( buff_y );
  FLA_free( buff_z );
  FLA_free( buff_v );
  FLA_free( buff_w );

  return FLA_SUCCESS;
}

