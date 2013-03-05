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

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_LU_piv_opt_var4( FLA_Obj A, FLA_Obj p )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_p;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      int*   buff_p = FLA_INT_PTR( p );

      FLA_LU_piv_ops_var4( m_A,
                           n_A,
                           buff_A, rs_A, cs_A,
                           buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      int*    buff_p = FLA_INT_PTR( p );

      FLA_LU_piv_opd_var4( m_A,
                           n_A,
                           buff_A, rs_A, cs_A,
                           buff_p, inc_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      int*      buff_p = FLA_INT_PTR( p );

      FLA_LU_piv_opc_var4( m_A,
                           n_A,
                           buff_A, rs_A, cs_A,
                           buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      int*      buff_p = FLA_INT_PTR( p );

      FLA_LU_piv_opz_var4( m_A,
                           n_A,
                           buff_A, rs_A, cs_A,
                           buff_p, inc_p );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_piv_ops_var4( int m_A,
                               int n_A,
                               float*    buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_sdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bli_sgemv( BLIS_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bli_samax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
    FLA_Apply_pivots_ln_ops_var1( 1,
                                  alpha11, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a10t,
    //                A20,      &AB0 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
    FLA_Apply_pivots_ln_ops_var1( mn_behind,
                                  a10t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a12t,
    //                A22,      &AB2 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
    FLA_Apply_pivots_ln_ops_var1( n_ahead,
                                  a12t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bli_sgemv( BLIS_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( alpha11, a21 );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_piv_opd_var4( int m_A,
                               int n_A,
                               double*   buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_ddots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bli_dgemv( BLIS_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bli_damax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
    FLA_Apply_pivots_ln_opd_var1( 1,
                                  alpha11, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a10t,
    //                A20,      &AB0 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
    FLA_Apply_pivots_ln_opd_var1( mn_behind,
                                  a10t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a12t,
    //                A22,      &AB2 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
    FLA_Apply_pivots_ln_opd_var1( n_ahead,
                                  a12t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bli_dgemv( BLIS_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( alpha11, a21 );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_piv_opc_var4( int m_A,
                               int n_A,
                               scomplex* buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_cdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bli_cgemv( BLIS_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bli_camax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
    FLA_Apply_pivots_ln_opc_var1( 1,
                                  alpha11, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a10t,
    //                A20,      &AB0 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
    FLA_Apply_pivots_ln_opc_var1( mn_behind,
                                  a10t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a12t,
    //                A22,      &AB2 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
    FLA_Apply_pivots_ln_opc_var1( n_ahead,
                                  a12t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bli_cgemv( BLIS_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( alpha11, a21 );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_piv_opz_var4( int m_A,
                               int n_A,
                               dcomplex* buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       n_ahead   = n_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_zdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bli_zgemv( BLIS_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bli_zamax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
    FLA_Apply_pivots_ln_opz_var1( 1,
                                  alpha11, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a10t,
    //                A20,      &AB0 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
    FLA_Apply_pivots_ln_opz_var1( mn_behind,
                                  a10t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Merge_2x1( a12t,
    //                A22,      &AB2 );

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
    FLA_Apply_pivots_ln_opz_var1( n_ahead,
                                  a12t, rs_A, cs_A,
                                  0,
                                  0,
                                  pi1, inc_p );

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bli_zgemv( BLIS_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    // FLA_Inv_scal_external( alpha11, a21 );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   alpha11,
                   a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

#endif
