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

FLA_Error FLA_LU_nopiv_opt_var1( FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      FLA_LU_nopiv_ops_var1( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      FLA_LU_nopiv_opd_var1( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      FLA_LU_nopiv_opc_var1( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      FLA_LU_nopiv_opz_var1( m_A,
                             n_A,
                             buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_ops_var1( int m_A,
                                 int n_A,
                                 float* buff_A, int rs_A, int cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bli_strsv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bli_strsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_sdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  if ( m_A > n_A )
  {
    float*    ATL = buff_A;
    float*    ABL = buff_A + n_A*rs_A;

    // FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    FLA_ONE, ATL, ABL );
    bli_strsm( BLIS_RIGHT,
               BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_A - n_A,
               n_A,
               buff_1,
               ATL, rs_A, cs_A,
               ABL, rs_A, cs_A );
  }
  else if ( m_A < n_A )
  {
    float*    ATL = buff_A;
    float*    ATR = buff_A + m_A*cs_A;

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bli_strsm( BLIS_LEFT,
               BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opd_var1( int m_A,
                                 int n_A,
                                 double* buff_A, int rs_A, int cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bli_dtrsv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bli_dtrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_ddots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  if ( m_A > n_A )
  {
    double*   ATL = buff_A;
    double*   ABL = buff_A + n_A*rs_A;

    // FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    FLA_ONE, ATL, ABL );
    bli_dtrsm( BLIS_RIGHT,
               BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_A - n_A,
               n_A,
               buff_1,
               ATL, rs_A, cs_A,
               ABL, rs_A, cs_A );
  }
  else if ( m_A < n_A )
  {
    double*   ATL = buff_A;
    double*   ATR = buff_A + m_A*cs_A;

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bli_dtrsm( BLIS_LEFT,
               BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opc_var1( int m_A,
                                 int n_A,
                                 scomplex* buff_A, int rs_A, int cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bli_ctrsv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bli_ctrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_cdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  if ( m_A > n_A )
  {
    scomplex* ATL = buff_A;
    scomplex* ABL = buff_A + n_A*rs_A;

    // FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    FLA_ONE, ATL, ABL );
    bli_ctrsm( BLIS_RIGHT,
               BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_A - n_A,
               n_A,
               buff_1,
               ATL, rs_A, cs_A,
               ABL, rs_A, cs_A );
  }
  else if ( m_A < n_A )
  {
    scomplex* ATL = buff_A;
    scomplex* ATR = buff_A + m_A*cs_A;

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bli_ctrsm( BLIS_LEFT,
               BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_LU_nopiv_opz_var1( int m_A,
                                 int n_A,
                                 dcomplex* buff_A, int rs_A, int cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;

    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bli_ztrsv( BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG, A00, a10t );
    bli_ztrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a10t, cs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bli_zdots( BLIS_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    /*------------------------------------------------------------*/

  }

  if ( m_A > n_A )
  {
    dcomplex* ATL = buff_A;
    dcomplex* ABL = buff_A + n_A*rs_A;

    // FLA_Trsm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    FLA_ONE, ATL, ABL );
    bli_ztrsm( BLIS_RIGHT,
               BLIS_UPPER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_A - n_A,
               n_A,
               buff_1,
               ATL, rs_A, cs_A,
               ABL, rs_A, cs_A );
  }
  else if ( m_A < n_A )
  {
    dcomplex* ATL = buff_A;
    dcomplex* ATR = buff_A + m_A*cs_A;

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bli_ztrsm( BLIS_LEFT,
               BLIS_LOWER_TRIANGULAR,
               BLIS_NO_TRANSPOSE,
               BLIS_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return FLA_SUCCESS;
}

#endif
