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

FLA_Error FLA_Accum_T_UT_fr_opt_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_t;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_t    = FLA_Obj_vector_inc( t );

  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_t = FLA_FLOAT_PTR( t );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_Accum_T_UT_fr_ops_var1( m_A,
                                  n_A,
                                  buff_A, rs_A, cs_A,
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_t = FLA_DOUBLE_PTR( t );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_Accum_T_UT_fr_opd_var1( m_A,
                                  n_A,
                                  buff_A, rs_A, cs_A,
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_t = FLA_COMPLEX_PTR( t );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_Accum_T_UT_fr_opc_var1( m_A,
                                  n_A,
                                  buff_A, rs_A, cs_A,
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_t = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_Accum_T_UT_fr_opz_var1( m_A,
                                  n_A,
                                  buff_A, rs_A, cs_A,
                                  buff_t, inc_t,
                                  buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fr_ops_var1( int m_A,
                                      int n_A,
                                      float* buff_A, int rs_A, int cs_A,
                                      float* buff_t, int inc_t,
                                      float* buff_T, int rs_T, int cs_T )
{
  float* buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  int    i;

  for ( i = 0; i < m_A; ++i )
  {
    float* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    float* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float* tau1     = buff_t + (i  )*inc_t;

    float* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    float* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int    n_ahead  = n_A - i - 1;
    int    m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bli_scopyv( BLIS_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );
  
    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bli_sgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fr_opd_var1( int m_A,
                                      int n_A,
                                      double* buff_A, int rs_A, int cs_A,
                                      double* buff_t, int inc_t,
                                      double* buff_T, int rs_T, int cs_T )
{
  double* buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  int     i;

  for ( i = 0; i < m_A; ++i )
  {
    double* a01     = buff_A + (i  )*cs_A + (0  )*rs_A;
    double* A02     = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double* a12t    = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double* tau1    = buff_t + (i  )*inc_t;

    double* tau11   = buff_T + (i  )*cs_T + (i  )*rs_T;
    double* t01     = buff_T + (i  )*cs_T + (0  )*rs_T;

    int    n_ahead  = n_A - i - 1;
    int    m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bli_dcopyv( BLIS_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bli_dgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fr_opc_var1( int m_A,
                                      int n_A,
                                      scomplex* buff_A, int rs_A, int cs_A,
                                      scomplex* buff_t, int inc_t,
                                      scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < m_A; ++i )
  {
    scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    scomplex* tau1     = buff_t + (i  )*inc_t;

    scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bli_ccopyv( BLIS_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bli_cgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Accum_T_UT_fr_opz_var1( int m_A,
                                      int n_A,
                                      dcomplex* buff_A, int rs_A, int cs_A,
                                      dcomplex* buff_t, int inc_t,
                                      dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  int       i;

  for ( i = 0; i < m_A; ++i )
  {
    dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* A02      = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    dcomplex* tau1     = buff_t + (i  )*inc_t;

    dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* t01      = buff_T + (i  )*cs_T + (0  )*rs_T;

    int       n_ahead  = n_A - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( tau1, tau11 );
    *tau11 = *tau1;

    // FLA_Copyt_external( FLA_CONJ_NO_TRANSPOSE, a01, t01 );
    bli_zcopyv( BLIS_CONJUGATE,
                m_behind,
                a01, rs_A,
                t01, rs_T );

    // FLA_Gemvc_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_CONJUGATE,
    //                     FLA_ONE, A02, a12t, FLA_ONE, t01 );
    bli_zgemv( BLIS_CONJ_NO_TRANSPOSE,
               BLIS_NO_CONJUGATE,
               m_behind,
               n_ahead,
               buff_1,
               A02, rs_A, cs_A,
               a12t, cs_A,
               buff_1,
               t01, rs_T );

    /*------------------------------------------------------------*/

  }


  return FLA_SUCCESS;
}

