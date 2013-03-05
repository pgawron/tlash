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

FLA_Error FLA_Eig_gest_iu_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
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

      FLA_Eig_gest_iu_ops_var5( m_AB,
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

      FLA_Eig_gest_iu_opd_var5( m_AB,
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

      FLA_Eig_gest_iu_opc_var5( m_AB,
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

      FLA_Eig_gest_iu_opz_var5( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_ops_var5( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_m1  = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float*    buff_m1h = FLA_FLOAT_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    float*    b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    float*    B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    float     psi11;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_sinvscals( beta11, alpha11 );
    bli_sinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bli_smult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_sinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bli_sher2( BLIS_UPPER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bli_strsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opd_var5( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_m1  = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  double*   buff_m1h = FLA_DOUBLE_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    double*   b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    double*   B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    double    psi11;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_dinvscals( beta11, alpha11 );
    bli_dinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bli_dmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_dinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bli_dher2( BLIS_UPPER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bli_dtrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opc_var5( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_m1  = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex* buff_m1h = FLA_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    scomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    scomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    scomplex  psi11;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_cinvscals( beta11, alpha11 );
    bli_cinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bli_cmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_cinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bli_cher2( BLIS_UPPER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bli_ctrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_iu_opz_var5( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_m1  = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  dcomplex* buff_m1h = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
    dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;
    dcomplex* b12t     = buff_B + (i+1)*cs_B + (i  )*rs_B;
    dcomplex* B22      = buff_B + (i+1)*cs_B + (i+1)*rs_B;

    dcomplex  psi11;

    int       m_ahead  = m_AB - i - 1;

    /*------------------------------------------------------------*/

    // FLA_Inv_scal_external( beta11, alpha11 );
    // FLA_Inv_scal_external( beta11, alpha11 );
    bli_zinvscals( beta11, alpha11 );
    bli_zinvscals( beta11, alpha11 );

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_MINUS_ONE_HALF, psi11 );
    bli_zmult3( buff_m1h, alpha11, &psi11 );

    // FLA_Inv_scal_external( beta11, a12t );
    bli_zinvscalv( BLIS_NO_CONJUGATE,
                   m_ahead,
                   beta11,
                   a12t, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Her2c_external( FLA_UPPER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_MINUS_ONE, a12t, b12t, A22 );
    bli_zher2( BLIS_UPPER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_ahead,
               buff_m1,
               a12t, cs_A,
               b12t, cs_B,
               A22,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b12t, a12t );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_ahead,
                &psi11,
                b12t, cs_B,
                a12t, cs_A );

    // FLA_Trsv_external( FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
    //                    B22, a12t );
    bli_ztrsv( BLIS_UPPER_TRIANGULAR,
               BLIS_TRANSPOSE,
               BLIS_NONUNIT_DIAG,
               m_ahead,
               B22,  rs_B, cs_B,
               a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

