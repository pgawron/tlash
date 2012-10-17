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

FLA_Error FLA_Eig_gest_nl_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_AB;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          inc_y;
  FLA_Obj      yT, yB;

  datatype = FLA_Obj_datatype( A );

  m_AB     = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
 
  FLA_Part_2x1( Y,    &yT,
                      &yB,     1, FLA_TOP );

  inc_y    = FLA_Obj_vector_inc( yT );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_y = FLA_FLOAT_PTR( yT );
      float* buff_B = FLA_FLOAT_PTR( B );

      FLA_Eig_gest_nl_ops_var4( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_y = FLA_DOUBLE_PTR( yT );
      double* buff_B = FLA_DOUBLE_PTR( B );

      FLA_Eig_gest_nl_opd_var4( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_y = FLA_COMPLEX_PTR( yT );
      scomplex* buff_B = FLA_COMPLEX_PTR( B );

      FLA_Eig_gest_nl_opc_var4( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( yT );
      dcomplex* buff_B = FLA_DOUBLE_COMPLEX_PTR( B );

      FLA_Eig_gest_nl_opz_var4( m_AB,
                                buff_A, rs_A, cs_A,
                                buff_y, inc_y,
                                buff_B, rs_B, cs_B );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_ops_var4( int m_AB,
                                    float* buff_A, int rs_A, int cs_A, 
                                    float* buff_y, int inc_y, 
                                    float* buff_B, int rs_B, int cs_B )
{
  float*    buff_1   = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_1h  = FLA_FLOAT_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    float*    A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float*    b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    float*    beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    float     psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bli_smult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bli_sher2( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_saxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bli_sscalv( BLIS_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_sscals( beta11, alpha11 );
    bli_sscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, a21, b10t, A20 );
    bli_sger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              a21,  rs_A,
              b10t, cs_B,
              A20,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a21 );
    bli_sscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opd_var4( int m_AB,
                                    double* buff_A, int rs_A, int cs_A, 
                                    double* buff_y, int inc_y, 
                                    double* buff_B, int rs_B, int cs_B )
{
  double*   buff_1   = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_1h  = FLA_DOUBLE_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    double*   A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double*   b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    double*   beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    double    psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bli_dmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bli_dher2( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bli_dscalv( BLIS_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_dscals( beta11, alpha11 );
    bli_dscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, a21, b10t, A20 );
    bli_dger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              a21,  rs_A,
              b10t, cs_B,
              A20,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a21 );
    bli_dscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opc_var4( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B )
{
  scomplex* buff_1   = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_1h  = FLA_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    scomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    scomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    scomplex  psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bli_cmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bli_cher2( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_caxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bli_cscalv( BLIS_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_cscals( beta11, alpha11 );
    bli_cscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, a21, b10t, A20 );
    bli_cger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              a21,  rs_A,
              b10t, cs_B,
              A20,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a21 );
    bli_cscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Eig_gest_nl_opz_var4( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B )
{
  dcomplex* buff_1   = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_1h  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  int       i;

  for ( i = 0; i < m_AB; ++i )
  {
    dcomplex* A00      = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t     = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20      = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* b10t     = buff_B + (0  )*cs_B + (i  )*rs_B;
    dcomplex* beta11   = buff_B + (i  )*cs_B + (i  )*rs_B;

    dcomplex  psi11;

    int       m_ahead  = m_AB - i - 1;
    int       m_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Copy_external( alpha11, psi11 );
    // FLA_Scal_external( FLA_ONE_HALF, psi11 );
    bli_zmult3( buff_1h, alpha11, &psi11 );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Her2c_external( FLA_LOWER_TRIANGULAR, FLA_CONJUGATE,
    //                     FLA_ONE, a10t, b10t, A00 );
    bli_zher2( BLIS_LOWER_TRIANGULAR,
               BLIS_CONJUGATE,
               m_behind,
               buff_1,
               a10t, cs_A,
               b10t, cs_B,
               A00,  rs_A, cs_A );

    // FLA_Axpy_external( psi11, b10t, a10t );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_behind,
                &psi11,
                b10t, cs_B,
                a10t, cs_A );

    // FLA_Scal_external( beta11, a10t );
    bli_zscalv( BLIS_NO_CONJUGATE,
                m_behind,
                beta11,
                a10t, cs_A );

    // FLA_Scal_external( beta11, alpha11 );
    // FLA_Scal_external( beta11, alpha11 );
    bli_zscals( beta11, alpha11 );
    bli_zscals( beta11, alpha11 );

    // FLA_Ger_external( FLA_ONE, a21, b10t, A20 );
    bli_zger( BLIS_NO_CONJUGATE,
              BLIS_NO_CONJUGATE,
              m_ahead,
              m_behind,
              buff_1,
              a21,  rs_A,
              b10t, cs_B,
              A20,  rs_A, cs_A );

    // FLA_Scal_external( beta11, a21 );
    bli_zscalv( BLIS_NO_CONJUGATE,
                m_ahead,
                beta11,
                a21, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

