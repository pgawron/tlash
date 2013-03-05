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

FLA_Error FLA_UDdate_UT_opt_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T )
{
  FLA_Datatype datatype;
  int          mn_RT, m_C, m_D;
  int          rs_R, cs_R;
  int          rs_C, cs_C;
  int          rs_D, cs_D;
  int          rs_T, cs_T;

  datatype = FLA_Obj_datatype( R );

  mn_RT    = FLA_Obj_length( R );
  m_C      = FLA_Obj_length( C );
  m_D      = FLA_Obj_length( D );

  rs_R     = FLA_Obj_row_stride( R );
  cs_R     = FLA_Obj_col_stride( R );
  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );
  rs_D     = FLA_Obj_row_stride( D );
  cs_D     = FLA_Obj_col_stride( D );
  rs_T     = FLA_Obj_row_stride( T );
  cs_T     = FLA_Obj_col_stride( T );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_R = FLA_FLOAT_PTR( R );
      float* buff_C = FLA_FLOAT_PTR( C );
      float* buff_D = FLA_FLOAT_PTR( D );
      float* buff_T = FLA_FLOAT_PTR( T );

      FLA_UDdate_UT_ops_var1( mn_RT,
                              m_C,
                              m_D,
                              buff_R, rs_R, cs_R,
                              buff_C, rs_C, cs_C,
                              buff_D, rs_D, cs_D,
                              buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_R = FLA_DOUBLE_PTR( R );
      double* buff_C = FLA_DOUBLE_PTR( C );
      double* buff_D = FLA_DOUBLE_PTR( D );
      double* buff_T = FLA_DOUBLE_PTR( T );

      FLA_UDdate_UT_opd_var1( mn_RT,
                              m_C,
                              m_D,
                              buff_R, rs_R, cs_R,
                              buff_C, rs_C, cs_C,
                              buff_D, rs_D, cs_D,
                              buff_T, rs_T, cs_T );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_R = FLA_COMPLEX_PTR( R );
      scomplex* buff_C = FLA_COMPLEX_PTR( C );
      scomplex* buff_D = FLA_COMPLEX_PTR( D );
      scomplex* buff_T = FLA_COMPLEX_PTR( T );

      FLA_UDdate_UT_opc_var1( mn_RT,
                              m_C,
                              m_D,
                              buff_R, rs_R, cs_R,
                              buff_C, rs_C, cs_C,
                              buff_D, rs_D, cs_D,
                              buff_T, rs_T, cs_T );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_R = FLA_DOUBLE_COMPLEX_PTR( R );
      dcomplex* buff_C = FLA_DOUBLE_COMPLEX_PTR( C );
      dcomplex* buff_D = FLA_DOUBLE_COMPLEX_PTR( D );
      dcomplex* buff_T = FLA_DOUBLE_COMPLEX_PTR( T );

      FLA_UDdate_UT_opz_var1( mn_RT,
                              m_C,
                              m_D,
                              buff_R, rs_R, cs_R,
                              buff_C, rs_C, cs_C,
                              buff_D, rs_D, cs_D,
                              buff_T, rs_T, cs_T );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_UDdate_UT_ops_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  float* buff_R, int rs_R, int cs_R,
                                  float* buff_C, int rs_C, int cs_C,
                                  float* buff_D, int rs_D, int cs_D,
                                  float* buff_T, int rs_T, int cs_T )
{
  float*    buff_half = FLA_FLOAT_PTR( FLA_ONE_HALF );
  float*    buff_1    = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1   = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_RT; ++i )
  {
    float*    rho11     = buff_R + (i  )*cs_R + (i  )*rs_R;
    float*    r12t      = buff_R + (i+1)*cs_R + (i  )*rs_R;

    float*    c1        = buff_C + (i  )*cs_C + (0  )*rs_C;
    float*    C2        = buff_C + (i+1)*cs_C + (0  )*rs_C;

    float*    d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    float*    D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    float*    tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    float*    w12t      = buff_T + (i+1)*cs_T + (i  )*rs_T;

    int       mn_ahead  = mn_RT - i - 1;

    //------------------------------------------------------------//

    // FLA_Househ3UD_UT( rho11,
    //                   c1,
    //                   d1, tau11 );
    FLA_Househ3UD_UT_ops( m_C,
                          m_D,
                          rho11,
                          c1, rs_C,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_HUD_UT( FLA_LEFT,
    //                   tau11, r12t,
    //                   c1,    C2,
    //                   d1,    D2 );
    FLA_Apply_HUD_UT_l_ops_var1( m_C,
                                 m_D,
                                 mn_ahead,
                                 tau11,
                                 w12t, cs_T,
                                 r12t, cs_R,
                                 c1, rs_C,
                                 C2, rs_C, cs_C,
                                 d1, rs_D,
                                 D2, rs_D, cs_D );

    //------------------------------------------------------------//

  }

  // FLA_Set_to_identity( T );
  bli_sident( mn_RT, buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_ONE, C, FLA_ONE, T );
  bli_ssyrk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_C,
             buff_1,
             buff_C, rs_C, cs_C,
             buff_1,
             buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_MINUS_ONE, D, FLA_ONE, T );
  bli_ssyrk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_D,
             buff_m1,
             buff_D, rs_D, cs_D,
             buff_1,
             buff_T, rs_T, cs_T );
 
  // FLA_Scale_diag( FLA_NO_CONJUGATE, FLA_ONE_HALF, T );
  bli_sscalediag( BLIS_NO_CONJUGATE,
                  0,
                  mn_RT,
                  mn_RT,
                  buff_half,
                  buff_T, rs_T, cs_T );

  return FLA_SUCCESS;
}



FLA_Error FLA_UDdate_UT_opd_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  double* buff_R, int rs_R, int cs_R,
                                  double* buff_C, int rs_C, int cs_C,
                                  double* buff_D, int rs_D, int cs_D,
                                  double* buff_T, int rs_T, int cs_T )
{
  double*   buff_half = FLA_DOUBLE_PTR( FLA_ONE_HALF );
  double*   buff_1    = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1   = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_RT; ++i )
  {
    double*   rho11     = buff_R + (i  )*cs_R + (i  )*rs_R;
    double*   r12t      = buff_R + (i+1)*cs_R + (i  )*rs_R;

    double*   c1        = buff_C + (i  )*cs_C + (0  )*rs_C;
    double*   C2        = buff_C + (i+1)*cs_C + (0  )*rs_C;

    double*   d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    double*   D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    double*   tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    double*   w12t      = buff_T + (i+1)*cs_T + (i  )*rs_T;

    int       mn_ahead  = mn_RT - i - 1;

    //------------------------------------------------------------//

    // FLA_Househ3UD_UT( rho11,
    //                   c1,
    //                   d1, tau11 );
    FLA_Househ3UD_UT_opd( m_C,
                          m_D,
                          rho11,
                          c1, rs_C,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_HUD_UT( FLA_LEFT,
    //                   tau11, r12t,
    //                   c1,    C2,
    //                   d1,    D2 );
    FLA_Apply_HUD_UT_l_opd_var1( m_C,
                                 m_D,
                                 mn_ahead,
                                 tau11,
                                 w12t, cs_T,
                                 r12t, cs_R,
                                 c1, rs_C,
                                 C2, rs_C, cs_C,
                                 d1, rs_D,
                                 D2, rs_D, cs_D );

    //------------------------------------------------------------//

  }

  // FLA_Set_to_identity( T );
  bli_dident( mn_RT, buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_ONE, C, FLA_ONE, T );
  bli_dsyrk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_C,
             buff_1,
             buff_C, rs_C, cs_C,
             buff_1,
             buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_MINUS_ONE, D, FLA_ONE, T );
  bli_dsyrk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_D,
             buff_m1,
             buff_D, rs_D, cs_D,
             buff_1,
             buff_T, rs_T, cs_T );
 
  // FLA_Scale_diag( FLA_NO_CONJUGATE, FLA_ONE_HALF, T );
  bli_dscalediag( BLIS_NO_CONJUGATE,
                  0,
                  mn_RT,
                  mn_RT,
                  buff_half,
                  buff_T, rs_T, cs_T );

  return FLA_SUCCESS;
}



FLA_Error FLA_UDdate_UT_opc_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  scomplex* buff_R, int rs_R, int cs_R,
                                  scomplex* buff_C, int rs_C, int cs_C,
                                  scomplex* buff_D, int rs_D, int cs_D,
                                  scomplex* buff_T, int rs_T, int cs_T )
{
  scomplex* buff_half = FLA_COMPLEX_PTR( FLA_ONE_HALF );
  float*    buff_1    = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1   = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_RT; ++i )
  {
    scomplex* rho11     = buff_R + (i  )*cs_R + (i  )*rs_R;
    scomplex* r12t      = buff_R + (i+1)*cs_R + (i  )*rs_R;

    scomplex* c1        = buff_C + (i  )*cs_C + (0  )*rs_C;
    scomplex* C2        = buff_C + (i+1)*cs_C + (0  )*rs_C;

    scomplex* d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    scomplex* D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    scomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    scomplex* w12t      = buff_T + (i+1)*cs_T + (i  )*rs_T;

    int       mn_ahead  = mn_RT - i - 1;

    //------------------------------------------------------------//

    // FLA_Househ3UD_UT( rho11,
    //                   c1,
    //                   d1, tau11 );
    FLA_Househ3UD_UT_opc( m_C,
                          m_D,
                          rho11,
                          c1, rs_C,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_HUD_UT( FLA_LEFT,
    //                   tau11, r12t,
    //                   c1,    C2,
    //                   d1,    D2 );
    FLA_Apply_HUD_UT_l_opc_var1( m_C,
                                 m_D,
                                 mn_ahead,
                                 tau11,
                                 w12t, cs_T,
                                 r12t, cs_R,
                                 c1, rs_C,
                                 C2, rs_C, cs_C,
                                 d1, rs_D,
                                 D2, rs_D, cs_D );

    //------------------------------------------------------------//

  }

  // FLA_Set_to_identity( T );
  bli_cident( mn_RT, buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_ONE, C, FLA_ONE, T );
  bli_cherk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_C,
             buff_1,
             buff_C, rs_C, cs_C,
             buff_1,
             buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_MINUS_ONE, D, FLA_ONE, T );
  bli_cherk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_D,
             buff_m1,
             buff_D, rs_D, cs_D,
             buff_1,
             buff_T, rs_T, cs_T );
 
  // FLA_Scale_diag( FLA_NO_CONJUGATE, FLA_ONE_HALF, T );
  bli_cscalediag( BLIS_NO_CONJUGATE,
                  0,
                  mn_RT,
                  mn_RT,
                  buff_half,
                  buff_T, rs_T, cs_T );

  return FLA_SUCCESS;
}



FLA_Error FLA_UDdate_UT_opz_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  dcomplex* buff_R, int rs_R, int cs_R,
                                  dcomplex* buff_C, int rs_C, int cs_C,
                                  dcomplex* buff_D, int rs_D, int cs_D,
                                  dcomplex* buff_T, int rs_T, int cs_T )
{
  dcomplex* buff_half = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  double*   buff_1    = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1   = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       i;

  for ( i = 0; i < mn_RT; ++i )
  {
    dcomplex* rho11     = buff_R + (i  )*cs_R + (i  )*rs_R;
    dcomplex* r12t      = buff_R + (i+1)*cs_R + (i  )*rs_R;

    dcomplex* c1        = buff_C + (i  )*cs_C + (0  )*rs_C;
    dcomplex* C2        = buff_C + (i+1)*cs_C + (0  )*rs_C;

    dcomplex* d1        = buff_D + (i  )*cs_D + (0  )*rs_D;
    dcomplex* D2        = buff_D + (i+1)*cs_D + (0  )*rs_D;

    dcomplex* tau11     = buff_T + (i  )*cs_T + (i  )*rs_T;
    dcomplex* w12t      = buff_T + (i+1)*cs_T + (i  )*rs_T;

    int       mn_ahead  = mn_RT - i - 1;

    //------------------------------------------------------------//

    // FLA_Househ3UD_UT( rho11,
    //                   c1,
    //                   d1, tau11 );
    FLA_Househ3UD_UT_opz( m_C,
                          m_D,
                          rho11,
                          c1, rs_C,
                          d1, rs_D,
                          tau11 );

    // FLA_Apply_HUD_UT( FLA_LEFT,
    //                   tau11, r12t,
    //                   c1,    C2,
    //                   d1,    D2 );
    FLA_Apply_HUD_UT_l_opz_var1( m_C,
                                 m_D,
                                 mn_ahead,
                                 tau11,
                                 w12t, cs_T,
                                 r12t, cs_R,
                                 c1, rs_C,
                                 C2, rs_C, cs_C,
                                 d1, rs_D,
                                 D2, rs_D, cs_D );

    //------------------------------------------------------------//

  }

  // FLA_Set_to_identity( T );
  bli_zident( mn_RT, buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_ONE, C, FLA_ONE, T );
  bli_zherk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_C,
             buff_1,
             buff_C, rs_C, cs_C,
             buff_1,
             buff_T, rs_T, cs_T );

  // FLA_Herk( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
  //           FLA_MINUS_ONE, D, FLA_ONE, T );
  bli_zherk( BLIS_UPPER_TRIANGULAR,
             BLIS_CONJ_TRANSPOSE,
             mn_RT,
             m_D,
             buff_m1,
             buff_D, rs_D, cs_D,
             buff_1,
             buff_T, rs_T, cs_T );
 
  // FLA_Scale_diag( FLA_NO_CONJUGATE, FLA_ONE_HALF, T );
  bli_zscalediag( BLIS_NO_CONJUGATE,
                  0,
                  mn_RT,
                  mn_RT,
                  buff_half,
                  buff_T, rs_T, cs_T );

  return FLA_SUCCESS;
}

