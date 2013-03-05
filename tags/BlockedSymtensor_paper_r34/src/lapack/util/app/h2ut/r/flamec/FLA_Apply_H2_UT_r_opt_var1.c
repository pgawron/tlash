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

FLA_Error FLA_Apply_H2_UT_r_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
  FLA_Datatype datatype;
  int          n_u2_A2;
  int          m_a1;
  int          inc_u2;
  int          inc_a1;
  int          rs_A2, cs_A2;

  if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A2 );

  n_u2_A2  = FLA_Obj_width( A2 );
  m_a1     = FLA_Obj_length( a1 );
  inc_u2   = FLA_Obj_vector_inc( u2 );
  inc_a1   = FLA_Obj_vector_inc( a1 );
  rs_A2    = FLA_Obj_row_stride( A2 );
  cs_A2    = FLA_Obj_col_stride( A2 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* tau_p = ( float* ) FLA_FLOAT_PTR( tau );
      float* u2_p  = ( float* ) FLA_FLOAT_PTR( u2 );
      float* a1_p  = ( float* ) FLA_FLOAT_PTR( a1 );
      float* A2_p  = ( float* ) FLA_FLOAT_PTR( A2 );

      FLA_Apply_H2_UT_r_ops_var1( m_a1, n_u2_A2,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1_p, inc_a1,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE:
    {
      double* tau_p = ( double* ) FLA_DOUBLE_PTR( tau );
      double* u2_p  = ( double* ) FLA_DOUBLE_PTR( u2 );
      double* a1_p  = ( double* ) FLA_DOUBLE_PTR( a1 );
      double* A2_p  = ( double* ) FLA_DOUBLE_PTR( A2 );

      FLA_Apply_H2_UT_r_opd_var1( m_a1, n_u2_A2,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1_p, inc_a1,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* tau_p = ( scomplex* ) FLA_COMPLEX_PTR( tau );
      scomplex* u2_p  = ( scomplex* ) FLA_COMPLEX_PTR( u2 );
      scomplex* a1_p  = ( scomplex* ) FLA_COMPLEX_PTR( a1 );
      scomplex* A2_p  = ( scomplex* ) FLA_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_r_opc_var1( m_a1, n_u2_A2,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1_p, inc_a1,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* tau_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* u2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( u2 );
      dcomplex* a1_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( a1 );
      dcomplex* A2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_r_opz_var1( m_a1, n_u2_A2,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1_p, inc_a1,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_ops_var1( int m_a1,
                                      int n_u2_A2,
                                      float* tau,
                                      float* u2, int inc_u2,
                                      float* a1, int inc_a1,
                                      float* A2, int rs_A2, int cs_A2 )
{
  float*    one_p       = FLA_FLOAT_PTR( FLA_ONE );
  float*    minus_one_p = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  float*    w1;

  // if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;
  if ( m_a1 == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( float* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bli_scopyv( BLIS_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bli_sgemv( BLIS_NO_TRANSPOSE,
             BLIS_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bli_sinvscalv( BLIS_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bli_saxpyv( BLIS_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bli_sger( BLIS_NO_CONJUGATE,
            BLIS_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opd_var1( int m_a1,
                                      int n_u2_A2,
                                      double* tau,
                                      double* u2, int inc_u2,
                                      double* a1, int inc_a1,
                                      double* A2, int rs_A2, int cs_A2 )
{
  double*   one_p       = FLA_DOUBLE_PTR( FLA_ONE );
  double*   minus_one_p = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  double*   w1;

  // if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;
  if ( m_a1 == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( double* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bli_dcopyv( BLIS_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bli_dgemv( BLIS_NO_TRANSPOSE,
             BLIS_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bli_dinvscalv( BLIS_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bli_daxpyv( BLIS_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bli_dger( BLIS_NO_CONJUGATE,
            BLIS_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opc_var1( int m_a1,
                                      int n_u2_A2,
                                      scomplex* tau,
                                      scomplex* u2, int inc_u2,
                                      scomplex* a1, int inc_a1,
                                      scomplex* A2, int rs_A2, int cs_A2 )
{
  scomplex* one_p       = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* minus_one_p = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  scomplex* w1;

  // if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;
  if ( m_a1 == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( scomplex* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bli_ccopyv( BLIS_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bli_cgemv( BLIS_NO_TRANSPOSE,
             BLIS_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bli_cinvscalv( BLIS_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bli_caxpyv( BLIS_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bli_cger( BLIS_NO_CONJUGATE,
            BLIS_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_r_opz_var1( int m_a1,
                                      int n_u2_A2,
                                      dcomplex* tau,
                                      dcomplex* u2, int inc_u2,
                                      dcomplex* a1, int inc_a1,
                                      dcomplex* A2, int rs_A2, int cs_A2 )
{
  dcomplex* one_p       = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* minus_one_p = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1;

  // FLA_Obj w1;
  dcomplex* w1;

  // if ( FLA_Obj_has_zero_dim( a1 ) ) return FLA_SUCCESS;
  if ( m_a1 == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1, &w1 );
  w1 = ( dcomplex* ) FLA_malloc( m_a1 * sizeof( *a1 ) );
  inc_w1 = 1;

  // // w1 = a1;
  // FLA_Copy_external( a1, w1 );
  bli_zcopyv( BLIS_NO_CONJUGATE,
              m_a1,
              a1, inc_a1, 
              w1, inc_w1 ); 

  // // w1 = w1 + A2 * u2;
  // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A2, u2, FLA_ONE, w1 );
  bli_zgemv( BLIS_NO_TRANSPOSE,
             BLIS_NO_CONJUGATE,
             m_a1,
             n_u2_A2,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1, inc_w1 );

  // // w1 = w1 / tau;
  // FLA_Inv_scal_external( tau, w1 );
  bli_zinvscalv( BLIS_NO_CONJUGATE,
                 m_a1,
                 tau,
                 w1, inc_w1 );

  // // a1 = a1 - w1;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1, a1 );
  bli_zaxpyv( BLIS_NO_CONJUGATE,
              m_a1,
              minus_one_p,
              w1, inc_w1,
              a1, inc_a1 );

  // // A2 = A2 - w1 * u2';
  // FLA_Gerc( FLA_NO_CONJUGATE, FLA_CONJUGATE, FLA_MINUS_ONE, w1, u2, A2 );
  bli_zger( BLIS_NO_CONJUGATE,
            BLIS_CONJUGATE,
            m_a1,
            n_u2_A2,
            minus_one_p,
            w1, inc_w1,
            u2, inc_u2,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1 );
  FLA_free( w1 );

  return FLA_SUCCESS;
}
