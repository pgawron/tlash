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

FLA_Error FLA_Apply_H2_UT_l_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 )
/*
  Compute:

    / a1t \  :=  / I - 1/tau / 1  \ ( 1  u2' ) \ / a1t \ 
    \ A2  /      \           \ u2 /            / \ A2  / 
 
  w = ( a1t + u2' * A2 ) / tau;

  a1t = a1t - w;
  A2  = A2  - u2 * w;
*/
{
  FLA_Datatype datatype;
  int          m_u2_A2;
  int          n_a1t;
  int          inc_u2;
  int          inc_a1t;
  int          rs_A2;
  int          cs_A2;

  if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A2 );

  m_u2_A2  = FLA_Obj_length( A2 );
  n_a1t    = FLA_Obj_width( a1t );
  inc_u2   = FLA_Obj_vector_inc( u2 );
  inc_a1t  = FLA_Obj_vector_inc( a1t );
  rs_A2    = FLA_Obj_row_stride( A2 );
  cs_A2    = FLA_Obj_col_stride( A2 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* tau_p = ( float* ) FLA_FLOAT_PTR( tau );
      float* u2_p  = ( float* ) FLA_FLOAT_PTR( u2 );
      float* a1t_p = ( float* ) FLA_FLOAT_PTR( a1t );
      float* A2_p  = ( float* ) FLA_FLOAT_PTR( A2 );

      FLA_Apply_H2_UT_l_ops_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE:
    {
      double* tau_p = ( double* ) FLA_DOUBLE_PTR( tau );
      double* u2_p  = ( double* ) FLA_DOUBLE_PTR( u2 );
      double* a1t_p = ( double* ) FLA_DOUBLE_PTR( a1t );
      double* A2_p  = ( double* ) FLA_DOUBLE_PTR( A2 );

      FLA_Apply_H2_UT_l_opd_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* tau_p = ( scomplex* ) FLA_COMPLEX_PTR( tau );
      scomplex* u2_p  = ( scomplex* ) FLA_COMPLEX_PTR( u2 );
      scomplex* a1t_p = ( scomplex* ) FLA_COMPLEX_PTR( a1t );
      scomplex* A2_p  = ( scomplex* ) FLA_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_l_opc_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* tau_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* u2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( u2 );
      dcomplex* a1t_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( a1t );
      dcomplex* A2_p  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A2 );

      FLA_Apply_H2_UT_l_opz_var1( m_u2_A2, n_a1t,
                                  tau_p,
                                  u2_p, inc_u2,
                                  a1t_p, inc_a1t,
                                  A2_p, rs_A2, cs_A2 );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_ops_var1( int m_u2_A2,
                                      int n_a1t,
                                      float* tau,
                                      float* u2, int inc_u2,
                                      float* a1t, int inc_a1t,
                                      float* A2, int rs_A2, int cs_A2 )
{
  float*    one_p       = FLA_FLOAT_PTR( FLA_ONE );
  float*    minus_one_p = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       inc_w1t;

  // FLA_Obj w1t;
  float*    w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( float* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bli_scopyv( BLIS_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bli_sgemv( BLIS_TRANSPOSE,
             BLIS_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bli_sinvscalv( BLIS_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bli_saxpyv( BLIS_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bli_sger( BLIS_NO_CONJUGATE,
            BLIS_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opd_var1( int m_u2_A2,
                                      int n_a1t,
                                      double* tau,
                                      double* u2, int inc_u2,
                                      double* a1t, int inc_a1t,
                                      double* A2, int rs_A2, int cs_A2 )
{
  double*   one_p       = FLA_DOUBLE_PTR( FLA_ONE );
  double*   minus_one_p = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       inc_w1t;

  // FLA_Obj w1t;
  double*   w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( double* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bli_dcopyv( BLIS_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bli_dgemv( BLIS_TRANSPOSE,
             BLIS_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bli_dinvscalv( BLIS_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bli_daxpyv( BLIS_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bli_dger( BLIS_NO_CONJUGATE,
            BLIS_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opc_var1( int m_u2_A2,
                                      int n_a1t,
                                      scomplex* tau,
                                      scomplex* u2, int inc_u2,
                                      scomplex* a1t, int inc_a1t,
                                      scomplex* A2, int rs_A2, int cs_A2 )
{
  scomplex* one_p       = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* minus_one_p = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1t;

  // FLA_Obj w1t;
  scomplex* w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( scomplex* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bli_ccopyv( BLIS_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bli_cgemv( BLIS_TRANSPOSE,
             BLIS_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bli_cinvscalv( BLIS_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bli_caxpyv( BLIS_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bli_cger( BLIS_NO_CONJUGATE,
            BLIS_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}



FLA_Error FLA_Apply_H2_UT_l_opz_var1( int m_u2_A2,
                                      int n_a1t,
                                      dcomplex* tau,
                                      dcomplex* u2, int inc_u2,
                                      dcomplex* a1t, int inc_a1t,
                                      dcomplex* A2, int rs_A2, int cs_A2 )
{
  dcomplex* one_p       = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* minus_one_p = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       inc_w1t;

  // FLA_Obj w1t;
  dcomplex* w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  w1t = ( dcomplex* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  bli_zcopyv( BLIS_NO_CONJUGATE,
              n_a1t,
              a1t, inc_a1t, 
              w1t, inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );
  bli_zgemv( BLIS_TRANSPOSE,
             BLIS_CONJUGATE,
             m_u2_A2,
             n_a1t,
             one_p,
             A2, rs_A2, cs_A2,
             u2, inc_u2,
             one_p,
             w1t, inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  bli_zinvscalv( BLIS_NO_CONJUGATE,
                 n_a1t,
                 tau,
                 w1t, inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  bli_zaxpyv( BLIS_NO_CONJUGATE,
              n_a1t,
              minus_one_p,
              w1t, inc_w1t,
              a1t, inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  bli_zger( BLIS_NO_CONJUGATE,
            BLIS_NO_CONJUGATE,
            m_u2_A2,
            n_a1t,
            minus_one_p,
            u2, inc_u2,
            w1t, inc_w1t,
            A2, rs_A2, cs_A2 );

  // FLA_Obj_free( &w1t );
  FLA_free( w1t );

  return FLA_SUCCESS;
}
