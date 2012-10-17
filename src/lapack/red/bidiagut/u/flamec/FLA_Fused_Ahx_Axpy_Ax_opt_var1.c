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

FLA_Error FLA_Fused_Ahx_Axpy_Ax_opt_var1( FLA_Obj A, FLA_Obj u, FLA_Obj tau, FLA_Obj a, FLA_Obj beta, FLA_Obj y, FLA_Obj w )
{
/*
   Effective computation:
   y = beta * y + A' * u;
   a = a - conj(y) / tau;
   w = A * conj(a);
*/
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_u, inc_a, inc_y, inc_w;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_u    = FLA_Obj_vector_inc( u );

  inc_a    = FLA_Obj_vector_inc( a );

  inc_y    = FLA_Obj_vector_inc( y );

  inc_w    = FLA_Obj_vector_inc( w );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A   = FLA_FLOAT_PTR( A );
      float* buff_u   = FLA_FLOAT_PTR( u );
      float* buff_a   = FLA_FLOAT_PTR( a );
      float* buff_y   = FLA_FLOAT_PTR( y );
      float* buff_w   = FLA_FLOAT_PTR( w );
      float* buff_tau = FLA_FLOAT_PTR( tau );
      float* buff_beta = FLA_FLOAT_PTR( beta );

      FLA_Fused_Ahx_Axpy_Ax_ops_var1( m_A,
                                      n_A,
                                      buff_tau,
                                      buff_beta,
                                      buff_A, rs_A, cs_A,
                                      buff_u, inc_u,
                                      buff_a, inc_a,
                                      buff_y, inc_y,
                                      buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A   = FLA_DOUBLE_PTR( A );
      double* buff_u   = FLA_DOUBLE_PTR( u );
      double* buff_a   = FLA_DOUBLE_PTR( a );
      double* buff_y   = FLA_DOUBLE_PTR( y );
      double* buff_w   = FLA_DOUBLE_PTR( w );
      double* buff_tau = FLA_DOUBLE_PTR( tau );
      double* buff_beta = FLA_DOUBLE_PTR( beta );

      FLA_Fused_Ahx_Axpy_Ax_opd_var1( m_A,
                                      n_A,
                                      buff_tau,
                                      buff_beta,
                                      buff_A, rs_A, cs_A,
                                      buff_u, inc_u,
                                      buff_a, inc_a,
                                      buff_y, inc_y,
                                      buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A   = FLA_COMPLEX_PTR( A );
      scomplex* buff_u   = FLA_COMPLEX_PTR( u );
      scomplex* buff_a   = FLA_COMPLEX_PTR( a );
      scomplex* buff_y   = FLA_COMPLEX_PTR( y );
      scomplex* buff_w   = FLA_COMPLEX_PTR( w );
      scomplex* buff_tau = FLA_COMPLEX_PTR( tau );
      scomplex* buff_beta = FLA_COMPLEX_PTR( beta );

      FLA_Fused_Ahx_Axpy_Ax_opc_var1( m_A,
                                      n_A,
                                      buff_tau,
                                      buff_beta,
                                      buff_A, rs_A, cs_A,
                                      buff_u, inc_u,
                                      buff_a, inc_a,
                                      buff_y, inc_y,
                                      buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A   = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_u   = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_a   = FLA_DOUBLE_COMPLEX_PTR( a );
      dcomplex* buff_y   = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_w   = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_tau = FLA_DOUBLE_COMPLEX_PTR( tau );
      dcomplex* buff_beta = FLA_DOUBLE_COMPLEX_PTR( beta );

      FLA_Fused_Ahx_Axpy_Ax_opz_var1( m_A,
                                      n_A,
                                      buff_tau,
                                      buff_beta,
                                      buff_A, rs_A, cs_A,
                                      buff_u, inc_u,
                                      buff_a, inc_a,
                                      buff_y, inc_y,
                                      buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Axpy_Ax_ops_var1( int m_A,
                                          int n_A,
                                          float* buff_tau, 
                                          float* buff_beta, 
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_u, int inc_u, 
                                          float* buff_a, int inc_a, 
                                          float* buff_y, int inc_y, 
                                          float* buff_w, int inc_w )
{
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  float     minus_inv_tau;
  float     rho;
  int       i;

  bli_ssetv( m_A,
             buff_0,
             buff_w, inc_w );

  minus_inv_tau = *buff_m1 / *buff_tau;

  for ( i = 0; i < n_A; ++i )
  {
    float*    a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    psi1     = buff_y + (i  )*inc_y;
    float*    alpha1   = buff_a + (i  )*inc_a;
    float*    u        = buff_u;
    float*    w        = buff_w;

    /*------------------------------------------------------------*/

    *psi1 = *buff_beta * *psi1;
    // bli_sdot( BLIS_CONJUGATE,
    //           m_A,
    //           a1, rs_A,
    //           u,  inc_u,
    //           psi1 );
    rho = F77_sdot( &m_A,
                    a1, &rs_A,
                    u,  &inc_u );
    *psi1 = *psi1 + rho;

    // bli_dmult4( &minus_inv_tau, conj_psi1, alpha1, alpha1 );
    *alpha1 = *alpha1 + minus_inv_tau * *psi1;

    // bli_saxpyv( BLIS_NO_CONJUGATE,
    //             m_A,
    //             conj_alpha1,
    //             a1, rs_A,
    //             w,  inc_w );
    F77_saxpy( &m_A,
               alpha1,
               a1, &rs_A,
               w,  &inc_w );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Axpy_Ax_opd_var1( int m_A,
                                          int n_A,
                                          double* buff_tau, 
                                          double* buff_beta, 
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_u, int inc_u, 
                                          double* buff_a, int inc_a, 
                                          double* buff_y, int inc_y, 
                                          double* buff_w, int inc_w )
{
  double    zero      = bli_d0();
  double    minus_one = bli_dm1();
  double*   restrict u = buff_u;
  double*   restrict w = buff_w;
  double*   restrict beta = buff_beta;
  double*   restrict a1;
  double*   restrict a2;
  double*   restrict psi1;
  double*   restrict psi2;
  double*   restrict alpha1;
  double*   restrict alpha2;

  double    minus_inv_tau;
  int       i;

  int       n_run    = n_A / 2;
  int       n_left   = n_A % 2;
  int       stepcs_A  = 2*cs_A;
  int       stepinc_y = 2*inc_y;
  int       stepinc_a = 2*inc_a;


  bli_dsetv( m_A,
             &zero,
             buff_w, inc_w );

  bli_ddiv3( &minus_one, buff_tau, &minus_inv_tau );

  a1     = buff_A;
  a2     = buff_A + cs_A;
  psi1   = buff_y;
  psi2   = buff_y + inc_y;
  alpha1 = buff_a;
  alpha2 = buff_a + inc_a;

  for ( i = 0; i < n_run; ++i )
  {
/*
   Effective computation:
   y = beta * y + A' * u;
   a = a - conj(y) / tau;
   w = A * conj(a);
*/
    /*------------------------------------------------------------*/

    bli_ddotsv2( BLIS_CONJUGATE,
                 m_A,
                 a1, rs_A,
                 a2, rs_A,
                 u,  inc_u,
                 beta,
                 psi1,
                 psi2 );

    bli_dmult4( &minus_inv_tau, psi1, alpha1, alpha1 );
    bli_dmult4( &minus_inv_tau, psi2, alpha2, alpha2 );

    bli_daxpyv2b( m_A,
                  alpha1,
                  alpha2,
                  a1, rs_A,
                  a2, rs_A,
                  w,  inc_w );

    /*------------------------------------------------------------*/

    a1     += stepcs_A;
    a2     += stepcs_A;
    psi1   += stepinc_y;
    psi2   += stepinc_y;
    alpha1 += stepinc_a;
    alpha2 += stepinc_a;
  }

  if ( n_left == 1 )
  //for ( i = 0; i < n_left; ++i )
  {
    double   rho1;

    bli_ddot( BLIS_CONJUGATE,
              m_A,
              a1, rs_A,
              u,  inc_u,
              &rho1 );
    bli_dscals( buff_beta, psi1 );
    bli_dadd3( psi1, &rho1, psi1 );

    bli_dmult4( &minus_inv_tau, psi1, alpha1, alpha1 );

    bli_daxpyv( BLIS_NO_CONJUGATE,
	            m_A,
                alpha1,
                a1, rs_A,
                w,  inc_w );

    //a1     += cs_A;
    //psi1   += inc_y;
    //alpha1 += inc_a;
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Axpy_Ax_opc_var1( int m_A,
                                          int n_A,
                                          scomplex* buff_tau, 
                                          scomplex* buff_beta, 
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_u, int inc_u, 
                                          scomplex* buff_a, int inc_a, 
                                          scomplex* buff_y, int inc_y, 
                                          scomplex* buff_w, int inc_w )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  scomplex  minus_inv_tau;
  scomplex  conj_psi1;
  scomplex  conj_alpha1;
  int       i;

  bli_csetv( m_A,
             buff_0,
             buff_w, inc_w );

  bli_cdiv3( buff_m1, buff_tau, &minus_inv_tau );

  for ( i = 0; i < n_A; ++i )
  {
    scomplex* a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* psi1     = buff_y + (i  )*inc_y;
    scomplex* alpha1   = buff_a + (i  )*inc_a;
    scomplex* u        = buff_u;
    scomplex* w        = buff_w;

    /*------------------------------------------------------------*/

    bli_cdots( BLIS_CONJUGATE,
               m_A,
               buff_1,
               a1, rs_A,
               u,  inc_u,
               buff_beta,
               psi1 );

    bli_ccopyconj( psi1, &conj_psi1 );
    bli_cmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );

    bli_ccopyconj( alpha1, &conj_alpha1 );

    // bli_caxpyv( BLIS_NO_CONJUGATE,
    //             m_A,
    //             conj_alpha1,
    //             a1, rs_A,
    //             w,  inc_w );
    F77_caxpy( &m_A,
               &conj_alpha1,
               a1, &rs_A,
               w,  &inc_w );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Ahx_Axpy_Ax_opz_var1( int m_A,
                                          int n_A,
                                          dcomplex* buff_tau, 
                                          dcomplex* buff_beta, 
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_u, int inc_u, 
                                          dcomplex* buff_a, int inc_a, 
                                          dcomplex* buff_y, int inc_y, 
                                          dcomplex* buff_w, int inc_w )
{
  dcomplex  zero      = bli_z0();
  dcomplex  minus_one = bli_zm1();
  dcomplex* restrict u = buff_u;
  dcomplex* restrict w = buff_w;
  dcomplex* restrict beta = buff_beta;
  dcomplex* restrict a1;
  dcomplex* restrict a2;
  dcomplex* restrict psi1;
  dcomplex* restrict psi2;
  dcomplex* restrict alpha1;
  dcomplex* restrict alpha2;

  dcomplex  minus_inv_tau;
  dcomplex  conj_psi1;
  dcomplex  conj_psi2;
  dcomplex  conj_alpha1;
  dcomplex  conj_alpha2;
  int       i;
  int       n_run    = n_A / 2;
  int       n_left   = n_A % 2;
  int       twocs_A  = 2*cs_A;
  int       twoinc_y = 2*inc_y;
  int       twoinc_a = 2*inc_a;


  bli_zsetv( m_A,
             &zero,
             buff_w, inc_w );

  bli_zdiv3( &minus_one, buff_tau, &minus_inv_tau );

  a1     = buff_A;
  a2     = buff_A + cs_A;
  psi1   = buff_y;
  psi2   = buff_y + inc_y;
  alpha1 = buff_a;
  alpha2 = buff_a + inc_a;

  for ( i = 0; i < n_run; ++i )
  {
/*
   Effective computation:
   y = beta * y + A' * u;
   a = a - conj(y) / tau;
   w = A * conj(a);
*/
    /*------------------------------------------------------------*/

    bli_zdotsv2( BLIS_CONJUGATE,
                 m_A,
                 a1, rs_A,
                 a2, rs_A,
                 u,  inc_u,
                 beta,
                 psi1,
                 psi2 );

    bli_zcopyconj( psi1, &conj_psi1 );
    bli_zcopyconj( psi2, &conj_psi2 );
    bli_zmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );
    bli_zmult4( &minus_inv_tau, &conj_psi2, alpha2, alpha2 );
    bli_zcopyconj( alpha1, &conj_alpha1 );
    bli_zcopyconj( alpha2, &conj_alpha2 );

    bli_zaxpyv2b( m_A,
                  &conj_alpha1,
                  &conj_alpha2,
                  a1, rs_A,
                  a2, rs_A,
                  w,  inc_w );

    /*------------------------------------------------------------*/

    a1     += twocs_A;
    a2     += twocs_A;
    psi1   += twoinc_y;
    psi2   += twoinc_y;
    alpha1 += twoinc_a;
    alpha2 += twoinc_a;
  }

  if ( n_left == 1 )
  {
    dcomplex rho1;

    bli_zdot( BLIS_CONJUGATE,
              m_A,
              a1, rs_A,
              u,  inc_u,
              &rho1 );
    bli_zscals( buff_beta, psi1 );
    bli_zadd3( psi1, &rho1, psi1 );

    bli_zcopyconj( psi1, &conj_psi1 );
    bli_zmult4( &minus_inv_tau, &conj_psi1, alpha1, alpha1 );
    bli_zcopyconj( alpha1, &conj_alpha1 );

    bli_zaxpyv( BLIS_NO_CONJUGATE,
	            m_A,
                &conj_alpha1,
                a1, rs_A,
                w,  inc_w );
  }

  return FLA_SUCCESS;
}

