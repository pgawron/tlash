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

FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj y, FLA_Obj z )
{
/*
   Effective computation:
   y = y + delta * ( Y ( U' u ) + U ( Z' u ) );
   z = z + delta * ( U ( Y' u ) + Z ( U' u ) );
   t = U' u;
*/
  FLA_Datatype datatype;
  int          m_U, n_U;
  int          rs_U, cs_U;
  int          rs_Y, cs_Y;
  int          rs_Z, cs_Z;
  int          inc_u, inc_y, inc_z, inc_t;

  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  n_U      = FLA_Obj_width( U );

  rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  rs_Y     = FLA_Obj_row_stride( Y );
  cs_Y     = FLA_Obj_col_stride( Y );

  rs_Z     = FLA_Obj_row_stride( Z );
  cs_Z     = FLA_Obj_col_stride( Z );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_y    = FLA_Obj_vector_inc( y );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_t    = FLA_Obj_vector_inc( t );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_U     = FLA_FLOAT_PTR( U );
      float*    buff_Y     = FLA_FLOAT_PTR( Y );
      float*    buff_Z     = FLA_FLOAT_PTR( Z );
      float*    buff_t     = FLA_FLOAT_PTR( t );
      float*    buff_u     = FLA_FLOAT_PTR( u );
      float*    buff_y     = FLA_FLOAT_PTR( y );
      float*    buff_z     = FLA_FLOAT_PTR( z );
      float*    buff_delta = FLA_FLOAT_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_ops_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_U     = FLA_DOUBLE_PTR( U );
      double*   buff_Y     = FLA_DOUBLE_PTR( Y );
      double*   buff_Z     = FLA_DOUBLE_PTR( Z );
      double*   buff_t     = FLA_DOUBLE_PTR( t );
      double*   buff_u     = FLA_DOUBLE_PTR( u );
      double*   buff_y     = FLA_DOUBLE_PTR( y );
      double*   buff_z     = FLA_DOUBLE_PTR( z );
      double*   buff_delta = FLA_DOUBLE_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opd_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_U     = FLA_COMPLEX_PTR( U );
      scomplex* buff_Y     = FLA_COMPLEX_PTR( Y );
      scomplex* buff_Z     = FLA_COMPLEX_PTR( Z );
      scomplex* buff_t     = FLA_COMPLEX_PTR( t );
      scomplex* buff_u     = FLA_COMPLEX_PTR( u );
      scomplex* buff_y     = FLA_COMPLEX_PTR( y );
      scomplex* buff_z     = FLA_COMPLEX_PTR( z );
      scomplex* buff_delta = FLA_COMPLEX_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opc_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_U     = FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_Y     = FLA_DOUBLE_COMPLEX_PTR( Y );
      dcomplex* buff_Z     = FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_t     = FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex* buff_u     = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_y     = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_z     = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_delta = FLA_DOUBLE_COMPLEX_PTR( delta );

      FLA_Fused_Uhu_Yhu_Zhu_opz_var1( m_U,
                                      n_U,
                                      buff_delta,
                                      buff_U, rs_U, cs_U,
                                      buff_Y, rs_Y, cs_Y,
                                      buff_Z, rs_Z, cs_Z,
                                      buff_t, inc_t,
                                      buff_u, inc_u,
                                      buff_y, inc_y,
                                      buff_z, inc_z );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_ops_var1( int m_U,
                                          int n_U,
                                          float* buff_delta, 
                                          float* buff_U, int rs_U, int cs_U, 
                                          float* buff_Y, int rs_Y, int cs_Y, 
                                          float* buff_Z, int rs_Z, int cs_Z, 
                                          float* buff_t, int inc_t, 
                                          float* buff_u, int inc_u, 
                                          float* buff_y, int inc_y, 
                                          float* buff_z, int inc_z ) 
{
  int       i;

  for ( i = 0; i < n_U; ++i )
  {
    float*    u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    float*    y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    float*    z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    float*    delta    = buff_delta;
    float*    tau1     = buff_t + (i  )*inc_t;
    float*    u        = buff_u;
    float*    y        = buff_y;
    float*    z        = buff_z;
    float     alpha;
    float     beta;
    float     gamma;

    /*------------------------------------------------------------*/

    // bli_sdot( BLIS_CONJUGATE,
    //           m_U,
    //           u1, rs_U,
    //           u,  inc_u,
    //           &alpha );
    alpha = F77_sdot( &m_U,
                      u1, &rs_U,
                      u,  &inc_u );

    // bli_sdot( BLIS_CONJUGATE,
    //           m_U,
    //           z1, rs_Z,
    //           u,  inc_u,
    //           &beta );
    beta  = F77_sdot( &m_U,
                      z1, &rs_Z,
                      u,  &inc_u );

    // bli_sdot( BLIS_CONJUGATE,
    //           m_U,
    //           y1, rs_Y,
    //           u,  inc_u,
    //           &gamma );
    gamma = F77_sdot( &m_U,
                      y1, &rs_Y,
                      u,  &inc_u );

    *tau1 = alpha;

    // bli_sscals( delta, &alpha );
    // bli_sscals( delta, &beta );
    // bli_sscals( delta, &gamma );
    alpha *= *delta;
    beta  *= *delta;
    gamma *= *delta;

    // bli_saxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             y1, rs_Y,
    //             y,  inc_y );
    F77_saxpy( &m_U,
               &alpha,
               y1, &rs_Y,
               y,  &inc_y );

    // bli_saxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &beta,
    //             u1, rs_U,
    //             y,  inc_y );
    F77_saxpy( &m_U,
               &beta,
               u1, &rs_U,
               y,  &inc_y );

    // bli_saxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             z1, rs_Z,
    //             z,  inc_z );
    F77_saxpy( &m_U,
               &alpha,
               z1, &rs_Z,
               z,  &inc_z );

    // bli_saxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &gamma,
    //             u1, rs_U,
    //             z,  inc_z );
    F77_saxpy( &m_U,
               &gamma,
               u1, &rs_U,
               z,  &inc_z );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opd_var1( int m_U,
                                          int n_U,
                                          double* buff_delta, 
                                          double* buff_U, int rs_U, int cs_U, 
                                          double* buff_Y, int rs_Y, int cs_Y, 
                                          double* buff_Z, int rs_Z, int cs_Z, 
                                          double* buff_t, int inc_t, 
                                          double* buff_u, int inc_u, 
                                          double* buff_y, int inc_y, 
                                          double* buff_z, int inc_z ) 
{
  double             zero  = bli_d0();

  double*   restrict delta = buff_delta;
  double*   restrict u     = buff_u;
  double*   restrict y     = buff_y;
  double*   restrict z     = buff_z;

  double*   restrict u1;
  double*   restrict y1;
  double*   restrict z1;
  double*   restrict upsilon1;
  double*   restrict tau1;

  double    alpha;
  double    beta;
  double    gamma;

  int       i;

  int       n_run         = n_U / 1;
  //int       n_left        = n_U % 1;
  int       step_u1       = 1*cs_U;
  int       step_y1       = 1*cs_Y;
  int       step_z1       = 1*cs_Z;
  int       step_upsilon1 = 1*inc_u;
  int       step_tau1     = 1*inc_t;

  u1       = buff_U;
  y1       = buff_Y;
  z1       = buff_Z;
  upsilon1 = buff_u;
  tau1     = buff_t;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/

/*
    bli_ddotsv3( BLIS_CONJUGATE,
                 m_U,
                 u1, rs_U,
                 z1, rs_Z,
                 y1, rs_Y,
                 u,  inc_u,
                 &zero,
                 &alpha,
                 &beta,
                 &gamma );

    *tau1 = alpha;

    bli_dscals( delta, &alpha );
    bli_dscals( delta, &beta );
    bli_dscals( delta, &gamma );

    bli_daxpyv2b( m_U,
                  &alpha,
                  &beta,
                  y1, rs_Y,
                  u1, rs_U,
                  y, inc_y );
    bli_daxpyv2b( m_U,
                  &alpha,
                  &gamma,
                  z1, rs_Z,
                  u1, rs_U,
                  z, inc_z );
*/

    bli_ddotsv2( BLIS_CONJUGATE,
                 m_U,
                 y1, rs_Y,
                 z1, rs_Z,
                 u,  inc_u,
                 &zero,
                 &beta,
                 &gamma );

    bli_ddotaxmyv2( m_U,
                    &gamma,
                    &beta,
                    u1, rs_U,
                    u,  inc_u,
                    &alpha,
                    y,  inc_y,
                    z,  inc_z );

    *tau1 = alpha;

    bli_dscals( delta, &alpha );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    bli_daxpyv( BLIS_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );


    /*------------------------------------------------------------*/

    u1       += step_u1;
    y1       += step_y1;
    z1       += step_z1;
    upsilon1 += step_upsilon1;
    tau1     += step_tau1;
  }


  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opc_var1( int m_U,
                                          int n_U,
                                          scomplex* buff_delta, 
                                          scomplex* buff_U, int rs_U, int cs_U, 
                                          scomplex* buff_Y, int rs_Y, int cs_Y, 
                                          scomplex* buff_Z, int rs_Z, int cs_Z, 
                                          scomplex* buff_t, int inc_t, 
                                          scomplex* buff_u, int inc_u, 
                                          scomplex* buff_y, int inc_y, 
                                          scomplex* buff_z, int inc_z ) 
{
  int       i;

  for ( i = 0; i < n_U; ++i )
  {
    scomplex* u1       = buff_U + (i  )*cs_U + (0  )*rs_U;
    scomplex* y1       = buff_Y + (i  )*cs_Y + (0  )*rs_Y;
    scomplex* z1       = buff_Z + (i  )*cs_Z + (0  )*rs_Z;
    scomplex* delta    = buff_delta;
    scomplex* tau1     = buff_t + (i  )*inc_t;
    scomplex* u        = buff_u;
    scomplex* y        = buff_y;
    scomplex* z        = buff_z;
    scomplex  alpha;
    scomplex  beta;
    scomplex  gamma;

    /*------------------------------------------------------------*/

    bli_cdot( BLIS_CONJUGATE,
              m_U,
              u1, rs_U,
              u,  inc_u,
              &alpha );

    bli_cdot( BLIS_CONJUGATE,
              m_U,
              z1, rs_Z,
              u,  inc_u,
              &beta );

    bli_cdot( BLIS_CONJUGATE,
              m_U,
              y1, rs_Y,
              u,  inc_u,
              &gamma );

    *tau1 = alpha;

    bli_cscals( delta, &alpha );
    bli_cscals( delta, &beta );
    bli_cscals( delta, &gamma );

    // bli_caxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             y1, rs_Y,
    //             y,  inc_y );
    F77_caxpy( &m_U,
               &alpha,
               y1, &rs_Y,
               y,  &inc_y );

    // bli_caxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &beta,
    //             u1, rs_U,
    //             y,  inc_y );
    F77_caxpy( &m_U,
               &beta,
               u1, &rs_U,
               y,  &inc_y );

    // bli_caxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &alpha,
    //             z1, rs_Z,
    //             z,  inc_z );
    F77_caxpy( &m_U,
               &alpha,
               z1, &rs_Z,
               z,  &inc_z );

    // bli_caxpyv( BLIS_NO_CONJUGATE,
    //             m_U,
    //             &gamma,
    //             u1, rs_U,
    //             z,  inc_z );
    F77_caxpy( &m_U,
               &gamma,
               u1, &rs_U,
               z,  &inc_z );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opz_var1( int m_U,
                                          int n_U,
                                          dcomplex* buff_delta, 
                                          dcomplex* buff_U, int rs_U, int cs_U, 
                                          dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                          dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                          dcomplex* buff_t, int inc_t, 
                                          dcomplex* buff_u, int inc_u, 
                                          dcomplex* buff_y, int inc_y, 
                                          dcomplex* buff_z, int inc_z ) 
{
  dcomplex           zero  = bli_z0();

  dcomplex* restrict delta = buff_delta;
  dcomplex* restrict u     = buff_u;
  dcomplex* restrict y     = buff_y;
  dcomplex* restrict z     = buff_z;

  dcomplex* restrict u1;
  dcomplex* restrict y1;
  dcomplex* restrict z1;
  dcomplex* restrict upsilon1;
  dcomplex* restrict tau1;

  dcomplex  alpha;
  dcomplex  beta;
  dcomplex  gamma;

  int       i;

  int       n_run         = n_U / 1;
  //int       n_left        = n_U % 1;
  int       step_u1       = 1*cs_U;
  int       step_y1       = 1*cs_Y;
  int       step_z1       = 1*cs_Z;
  int       step_upsilon1 = 1*inc_u;
  int       step_tau1     = 1*inc_t;

  u1       = buff_U;
  y1       = buff_Y;
  z1       = buff_Z;
  upsilon1 = buff_u;
  tau1     = buff_t;

  for ( i = 0; i < n_run; ++i )
  {
    /*------------------------------------------------------------*/


    bli_zdotsv3( BLIS_CONJUGATE,
                 m_U,
                 u1, rs_U,
                 z1, rs_Z,
                 y1, rs_Y,
                 u,  inc_u,
                 &zero,
                 &alpha,
                 &beta,
                 &gamma );

    *tau1 = alpha;

    bli_zscals( delta, &alpha );
    bli_zscals( delta, &beta );
    bli_zscals( delta, &gamma );

    bli_zaxpyv2b( m_U,
                  &alpha,
                  &beta,
                  y1, rs_Y,
                  u1, rs_U,
                  y, inc_y );
    bli_zaxpyv2b( m_U,
                  &alpha,
                  &gamma,
                  z1, rs_Z,
                  u1, rs_U,
                  z, inc_z );


/*
    bli_zdotsv2( BLIS_CONJUGATE,
                 m_U,
                 y1, rs_Y,
                 z1, rs_Z,
                 u,  inc_u,
                 &zero,
                 &beta,
                 &gamma );

    bli_zdotaxmyv2( m_U,
                    &gamma,
                    &beta,
                    u1, rs_U,
                    u,  inc_u,
                    &alpha,
                    y,  inc_y,
                    z,  inc_z );

    *tau1 = alpha;

    bli_zscals( delta, &alpha );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_U,
                &alpha,
                y1, rs_Y,
                y,  inc_y );
    bli_zaxpyv( BLIS_NO_CONJUGATE,
                m_U,
                &alpha,
                z1, rs_Z,
                z,  inc_z );
*/

    /*------------------------------------------------------------*/

    u1       += step_u1;
    y1       += step_y1;
    z1       += step_z1;
    upsilon1 += step_upsilon1;
    tau1     += step_tau1;
  }

  return FLA_SUCCESS;
}

