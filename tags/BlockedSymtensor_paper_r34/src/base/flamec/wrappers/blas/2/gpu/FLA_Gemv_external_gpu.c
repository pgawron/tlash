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

#ifdef FLA_ENABLE_GPU

#include "cublas.h"

FLA_Error FLA_Gemv_external_gpu( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          ldim_A;
  int          inc_x;
  int          inc_y;
  char         blas_transa;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemv_check( transa, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  ldim_A   = FLA_Obj_length( A );

  inc_x    = 1;
  inc_y    = 1;

  FLA_Param_map_flame_to_netlib_trans( transa, &blas_transa );


  switch( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    cublasSgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) x_gpu, inc_x,
                 *buff_beta,
                 ( float * ) y_gpu, inc_y );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    cublasDgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) x_gpu, inc_x,
                 *buff_beta,
                 ( double * ) y_gpu, inc_y );

    break;
  }

  case FLA_COMPLEX:
  {
    cuComplex *buff_alpha = ( cuComplex * ) FLA_COMPLEX_PTR( alpha );
    cuComplex *buff_beta  = ( cuComplex * ) FLA_COMPLEX_PTR( beta );

    cublasCgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) x_gpu, inc_x,
                 *buff_beta,
                 ( cuComplex * ) y_gpu, inc_y );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cuDoubleComplex *buff_alpha = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    cuDoubleComplex *buff_beta  = ( cuDoubleComplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    cublasZgemv( blas_transa,
                 m_A,
                 n_A, 
                 *buff_alpha,                   
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) x_gpu, inc_x,
                 *buff_beta,
                 ( cuDoubleComplex * ) y_gpu, inc_y );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

#endif
