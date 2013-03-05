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

FLA_Error FLA_Trsv_external_gpu( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu ) 
{
  FLA_Datatype datatype;
  int          m_A;
  int          ldim_A;
  int          inc_x;
  char         blas_uplo;
  char         blas_trans;
  char         blas_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  ldim_A   = FLA_Obj_length( A );

  inc_x    = 1;

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );
  FLA_Param_map_flame_to_netlib_diag( diag, &blas_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    cublasStrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( float * ) A_gpu, ldim_A,
                 ( float * ) x_gpu, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    cublasDtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( double * ) A_gpu, ldim_A,
                 ( double * ) x_gpu, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    cublasCtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( cuComplex * ) A_gpu, ldim_A,
                 ( cuComplex * ) x_gpu, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    cublasZtrsv( blas_uplo,
                 blas_trans,
                 blas_diag,
                 m_A,
                 ( cuDoubleComplex * ) A_gpu, ldim_A,
                 ( cuDoubleComplex * ) x_gpu, inc_x );

    break;
  }

  }

  return FLA_SUCCESS;
}

#endif
