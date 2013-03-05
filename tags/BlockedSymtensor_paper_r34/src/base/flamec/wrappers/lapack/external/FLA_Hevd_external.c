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

FLA_Error FLA_Hevd_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj e )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          n_A, cs_A;
  int          lwork, lrwork;
  FLA_Obj      work, rwork;
  char         blas_jobz;
  char         blas_uplo;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Hevd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  // Allocate the rwork array up front since its size is not dependent on
  // internal block sizes.
  lrwork   = max( 1, 3 * n_A - 2 );
  FLA_Obj_create( dt_real, lrwork, 1, 0, 0, &rwork );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );
  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size based on an internal block size.
  lwork = -1;
  FLA_Obj_create( datatype, 1, 1, 0, 0, &work );

  for ( i = 0; i < 2; ++i )
  {
    if ( i == 1 )
    {
      // Grab the queried ideal workspace size from the work array, free the
      // work object, and then re-allocate the workspace with the ideal size.
      if      ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
        lwork = ( int ) *FLA_FLOAT_PTR( work );
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
        lwork = ( int ) *FLA_DOUBLE_PTR( work );

      FLA_Obj_free( &work );
      FLA_Obj_create( datatype, lwork, 1, 0, 0, &work );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float* buff_A     = ( float * ) FLA_FLOAT_PTR( A );
      float* buff_e     = ( float * ) FLA_FLOAT_PTR( e );
      float* buff_work  = ( float * ) FLA_FLOAT_PTR( work );
      float* buff_rwork = ( float * ) FLA_FLOAT_PTR( rwork );

      F77_ssyev( &blas_jobz,
                 &blas_uplo,
                 &n_A,
                 buff_A,     &cs_A,
                 buff_e,
                 buff_work,  &lwork,
                 buff_rwork,
                 &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
      double* buff_e     = ( double * ) FLA_DOUBLE_PTR( e );
      double* buff_work  = ( double * ) FLA_DOUBLE_PTR( work );
      double* buff_rwork = ( double * ) FLA_DOUBLE_PTR( rwork );
  
      F77_dsyev( &blas_jobz,
                 &blas_uplo,
                 &n_A,
                 buff_A,     &cs_A,
                 buff_e,
                 buff_work,  &lwork,
                 buff_rwork,
                 &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      scomplex* buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float*    buff_e     = ( float    * ) FLA_FLOAT_PTR( e );
      scomplex* buff_work  = ( scomplex * ) FLA_COMPLEX_PTR( work );
      float*    buff_rwork = ( float    * ) FLA_FLOAT_PTR( rwork );
  
      F77_cheev( &blas_jobz,
                 &blas_uplo,
                 &n_A,
                 buff_A,     &cs_A,
                 buff_e,
                 buff_work,  &lwork,
                 buff_rwork,
                 &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_e     = ( double*   ) FLA_DOUBLE_PTR( e );
      dcomplex* buff_work  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( work );
      double*   buff_rwork = ( double*   ) FLA_DOUBLE_PTR( rwork );
  
      F77_zheev( &blas_jobz,
                 &blas_uplo,
                 &n_A,
                 buff_A,     &cs_A,
                 buff_e,
                 buff_work,  &lwork,
                 buff_rwork,
                 &info );
  
      break;
    } 

    }
  }

  FLA_Obj_free( &work );
  FLA_Obj_free( &rwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

