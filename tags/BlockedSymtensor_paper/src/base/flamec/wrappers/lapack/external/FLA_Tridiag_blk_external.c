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

FLA_Error FLA_Tridiag_blk_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          n_A, cs_A;
  int          lwork;
  FLA_Obj      d, e, work_obj;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Tridiag_check( uplo, A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), n_A,     1, 0, 0, &d );
  FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), n_A - 1, 1, 0, 0, &e );

  lwork    = n_A * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float* buff_d    = ( float * ) FLA_FLOAT_PTR( d );
    float* buff_e    = ( float * ) FLA_FLOAT_PTR( e );
    float* buff_t    = ( float * ) FLA_FLOAT_PTR( t );
    float* buff_work = ( float * ) FLA_FLOAT_PTR( work_obj );

    F77_ssytrd( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double* buff_d    = ( double * ) FLA_DOUBLE_PTR( d );
    double* buff_e    = ( double * ) FLA_DOUBLE_PTR( e );
    double* buff_t    = ( double * ) FLA_DOUBLE_PTR( t );
    double* buff_work = ( double * ) FLA_DOUBLE_PTR( work_obj );

    F77_dsytrd( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex* buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    float*    buff_d    = ( float    * ) FLA_FLOAT_PTR( d );
    float*    buff_e    = ( float    * ) FLA_FLOAT_PTR( e );
    scomplex* buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );
    scomplex* buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work_obj );

    F77_chetrd( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double*   buff_d    = ( double   * ) FLA_DOUBLE_PTR( d );
    double*   buff_e    = ( double   * ) FLA_DOUBLE_PTR( e );
    dcomplex* buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
    dcomplex* buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

    F77_zhetrd( &blas_uplo,
                &n_A,
                buff_A, &cs_A,
                buff_d,
                buff_e,
                buff_t,
                buff_work,
                &lwork,
                &info );

    break;
  } 

  }

  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &work_obj );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

FLA_Error FLA_Tridiag_blk_ext( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t )
{
  return FLA_Tridiag_blk_external( uplo, A, t );
}
