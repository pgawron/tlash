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

FLA_Error FLA_Bsvd_external( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          m_U, cs_U;
  int          n_V, cs_V;
  int          n_C, cs_C;
  int          min_m_n;
  int          inc_d, inc_e;
  int          lrwork;
  FLA_Obj      rwork;
  char         blas_uplo;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Hevd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( d ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( U );
  dt_real  = FLA_Obj_datatype_proj_to_real( U );

  m_U      = FLA_Obj_length( U );
  cs_U     = FLA_Obj_col_stride( U );

  n_V      = FLA_Obj_length( V );
  cs_V     = FLA_Obj_col_stride( V );

  n_C      = 0;
  cs_C     = 1;

  min_m_n  = FLA_Obj_vector_dim( d );

  inc_d    = FLA_Obj_vector_inc( d );
  inc_e    = FLA_Obj_vector_inc( e );

  lrwork   = max( 1, 4 * min_m_n - 4 );
  FLA_Obj_create( dt_real, lrwork, 1, 0, 0, &rwork );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d     = ( float * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float * ) FLA_FLOAT_PTR( e );
      float*    buff_U     = ( float * ) FLA_FLOAT_PTR( U );
      float*    buff_V     = ( float * ) FLA_FLOAT_PTR( V );
      float*    buff_C     = ( float * ) NULL;
      float*    buff_rwork = ( float * ) FLA_FLOAT_PTR( rwork );
  
      F77_sbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d     = ( double * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double * ) FLA_DOUBLE_PTR( e );
      double*   buff_U     = ( double * ) FLA_DOUBLE_PTR( U );
      double*   buff_V     = ( double * ) FLA_DOUBLE_PTR( V );
      double*   buff_C     = ( double * ) NULL;
      double*   buff_rwork = ( double * ) FLA_DOUBLE_PTR( rwork );
  
      F77_dbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d     = ( float    * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float    * ) FLA_FLOAT_PTR( e );
      scomplex* buff_U     = ( scomplex * ) FLA_COMPLEX_PTR( U );
      scomplex* buff_V     = ( scomplex * ) FLA_COMPLEX_PTR( V );
      scomplex* buff_C     = ( scomplex * ) NULL;
      float*    buff_rwork = ( float    * ) FLA_FLOAT_PTR( rwork );
  
      F77_cbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d     = ( double   * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double   * ) FLA_DOUBLE_PTR( e );
      dcomplex* buff_U     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_V     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( V );
      dcomplex* buff_C     = ( dcomplex * ) NULL;
      double*   buff_rwork = ( double   * ) FLA_DOUBLE_PTR( rwork );
  
      F77_zbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 

    }

  FLA_Obj_free( &rwork );

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

