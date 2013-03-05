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

FLA_Error FLA_Tevdr_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj A )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          n_A, cs_A;
  int          lisuppz, lwork, liwork;
  FLA_Obj      isuppz, work, iwork;
  char         blas_jobz;
  char         blas_range;
  int          i;
  int          vl, vu;
  int          il, iu;
  int          nzc;
  int          try_rac;
  int          n_eig_found;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Tevdd_check( jobz, d, e, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );

  // Hard-code some parameters.
  blas_range = 'A';
  nzc        = n_A;
  try_rac    = TRUE;

  // Allocate space for the isuppz array.
  lisuppz = 2 * n_A;
  FLA_Obj_create( FLA_INT, lisuppz, 1, 0, 0, &isuppz );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size.
  lwork = -1;
  liwork = -1;
  FLA_Obj_create( dt_real, 1, 1, 0, 0, &work );
  FLA_Obj_create( FLA_INT, 1, 1, 0, 0, &iwork );

  for ( i = 0; i < 2; ++i )
  {
    if ( i == 1 )
    {
      // Grab the queried ideal workspace size from the work arrays, free the
      // work object, and then re-allocate the workspace with the ideal size.
      if      ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
      {
        lwork  = ( int ) *FLA_FLOAT_PTR( work );
        liwork = ( int ) *FLA_INT_PTR( iwork );
      }
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
      {
        lwork  = ( int ) *FLA_DOUBLE_PTR( work );
        liwork = ( int ) *FLA_INT_PTR( iwork );
      }
//printf( "ideal workspace for n = %d\n", n_A );
//printf( "                lwork = %d\n", lwork );
//printf( "               liwork = %d\n", liwork );
      FLA_Obj_free( &work );
      FLA_Obj_free( &iwork );
      FLA_Obj_create( dt_real, lwork,  1, 0, 0, &work );
      FLA_Obj_create( FLA_INT, liwork, 1, 0, 0, &iwork );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d      = ( float*    ) FLA_FLOAT_PTR( d );
      float*    buff_e      = ( float*    ) FLA_FLOAT_PTR( e );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      float*    buff_A      = ( float*    ) FLA_FLOAT_PTR( A );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      float*    buff_work   = ( float*    ) FLA_FLOAT_PTR( work );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );
 
      F77_sstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vl, &vu,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d      = ( double*   ) FLA_DOUBLE_PTR( d );
      double*   buff_e      = ( double*   ) FLA_DOUBLE_PTR( e );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      double*   buff_A      = ( double*   ) FLA_DOUBLE_PTR( A );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      double*   buff_work   = ( double*   ) FLA_DOUBLE_PTR( work );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );

      F77_dstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vl, &vu,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d      = ( float*    ) FLA_FLOAT_PTR( d );
      float*    buff_e      = ( float*    ) FLA_FLOAT_PTR( e );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      scomplex* buff_A      = ( scomplex* ) FLA_COMPLEX_PTR( A );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      float*    buff_work   = ( float*    ) FLA_FLOAT_PTR( work );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );

      F77_cstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vl, &vu,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d      = ( double*   ) FLA_DOUBLE_PTR( d );
      double*   buff_e      = ( double*   ) FLA_DOUBLE_PTR( e );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      dcomplex* buff_A      = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      double*   buff_work   = ( double*   ) FLA_DOUBLE_PTR( work );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );

      F77_zstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vl, &vu,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 

    }
  }

  FLA_Obj_free( &isuppz );
  FLA_Obj_free( &work );
  FLA_Obj_free( &iwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

