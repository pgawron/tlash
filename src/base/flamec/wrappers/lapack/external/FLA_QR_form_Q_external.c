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

FLA_Error FLA_QR_form_Q_external( FLA_Obj A, FLA_Obj t )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          m_A, n_A, k_A;
  int          cs_A;
  int          lwork;
  FLA_Obj      work;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_QR_form_Q_check( A, t );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  k_A      = FLA_Obj_vector_dim( t );

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

    switch( datatype ){

    case FLA_FLOAT:
    {
      float*    buff_A    = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buff_t    = ( float    * ) FLA_FLOAT_PTR( t );
      float*    buff_work = ( float    * ) FLA_FLOAT_PTR( work );

      F77_sorgqr( &m_A,
                  &n_A,
                  &k_A,
                  buff_A, &cs_A,
                  buff_t,
                  buff_work, &lwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A    = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buff_t    = ( double   * ) FLA_DOUBLE_PTR( t );
      double*   buff_work = ( double   * ) FLA_DOUBLE_PTR( work );

      F77_dorgqr( &m_A,
                  &n_A,
                  &k_A,
                  buff_A, &cs_A,
                  buff_t,
                  buff_work, &lwork,
                  &info );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buff_t    = ( scomplex * ) FLA_COMPLEX_PTR( t );
      scomplex* buff_work = ( scomplex * ) FLA_COMPLEX_PTR( work );

      F77_cungqr( &m_A,
                  &n_A,
                  &k_A,
                  buff_A, &cs_A,
                  buff_t,
                  buff_work, &lwork,
                  &info );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
      dcomplex *buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work );

      F77_zungqr( &m_A,
                  &n_A,
                  &k_A,
                  buff_A, &cs_A,
                  buff_t,
                  buff_work, &lwork,
                  &info );

      break;
    }

    }

  }

  FLA_Obj_free( &work );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

