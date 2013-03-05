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

FLA_Error FLA_LU_piv_unb_external( FLA_Obj A, FLA_Obj p )
{
  FLA_Error    r_val = FLA_SUCCESS;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  int          info;
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    int   *buff_p = ( int   * ) FLA_INT_PTR( p );

    F77_sgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    int    *buff_p = ( int    * ) FLA_INT_PTR( p );

    F77_dgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    int      *buff_p = ( int      * ) FLA_INT_PTR( p );

    F77_cgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    int      *buff_p = ( int      * ) FLA_INT_PTR( p );

    F77_zgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  }

  FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );

  // Convert to zero-based indexing, if an index was reported.
  if ( info > 0 ) r_val = info - 1;
  else            r_val = FLA_SUCCESS;

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return r_val;
}

FLA_Error FLA_LU_piv_unb_ext( FLA_Obj A, FLA_Obj p )
{
  return FLA_LU_piv_unb_external( A, p );
}

