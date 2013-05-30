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

FLA_Error FLA_Trsv_external( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x ) 
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_x;
  uplo_t       blis_uplo;
  trans_t      blis_trans;
  diag_t       blis_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );
  FLA_Param_map_flame_to_blis_diag( diag, &blis_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bli_strsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bli_dtrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bli_ctrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bli_ztrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  }

  return FLA_SUCCESS;
}

