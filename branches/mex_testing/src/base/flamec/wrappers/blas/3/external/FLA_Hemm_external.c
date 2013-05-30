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

FLA_Error FLA_Hemm_external( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Datatype datatype;
  int          m_C, n_C;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          rs_C, cs_C;
  side_t       blis_side;
  uplo_t       blis_uplo; 

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Hemm_check( side, uplo, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );

  FLA_Param_map_flame_to_blis_side( side, &blis_side );
  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  

  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_C     = ( float * ) FLA_FLOAT_PTR( C );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    bli_ssymm( blis_side,
               blis_uplo,
               m_C,
               n_C,
               buff_alpha,
               buff_A, rs_A, cs_A,
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_C     = ( double * ) FLA_DOUBLE_PTR( C );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    bli_dsymm( blis_side,
               blis_uplo,
               m_C,
               n_C,
               buff_alpha,
               buff_A, rs_A, cs_A,
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B     = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_C     = ( scomplex * ) FLA_COMPLEX_PTR( C );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );

    bli_chemm( blis_side,
               blis_uplo,
               m_C,
               n_C,
               buff_alpha,
               buff_A, rs_A, cs_A,
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_C     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( C );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bli_zhemm( blis_side,
               blis_uplo,
               m_C,
               n_C,
               buff_alpha,
               buff_A, rs_A, cs_A,
               buff_B, rs_B, cs_B,
               buff_beta,
               buff_C, rs_C, cs_C );

    break;
  }

  }
 
  return FLA_SUCCESS;
}

