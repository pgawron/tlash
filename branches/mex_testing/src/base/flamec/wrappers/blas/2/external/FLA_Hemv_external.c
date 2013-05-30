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

FLA_Error FLA_Hemv_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_x;
  int          inc_y;
  uplo_t       blis_uplo;
  conj_t       blis_conj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Hemv_check( uplo, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );
  inc_y    = FLA_Obj_vector_inc( y );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_conj( FLA_NO_CONJUGATE, &blis_conj );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x     = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_y     = ( float * ) FLA_FLOAT_PTR( y );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    bli_ssymv( blis_uplo,
               m_A,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,
               buff_y, inc_y );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x     = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_y     = ( double * ) FLA_DOUBLE_PTR( y );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    bli_dsymv( blis_uplo,
               m_A,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,
               buff_y, inc_y );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x     = ( scomplex * ) FLA_COMPLEX_PTR( x );
    scomplex *buff_y     = ( scomplex * ) FLA_COMPLEX_PTR( y );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );

    bli_chemv( blis_uplo,
               blis_conj,
               m_A,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,
               buff_y, inc_y );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    dcomplex *buff_y     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( y );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bli_zhemv( blis_uplo,
               blis_conj,
               m_A,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,
               buff_y, inc_y );

    break;
  }

  }

  return FLA_SUCCESS;
}

