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

FLA_Error FLA_Axpys_external( FLA_Obj alpha0, FLA_Obj alpha1, FLA_Obj A, FLA_Obj beta, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  trans_t      blis_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Axpys_check( alpha0, alpha1, A, beta, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_blis_trans( FLA_NO_TRANSPOSE, &blis_trans );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_alpha0 = ( float * ) FLA_FLOAT_PTR( alpha0 );
    float *buff_alpha1 = ( float * ) FLA_FLOAT_PTR( alpha1 );
    float *buff_beta   = ( float * ) FLA_FLOAT_PTR( beta );
    float *buff_A      = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B      = ( float * ) FLA_FLOAT_PTR( B );

    bli_saxpysmt( blis_trans,
                  m_B,
                  n_B,
                  buff_alpha0,
                  buff_alpha1,
                  buff_A, rs_A, cs_A,
                  buff_beta,
                  buff_B, rs_B, cs_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_alpha0 = ( double * ) FLA_DOUBLE_PTR( alpha0 );
    double *buff_alpha1 = ( double * ) FLA_DOUBLE_PTR( alpha1 );
    double *buff_beta   = ( double * ) FLA_DOUBLE_PTR( beta );
    double *buff_A      = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B      = ( double * ) FLA_DOUBLE_PTR( B );

    bli_daxpysmt( blis_trans,
                  m_B,
                  n_B,
                  buff_alpha0,
                  buff_alpha1,
                  buff_A, rs_A, cs_A,
                  buff_beta,
                  buff_B, rs_B, cs_B );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_alpha0 = ( scomplex * ) FLA_COMPLEX_PTR( alpha0 );
    scomplex *buff_alpha1 = ( scomplex * ) FLA_COMPLEX_PTR( alpha1 );
    scomplex *buff_beta   = ( scomplex * ) FLA_COMPLEX_PTR( beta );
    scomplex *buff_A      = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B      = ( scomplex * ) FLA_COMPLEX_PTR( B );

    bli_caxpysmt( blis_trans,
                  m_B,
                  n_B,
                  buff_alpha0,
                  buff_alpha1,
                  buff_A, rs_A, cs_A,
                  buff_beta,
                  buff_B, rs_B, cs_B );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_alpha0 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha0 );
    dcomplex *buff_alpha1 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha1 );
    dcomplex *buff_beta   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );
    dcomplex *buff_A      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    bli_zaxpysmt( blis_trans,
                  m_B,
                  n_B,
                  buff_alpha0,
                  buff_alpha1,
                  buff_A, rs_A, cs_A,
                  buff_beta,
                  buff_B, rs_B, cs_B );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

