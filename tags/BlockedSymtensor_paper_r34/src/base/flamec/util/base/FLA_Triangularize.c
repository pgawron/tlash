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

FLA_Error FLA_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  uplo_t       blis_uplo;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Triangularize_check( uplo, diag, A );

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  // We have to toggle the uplo parameter because we will use it to specify
  // which triangle to zero out.
  if ( uplo == FLA_LOWER_TRIANGULAR ) uplo = FLA_UPPER_TRIANGULAR;
  else                                uplo = FLA_LOWER_TRIANGULAR;

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_0 = ( float * ) FLA_FLOAT_PTR( FLA_ZERO );
    float *buff_1 = ( float * ) FLA_FLOAT_PTR( FLA_ONE );

    bli_ssetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bli_ssetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bli_ssetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_0 = ( double * ) FLA_DOUBLE_PTR( FLA_ZERO );
    double *buff_1 = ( double * ) FLA_DOUBLE_PTR( FLA_ONE );

    bli_dsetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bli_dsetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bli_dsetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_0 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_ZERO );
    scomplex *buff_1 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_ONE );

    bli_csetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bli_csetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bli_csetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_0 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
    dcomplex *buff_1 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );

    bli_zsetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bli_zsetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bli_zsetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  }

  return FLA_SUCCESS;
}

