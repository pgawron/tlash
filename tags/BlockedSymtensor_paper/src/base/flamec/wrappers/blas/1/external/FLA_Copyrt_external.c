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

FLA_Error FLA_Copyrt_external( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype dt_A;
  FLA_Datatype dt_B;
  int          m_B, n_B;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  uplo_t       blis_uplo;
  trans_t      blis_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Copyrt_check( uplo, trans, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  dt_A     = FLA_Obj_datatype( A );
  dt_B     = FLA_Obj_datatype( B );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );

  // If A is of type FLA_CONSTANT, then we have to proceed based on the
  // datatype of B.
  if      ( dt_A == FLA_CONSTANT )
  {
    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bli_scopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bli_dcopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bli_ccopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bli_zcopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
  }
/*
  else if ( dt_A == FLA_INT )
  {
    int*      buff_A = ( int * ) FLA_INT_PTR( A );
    int*      buff_B = ( int * ) FLA_INT_PTR( B );

    bli_icopymrt( blis_uplo,
                  blis_trans,
                  m_B,
                  n_B,
                  buff_A, rs_A, cs_A,
                  buff_B, rs_B, cs_B );
  }
*/
  else if ( dt_A == FLA_FLOAT )
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bli_scopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bli_sdcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bli_sccopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bli_szcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_DOUBLE )
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bli_dscopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bli_dcopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bli_dccopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bli_dzcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_COMPLEX )
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bli_cscopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bli_cdcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bli_ccopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bli_czcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_DOUBLE_COMPLEX )
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bli_zscopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bli_zdcopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bli_zccopymrt( blis_uplo,
                     blis_trans,
                     m_B,
                     n_B,
                     buff_A, rs_A, cs_A,
                     buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bli_zcopymrt( blis_uplo,
                    blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
  }
  
  return FLA_SUCCESS;
}

