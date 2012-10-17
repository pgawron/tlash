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

FLA_Error FLA_Eig_gest_blk_external( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Error    r_val = FLA_SUCCESS;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  int          itype;
  int          info;
  FLA_Datatype datatype;
  int          m_A, cs_A;
  int          cs_B;
  char         blas_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

//  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  if ( inv == FLA_INVERSE )
    itype  = 1;
  else
    itype  = 2;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

    F77_ssygst( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

    F77_dsygst( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );

    F77_chegst( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    F77_zhegst( &itype,
                &blas_uplo,
                &m_A,
                buff_A, &cs_A,
                buff_B, &cs_B,
                &info );

    break;
  } 

  }

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return r_val;
}

FLA_Error FLA_Eig_gest_il_blk_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_iu_blk_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_nl_blk_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
}

FLA_Error FLA_Eig_gest_nu_blk_ext( FLA_Obj A, FLA_Obj B )
{
  return FLA_Eig_gest_blk_external( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
}

