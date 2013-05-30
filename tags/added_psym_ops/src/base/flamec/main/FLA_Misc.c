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



FLA_Error FLA_Obj_copy_view( FLA_Obj A, FLA_Obj* B )
{
  FLA_Obj  A_view;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_copy_view_check( A, B );

  // Set the m_inner and n_inner fields of a temporary copy of A.
  A_view         = A;
  A_view.m_inner = FLASH_Obj_scalar_length( A );
  A_view.n_inner = FLASH_Obj_scalar_width( A ); 

  // Copy the modified view into B. 
  *B = A_view; 

  return FLA_SUCCESS;
}



void FLA_Obj_extract_real_scalar( FLA_Obj alpha, double* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
    *alpha_value = ( double ) *FLA_FLOAT_PTR( alpha );
  else
    *alpha_value = *FLA_DOUBLE_PTR( alpha );
}



void FLA_Obj_extract_complex_scalar( FLA_Obj alpha, dcomplex* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_complex_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
  {
    scomplex temp = *FLA_COMPLEX_PTR( alpha );
    alpha_value->real = ( double ) temp.real;
    alpha_value->imag = ( double ) temp.imag;
  }
  else
    *alpha_value = *FLA_DOUBLE_COMPLEX_PTR( alpha );
}



void FLA_Obj_extract_real_part( FLA_Obj alpha, FLA_Obj beta )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_part_check( alpha, beta );

  if ( FLA_Obj_is_real( alpha ) )
    FLA_Copy( alpha, beta );
  else // if ( FLA_Obj_is_complex( alpha ) )
  {
    if      ( FLA_Obj_datatype( alpha ) == FLA_COMPLEX )
    {
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );
      float*    buff_beta  = FLA_FLOAT_PTR( beta );

      *buff_beta = buff_alpha->real;
    }
    else if ( FLA_Obj_datatype( alpha ) == FLA_DOUBLE_COMPLEX )
    {
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );
      double*   buff_beta  = FLA_DOUBLE_PTR( beta );

      *buff_beta = buff_alpha->real;
    }
  }
}



void FLA_Obj_extract_imag_part( FLA_Obj alpha, FLA_Obj beta )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_imag_part_check( alpha, beta );

  if ( FLA_Obj_is_real( alpha ) )
    FLA_Set( FLA_ZERO, beta );
  else // if ( FLA_Obj_is_complex( alpha ) )
  {
    if      ( FLA_Obj_datatype( alpha ) == FLA_COMPLEX )
    {
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );
      float*    buff_beta  = FLA_FLOAT_PTR( beta );

      *buff_beta = buff_alpha->imag;
    }
    else if ( FLA_Obj_datatype( alpha ) == FLA_DOUBLE_COMPLEX )
    {
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );
      double*   buff_beta  = FLA_DOUBLE_PTR( beta );

      *buff_beta = buff_alpha->imag;
    }
  }
}



void FLA_Obj_set_real_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_real_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
  }
}



void FLA_Obj_set_imag_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_imag_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
  }
}



FLA_Error FLA_Obj_fshow( FILE* file, char *s1, FLA_Obj A, char *format, char *s2 )
{
  FLA_Datatype datatype;
  dim_t        i, j, m, n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_fshow_check( file, s1, A, format, s2 );

  datatype = FLA_Obj_datatype( A );
  m        = FLA_Obj_length( A );
  n        = FLA_Obj_width( A );
  rs       = FLA_Obj_row_stride( A );
  cs       = FLA_Obj_col_stride( A );

  fprintf( file, "%s\n", s1 );

  switch ( datatype ){

  case FLA_CONSTANT:
  {
    int*      consti = FLA_INT_PTR( A );
    float*    consts = FLA_FLOAT_PTR( A );
    double*   constd = FLA_DOUBLE_PTR( A );
    scomplex* constc = FLA_COMPLEX_PTR( A );
    dcomplex* constz = FLA_DOUBLE_COMPLEX_PTR( A );

    fprintf( file, "int      = %d\n", *consti );
    fprintf( file, "float    = %e\n", *consts );
    fprintf( file, "double   = %e\n", *constd );
    fprintf( file, "scomplex = %e + %e\n", constc->real, constc->imag );
    fprintf( file, "dcomplex = %e + %e\n", constz->real, constc->imag );

    break;
  }

  case FLA_FLOAT:
  {
    float *buffer = ( float * ) FLA_FLOAT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double *buffer = ( double * ) FLA_DOUBLE_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buffer = ( scomplex * ) FLA_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buffer = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_INT:
  {
    int *buffer = ( int * ) FLA_INT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  }

  fprintf( file, "%s\n", s2 );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_show( char *s1, FLA_Obj A, char *format, char *s2 )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_show_check( s1, A, format, s2 );

  return FLA_Obj_fshow( stdout, s1, A, format, s2 );
}

