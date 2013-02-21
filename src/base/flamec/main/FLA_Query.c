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

FLA_Datatype FLA_Obj_datatype( FLA_Obj obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_datatype_check( obj );

  return obj.base->datatype;
}



FLA_Datatype FLA_Obj_datatype_proj_to_real( FLA_Obj A )
{
	FLA_Datatype datatype;

	if ( FLA_Obj_is_single_precision( A ) )
		datatype = FLA_FLOAT;
	else
		datatype = FLA_DOUBLE;

	return datatype;
}



FLA_Datatype FLA_Obj_datatype_proj_to_complex( FLA_Obj A )
{
	FLA_Datatype datatype;

	if ( FLA_Obj_is_single_precision( A ) )
		datatype = FLA_COMPLEX;
	else
		datatype = FLA_DOUBLE_COMPLEX;

	return datatype;
}



FLA_Elemtype FLA_Obj_elemtype( FLA_Obj obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_elemtype_check( obj );

  return obj.base->elemtype;
}



dim_t FLA_Obj_datatype_size( FLA_Datatype datatype )
{
  dim_t datatype_size = 0;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_datatype_size_check( datatype );

  switch( datatype )
  {
    case FLA_INT: 
      datatype_size = sizeof( int );
      break;
    case FLA_FLOAT: 
      datatype_size = sizeof( float );
      break;
    case FLA_DOUBLE: 
      datatype_size = sizeof( double );
      break;
    case FLA_COMPLEX: 
      datatype_size = sizeof( scomplex );
      break;
    case FLA_DOUBLE_COMPLEX: 
      datatype_size = sizeof( dcomplex );
      break;
    case FLA_CONSTANT: 
      datatype_size = FLA_CONSTANT_SIZE;
      break;
  }

  return datatype_size;
}



dim_t FLA_Obj_elem_size( FLA_Obj obj )
{
  dim_t elem_size = 0;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_elem_size_check( obj );

  if ( FLA_Obj_elemtype( obj ) == FLA_TENSOR )
  {
    elem_size = sizeof( FLA_Obj );
  }
  else if ( FLA_Obj_elemtype( obj ) == FLA_MATRIX )
  {
    elem_size = sizeof( FLA_Obj );
  }
  else // if ( FLA_Obj_elemtype( obj ) == FLA_SCALAR )
  {
    elem_size = FLA_Obj_datatype_size( FLA_Obj_datatype( obj ) );
  }

  return elem_size;
}

dim_t FLA_Obj_length( FLA_Obj obj )
{
  return obj.m;
}

dim_t FLA_Obj_width( FLA_Obj obj )
{
  return obj.n;
}

FLA_Uplo FLA_Obj_structure( FLA_Obj obj )
{
  return obj.base->uplo;
}

dim_t FLA_Obj_vector_dim( FLA_Obj obj )
{
  return ( obj.m == 1 ? obj.n
                      : obj.m );
}



dim_t FLA_Obj_vector_inc( FLA_Obj obj )
{
  return ( obj.m == 1 ? (obj.base)->cs
                      : (obj.base)->rs );
}



dim_t FLA_Obj_min_dim( FLA_Obj obj )
{
  return min( obj.m, obj.n );
}



dim_t FLA_Obj_max_dim( FLA_Obj obj )
{
  return max( obj.m, obj.n );
}



dim_t FLA_Obj_row_stride( FLA_Obj obj )
{
  return (obj.base)->rs;
}

dim_t FLA_Obj_col_stride( FLA_Obj obj )
{
  return (obj.base)->cs;
}

dim_t FLA_Obj_row_offset( FLA_Obj obj )
{
	return obj.offm;
}


dim_t FLA_Obj_col_offset( FLA_Obj obj )
{
	return obj.offn;
}

dim_t FLA_Obj_base_length( FLA_Obj obj )
{
	return (obj.base)->m;
}


dim_t FLA_Obj_base_width( FLA_Obj obj )
{
	return (obj.base)->n;
}


dim_t FLA_Obj_num_elem_alloc( FLA_Obj obj )
{
  return (obj.base)->n_elem_alloc;
}


void* FLA_Obj_base_buffer( FLA_Obj obj )
{
  return (obj.base)->buffer;
}

void* FLA_Obj_buffer_at_view( FLA_Obj obj )
{
  char*  buffer;
  size_t elem_size, offm, offn, rs, cs;
  size_t byte_offset;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_buffer_at_view_check( obj );

  elem_size   = ( size_t ) FLA_Obj_elem_size( obj );
  rs          = ( size_t ) FLA_Obj_row_stride( obj );
  cs          = ( size_t ) FLA_Obj_col_stride( obj );
  offm        = ( size_t ) obj.offm;
  offn        = ( size_t ) obj.offn;

  byte_offset = elem_size * ( offn * cs + offm * rs );

  buffer      = ( char * ) (obj.base)->buffer;

  return ( void* ) ( buffer + byte_offset );
}


FLA_Bool FLA_Obj_buffer_is_null( FLA_Obj obj )
{
  FLA_Bool r_val;

  if      ( obj.base == NULL )
    r_val = TRUE;
  else if ( (obj.base)->buffer == NULL )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_int( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_INT )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_floating_point( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX ||
       datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_constant( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_CONSTANT )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_real( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_CONSTANT || datatype == FLA_FLOAT || datatype == FLA_DOUBLE )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_complex( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_CONSTANT || datatype == FLA_COMPLEX || datatype == FLA_DOUBLE_COMPLEX )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_single_precision( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_CONSTANT || datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_double_precision( FLA_Obj A )
{
  FLA_Datatype datatype;
  FLA_Bool     r_val;

  datatype = FLA_Obj_datatype( A );

  if ( datatype == FLA_CONSTANT || datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
    r_val = TRUE;
  else
    r_val = FALSE;

  return r_val;
}



FLA_Bool FLA_Obj_is_scalar( FLA_Obj A )
{
  FLA_Bool r_val = FALSE;

  if ( FLA_Obj_length( A ) == 1 &&
       FLA_Obj_width( A ) == 1 )
    r_val = TRUE;

  return r_val;
}

FLA_Bool FLA_Obj_is_vector( FLA_Obj A )
{
  FLA_Bool r_val = FALSE;

  if ( FLA_Obj_length( A ) == 1 || FLA_Obj_width( A ) == 1 )
    r_val = TRUE;

  return r_val;
}



FLA_Bool FLA_Obj_has_zero_dim( FLA_Obj A )
{
  FLA_Bool r_val = FALSE;

  if ( FLA_Obj_length( A ) == 0 || FLA_Obj_width( A ) == 0 )
    r_val = TRUE;

  return r_val;
}



FLA_Bool FLA_Obj_is_col_major( FLA_Obj A )
{
  FLA_Bool r_val = FALSE;

  // A row stride of 1 indicates column-major storage.
  if ( FLA_Obj_row_stride( A ) == 1 )
    r_val = TRUE;

  return r_val;
}



FLA_Bool FLA_Obj_is_row_major( FLA_Obj A )
{
  FLA_Bool r_val = FALSE;

  // A column stride of 1 indicates row-major storage.
  if ( FLA_Obj_col_stride( A ) == 1 )
    r_val = TRUE;

  return r_val;
}



FLA_Bool FLA_Obj_is_conformal_to( FLA_Trans trans, FLA_Obj A, FLA_Obj B )
{
  FLA_Bool r_val = TRUE;

  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    if ( FLA_Obj_length( A ) != FLA_Obj_length( B ) ||
         FLA_Obj_width( A )  != FLA_Obj_width( B ) )
      r_val = FALSE;
  }
  else
  {
    if ( FLA_Obj_width( A )  != FLA_Obj_length( B ) ||
         FLA_Obj_length( A ) != FLA_Obj_width( B ) )
      r_val = FALSE;
  }

  return r_val;
}



FLA_Bool FLA_Obj_is( FLA_Obj A, FLA_Obj B )
{
  FLA_Bool r_val = FALSE;

  if ( A.base == B.base )
    r_val = TRUE;

  return r_val;
}



FLA_Bool FLA_Obj_equals( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  FLA_Datatype datatype;
  dim_t        m, n;
  dim_t        rs_A, cs_A;
  dim_t        rs_B, cs_B;
  dim_t        i, j;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_equals_check( A, B );

  m      = FLA_Obj_length( A );
  n      = FLA_Obj_width( A );
  rs_A   = FLA_Obj_row_stride( A );
  cs_A   = FLA_Obj_col_stride( A );
  rs_B   = FLA_Obj_row_stride( B );
  cs_B   = FLA_Obj_col_stride( B );

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  // If A is a non-FLA_CONSTANT object, then we should proceed based on the
  // value of datatype_A. In such a situation, either datatype_B is an exact
  // match and we're fine, or datatype_B is FLA_CONSTANT, in which case we're
  // also covered since FLA_CONSTANT encompassas all numerical types.
  // If A is an FLA_CONSTANT object, then we should proceed based on the value
  // of datatype_B. In this case, datatype_B is either a non-FLA_CONSTANT type,
  // which mirrors the second sub-case above, or datatype_B is FLA_CONSTANT,
  // in which case both types are FLA_CONSTANT and therefore we have to handle
  // that case. Only if both are FLA_CONSTANTs does the FLA_CONSTANT case
  // statement below execute.
  if ( datatype_A != FLA_CONSTANT )
    datatype = datatype_A;
  else
    datatype = datatype_B;

  switch ( datatype )
  {
    case FLA_CONSTANT:
    {
      // We require ALL floating-point fields to be the same.
      float*    buffs_A = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buffs_B = ( float    * ) FLA_FLOAT_PTR( B );
      double*   buffd_A = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buffd_B = ( double   * ) FLA_DOUBLE_PTR( B );
      scomplex* buffc_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buffc_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      dcomplex* buffz_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buffz_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      if ( *buffs_A != *buffs_B ||
           *buffd_A != *buffd_B ||
           buffc_A->real != buffc_B->real ||
           buffc_A->imag != buffc_B->imag ||
           buffz_A->real != buffz_B->real ||
           buffz_A->imag != buffz_B->imag )
      {
        return FALSE;
      }

      break;
    }

    case FLA_INT:
    {
      int *buff_A = ( int * ) FLA_INT_PTR( A );
      int *buff_B = ( int * ) FLA_INT_PTR( B );

      for ( j = 0; j < n; j++ )
        for ( i = 0; i < m; i++ )
          if ( buff_A[ j * cs_A + i * rs_A ] != 
               buff_B[ j * cs_B + i * rs_B ] )
          {
            return FALSE;
          }

      break;
    }

    case FLA_FLOAT:
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

      for ( j = 0; j < n; j++ )
        for ( i = 0; i < m; i++ )
          if ( buff_A[ j * cs_A + i * rs_A ] != 
               buff_B[ j * cs_B + i * rs_B ] )
          {
            return FALSE;
          }

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

      for ( j = 0; j < n; j++ )
        for ( i = 0; i < m; i++ )
          if ( buff_A[ j * cs_A + i * rs_A ] != 
               buff_B[ j * cs_B + i * rs_B ] )
          {
            return FALSE;
          }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n; j++ )
        for ( i = 0; i < m; i++ )
          if ( buff_A[ j * cs_A + i * rs_A ].real != buff_B[ j * cs_B + i * rs_B ].real ||
               buff_A[ j * cs_A + i * rs_A ].imag != buff_B[ j * cs_B + i * rs_B ].imag )
          {
            return FALSE;
          }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n; j++ )
        for ( i = 0; i < m; i++ )
          if ( buff_A[ j * cs_A + i * rs_A ].real != buff_B[ j * cs_B + i * rs_B ].real ||
               buff_A[ j * cs_A + i * rs_A ].imag != buff_B[ j * cs_B + i * rs_B ].imag )
          {
            return FALSE;
          }

      break;
    }
  }

  return TRUE;
}



FLA_Bool FLA_Obj_gt( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_gt_check( A, B );

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  if ( datatype_A != FLA_CONSTANT ) datatype = datatype_A;
  else                              datatype = datatype_B;

  switch ( datatype )
  {
    case FLA_CONSTANT:
    {
      // We require ALL floating-point fields to be the same.
      float*    buffs_A = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buffs_B = ( float    * ) FLA_FLOAT_PTR( B );
      double*   buffd_A = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buffd_B = ( double   * ) FLA_DOUBLE_PTR( B );
      scomplex* buffc_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buffc_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      dcomplex* buffz_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buffz_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      if ( !( *buffs_A      > *buffs_B      &&
              *buffd_A      > *buffd_B      &&
              buffc_A->real > buffc_B->real &&
              buffc_A->imag > buffc_B->imag &&
              buffz_A->real > buffz_B->real &&
              buffz_A->imag > buffz_B->imag ) )
      {
        return FALSE;
      }

      break;
    }

    case FLA_INT:
    {
      int *buff_A = ( int * ) FLA_INT_PTR( A );
      int *buff_B = ( int * ) FLA_INT_PTR( B );

      if ( !( *buff_A > *buff_B ) ) return FALSE;

      break;
    }

    case FLA_FLOAT:
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

      if ( !( *buff_A > *buff_B ) ) return FALSE;

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

      if ( !( *buff_A > *buff_B ) ) return FALSE;

      break;
    }

  }

  return TRUE;
}


FLA_Bool FLA_Obj_ge( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_ge_check( A, B );

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  if ( datatype_A != FLA_CONSTANT ) datatype = datatype_A;
  else                              datatype = datatype_B;

  switch ( datatype )
  {
    case FLA_CONSTANT:
    {
      // We require ALL floating-point fields to be the same.
      float*    buffs_A = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buffs_B = ( float    * ) FLA_FLOAT_PTR( B );
      double*   buffd_A = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buffd_B = ( double   * ) FLA_DOUBLE_PTR( B );
      scomplex* buffc_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buffc_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      dcomplex* buffz_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buffz_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      if ( !( *buffs_A      >= *buffs_B      &&
              *buffd_A      >= *buffd_B      &&
              buffc_A->real >= buffc_B->real &&
              buffc_A->imag >= buffc_B->imag &&
              buffz_A->real >= buffz_B->real &&
              buffz_A->imag >= buffz_B->imag ) )
      {
        return FALSE;
      }

      break;
    }

    case FLA_INT:
    {
      int *buff_A = ( int * ) FLA_INT_PTR( A );
      int *buff_B = ( int * ) FLA_INT_PTR( B );

      if ( !( *buff_A >= *buff_B ) ) return FALSE;

      break;
    }

    case FLA_FLOAT:
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

      if ( !( *buff_A >= *buff_B ) ) return FALSE;

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

      if ( !( *buff_A >= *buff_B ) ) return FALSE;

      break;
    }

  }

  return TRUE;
}

FLA_Bool FLA_Obj_lt( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_lt_check( A, B );

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  if ( datatype_A != FLA_CONSTANT ) datatype = datatype_A;
  else                              datatype = datatype_B;

  switch ( datatype )
  {
    case FLA_CONSTANT:
    {
      // We require ALL floating-point fields to be the same.
      float*    buffs_A = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buffs_B = ( float    * ) FLA_FLOAT_PTR( B );
      double*   buffd_A = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buffd_B = ( double   * ) FLA_DOUBLE_PTR( B );
      scomplex* buffc_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buffc_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      dcomplex* buffz_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buffz_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      if ( !( *buffs_A      < *buffs_B      &&
              *buffd_A      < *buffd_B      &&
              buffc_A->real < buffc_B->real &&
              buffc_A->imag < buffc_B->imag &&
              buffz_A->real < buffz_B->real &&
              buffz_A->imag < buffz_B->imag ) )
      {
        return FALSE;
      }

      break;
    }

    case FLA_INT:
    {
      int *buff_A = ( int * ) FLA_INT_PTR( A );
      int *buff_B = ( int * ) FLA_INT_PTR( B );

      if ( !( *buff_A < *buff_B ) ) return FALSE;

      break;
    }

    case FLA_FLOAT:
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

      if ( !( *buff_A < *buff_B ) ) return FALSE;

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

      if ( !( *buff_A < *buff_B ) ) return FALSE;

      break;
    }

  }

  return TRUE;
}

FLA_Bool FLA_Obj_le( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  FLA_Datatype datatype;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_le_check( A, B );

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  if ( datatype_A != FLA_CONSTANT ) datatype = datatype_A;
  else                              datatype = datatype_B;

  switch ( datatype )
  {
    case FLA_CONSTANT:
    {
      // We require ALL floating-point fields to be the same.
      float*    buffs_A = ( float    * ) FLA_FLOAT_PTR( A );
      float*    buffs_B = ( float    * ) FLA_FLOAT_PTR( B );
      double*   buffd_A = ( double   * ) FLA_DOUBLE_PTR( A );
      double*   buffd_B = ( double   * ) FLA_DOUBLE_PTR( B );
      scomplex* buffc_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex* buffc_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      dcomplex* buffz_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buffz_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

      if ( !( *buffs_A      <= *buffs_B      &&
              *buffd_A      <= *buffd_B      &&
              buffc_A->real <= buffc_B->real &&
              buffc_A->imag <= buffc_B->imag &&
              buffz_A->real <= buffz_B->real &&
              buffz_A->imag <= buffz_B->imag ) )
      {
        return FALSE;
      }

      break;
    }

    case FLA_INT:
    {
      int *buff_A = ( int * ) FLA_INT_PTR( A );
      int *buff_B = ( int * ) FLA_INT_PTR( B );

      if ( !( *buff_A <= *buff_B ) ) return FALSE;

      break;
    }

    case FLA_FLOAT:
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );

      if ( !( *buff_A <= *buff_B ) ) return FALSE;

      break;
    }

    case FLA_DOUBLE:
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );

      if ( !( *buff_A <= *buff_B ) ) return FALSE;

      break;
    }

  }

  return TRUE;
}



void* FLA_Submatrix_at( FLA_Datatype datatype, void* buffer, dim_t i, dim_t j, dim_t rs, dim_t cs )
{
  void* r_val = buffer;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Submatrix_at_check( datatype, buffer, i, j, rs, cs );

  switch( datatype )
  {
    case FLA_INT:
      r_val = ( void* ) ( (      ( int* ) buffer ) + i * rs + j * cs );
      break;

    case FLA_FLOAT:
      r_val = ( void* ) ( (    ( float* ) buffer ) + i * rs + j * cs );
      break;

    case FLA_DOUBLE:
      r_val = ( void* ) ( (   ( double* ) buffer ) + i * rs + j * cs );
      break;

    case FLA_COMPLEX:
      r_val = ( void* ) ( ( ( scomplex* ) buffer ) + i * rs + j * cs );
      break;

    case FLA_DOUBLE_COMPLEX:
      r_val = ( void* ) ( ( ( dcomplex* ) buffer ) + i * rs + j * cs );
      break;
  }

  return r_val;
}

