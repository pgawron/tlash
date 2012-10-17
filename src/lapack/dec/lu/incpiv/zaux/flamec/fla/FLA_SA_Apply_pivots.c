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

FLA_Error FLA_SA_Apply_pivots( FLA_Obj C, FLA_Obj E, FLA_Obj p )
{
  FLA_Datatype datatype;
  int          m_C, n_C, cs_C;
  int                    cs_E;
  // int                    rs_C;
  // int                    rs_E;
  int          m_p;
  int          i;
  int*         buff_p;

  if ( FLA_Obj_has_zero_dim( C ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( C );

  m_C    = FLA_Obj_length( C );
  n_C    = FLA_Obj_width( C );
  cs_C   = FLA_Obj_col_stride( C );
  // rs_C   = FLA_Obj_row_stride( C );

  cs_E   = FLA_Obj_col_stride( E );
  // rs_E   = FLA_Obj_row_stride( E );

  m_p    = FLA_Obj_length( p );
  
  buff_p = ( int * ) FLA_INT_PTR( p );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_C = ( float * ) FLA_FLOAT_PTR( C );
    float* buff_E = ( float * ) FLA_FLOAT_PTR( E );

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        bli_sswap( n_C, 
                   buff_C + 0*cs_C + i,                         cs_C, 
                   buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_C = ( double * ) FLA_DOUBLE_PTR( C );
    double* buff_E = ( double * ) FLA_DOUBLE_PTR( E );

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        bli_dswap( n_C, 
                   buff_C + 0*cs_C + i,                         cs_C, 
                   buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex* buff_C = ( scomplex * ) FLA_COMPLEX_PTR( C );
    scomplex* buff_E = ( scomplex * ) FLA_COMPLEX_PTR( E );

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        bli_cswap( n_C, 
                   buff_C + 0*cs_C + i,                         cs_C, 
                   buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_C = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( C );
    dcomplex* buff_E = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( E );

    for ( i = 0; i < m_p; ++i )
    {
      if ( buff_p[ i ] != 0 ) 
        bli_zswap( n_C, 
                   buff_C + 0*cs_C + i,                         cs_C, 
                   buff_E + 0*cs_E + buff_p[ i ] - ( m_C - i ), cs_E );
    }
    break;
  }

  }

  return FLA_SUCCESS;
}
