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

FLA_Error FLA_Shift_pivots_to( FLA_Pivot_type ptype, FLA_Obj p )
{
  int  m_p, n_p;
  int* buff_p;
  int  i;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Shift_pivots_to_check( ptype, p );

  m_p    = FLA_Obj_length( p );
  n_p    = FLA_Obj_width( p );
  buff_p = FLA_INT_PTR( p );

  if ( m_p < 1 || n_p < 1 ) return FLA_SUCCESS;

  if ( ptype == FLA_LAPACK_PIVOTS )
  {
    // Shift FLAME pivots to LAPACK pivots.
    for ( i = 0; i < m_p; i++ )
      buff_p[ i ] += i + 1;
  }
  else
  {
    // Otherwise, shift LAPACK pivots back to FLAME.
    for ( i = 0; i < m_p; i++ )
      buff_p[ i ] -= i + 1;
  }

  return FLA_SUCCESS;
}

