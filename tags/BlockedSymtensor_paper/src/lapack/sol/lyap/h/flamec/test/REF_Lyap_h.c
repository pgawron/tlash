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

FLA_Error REF_Lyap_h( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  FLA_Trans transA;

  if ( FLA_Obj_is_complex( A ) ) transA = FLA_CONJ_TRANSPOSE;
  else                           transA = FLA_TRANSPOSE;

  if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
    FLA_Scal( FLA_MINUS_ONE, C );

  FLA_Sylv_blk_external( transA, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );
  //FLA_Sylv( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, A, C, scale );

  return FLA_SUCCESS;
}
