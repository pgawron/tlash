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

FLA_Error FLA_Svd( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error r_val      = FLA_SUCCESS;
  dim_t     n_iter_max = 30;
  dim_t     k_accum    = 32;
  dim_t     b_alg      = 512;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Svd_check( jobu, jobv, A, s, U, V );

  if ( jobu != FLA_SVD_VECTORS_ALL || jobv != FLA_SVD_VECTORS_ALL )
    FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

  if ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) )
    r_val = FLA_Svd_uv_unb_var1( n_iter_max, A, s, U, V, k_accum, b_alg );
  else
    FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

  return r_val;
}
