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

extern fla_axpyt_t* fla_axpyt_cntl_blas;

FLA_Error FLA_Axpyt( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;

#ifdef FLA_ENABLE_BLAS1_FRONT_END_CNTL_TREES
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpyt_check( trans, alpha, A, B );

  // Invoke FLA_Axpyt_internal() with flat control tree that simply calls
  // external wrapper.
  r_val = FLA_Axpyt_internal( trans, alpha, A, B, fla_axpyt_cntl_blas );

#else
  r_val = FLA_Axpyt_external( trans, alpha, A, B );
#endif

  return r_val;
}

