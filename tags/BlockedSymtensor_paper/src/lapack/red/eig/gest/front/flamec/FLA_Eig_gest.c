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

extern fla_eig_gest_t* fla_eig_gest_ix_cntl;
extern fla_eig_gest_t* fla_eig_gest_nx_cntl;

FLA_Error FLA_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj   Y;
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );

  // Invoke FLA_Eig_gest_internal() with the appropriate control tree.
  if ( inv == FLA_INVERSE )
    r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, fla_eig_gest_ix_cntl );
  else
    r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, fla_eig_gest_nx_cntl );

  FLA_Obj_free( &Y );

  return r_val;
}

