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

extern fla_eig_gest_t* flash_eig_gest_cntl;

FLA_Error FLASH_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj   Y;
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

  // The temporary matrix object Y must exist when execution occurs, NOT just
  // when enqueuing occurs. So if the SuperMatrix stack depth is positive, then
  // it means the user has declared a "parallel region" in his code, and thus
  // execution won't occur until sometime after FLASH_Eig_gest() returns, at
  // which time Y will have been deallocated. Thus, we disallow this scenario
  // for now, until we can think of a more general solution.
  if ( FLASH_Queue_stack_depth() > 0 )
  {
     FLA_Print_message( "FLASH_Eig_gest() MUST be invoked with standalone parallelism, and may not be called from within a user-level parallel region",
                        __FILE__, __LINE__ );
     FLA_Abort();
  }


  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, flash_eig_gest_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  FLASH_Obj_free( &Y );

  return r_val;
}

