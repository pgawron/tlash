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


static unsigned int   flash_queue_stack           = 0;
static FLA_Bool       flash_queue_enabled         = TRUE;

void FLASH_Queue_begin( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_begin

----------------------------------------------------------------------------*/
{
   // Push onto the stack.
   flash_queue_stack++;

   return;
}


void FLASH_Queue_end( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_end

----------------------------------------------------------------------------*/
{
   // Pop off the stack.
   flash_queue_stack--;

   return;
}

FLA_Bool FLASH_Queue_get_enabled( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_get_enabled

----------------------------------------------------------------------------*/
{
   return FALSE;
}

FLA_Error FLASH_Queue_enable( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_enable

----------------------------------------------------------------------------*/
{

   // Raise an exception when SuperMatrix is not configured.
   FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED );
   return FLA_FAILURE;
}


FLA_Error FLASH_Queue_disable( void )
/*----------------------------------------------------------------------------

   FLASH_Queue_disable

----------------------------------------------------------------------------*/
{
   // Allow disabling enqueuing even when SuperMatrix is not configured.
   flash_queue_enabled = FALSE;
   return FLA_SUCCESS;
}

