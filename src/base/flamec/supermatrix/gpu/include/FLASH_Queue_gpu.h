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

#ifndef FLASH_QUEUE_GPU_H
#define FLASH_QUEUE_GPU_H

#ifdef FLA_ENABLE_GPU


void           FLASH_Queue_init_gpu( void );
void           FLASH_Queue_finalize_gpu( void );

FLA_Error      FLASH_Queue_enable_gpu( void );
FLA_Error      FLASH_Queue_disable_gpu( void );
FLA_Bool       FLASH_Queue_get_enabled_gpu( void );


// --- helper functions -------------------------------------------------------

void           FLASH_Queue_set_gpu_num_blocks( dim_t n_blocks );
dim_t          FLASH_Queue_get_gpu_num_blocks( void );

FLA_Error      FLASH_Queue_bind_gpu( int thread );
FLA_Error      FLASH_Queue_alloc_gpu( dim_t size, FLA_Datatype datatype, void** buffer_gpu );
FLA_Error      FLASH_Queue_free_gpu( void* buffer_gpu );
FLA_Error      FLASH_Queue_write_gpu( FLA_Obj obj, void* buffer_gpu );
FLA_Error      FLASH_Queue_read_gpu( FLA_Obj obj, void* buffer_gpu );

void           FLASH_Queue_exec_task_gpu( FLASH_Task* t, void** input_arg, void** output_arg );


#endif

#endif // FLASH_QUEUE_GPU_H
