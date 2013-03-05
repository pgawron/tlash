/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2008, The University of Texas

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

#ifndef FLASH_QUEUE_MAIN_PROTOTYPES_H
#define FLASH_QUEUE_MAIN_PROTOTYPES_H


void           FLASH_Queue_begin( void );
void           FLASH_Queue_end( void );
unsigned int   FLASH_Queue_stack_depth( void );

FLA_Error      FLASH_Queue_enable( void );
FLA_Error      FLASH_Queue_disable( void );
FLA_Bool       FLASH_Queue_get_enabled( void );

void           FLASH_Queue_set_num_threads( unsigned int n_threads );
unsigned int   FLASH_Queue_get_num_threads( void );


#ifdef FLA_ENABLE_SUPERMATRIX


void           FLASH_Queue_init( void );
void           FLASH_Queue_finalize( void );

unsigned int   FLASH_Queue_get_num_tasks( void );

void           FLASH_Queue_set_verbose_output( FLASH_Verbose verbose );
FLASH_Verbose  FLASH_Queue_get_verbose_output( void );
void           FLASH_Queue_set_sorting( FLA_Bool sorting );
FLA_Bool       FLASH_Queue_get_sorting( void );
void           FLASH_Queue_set_caching( FLA_Bool caching );
FLA_Bool       FLASH_Queue_get_caching( void );
void           FLASH_Queue_set_work_stealing( FLA_Bool work_stealing );
FLA_Bool       FLASH_Queue_get_work_stealing( void );
void           FLASH_Queue_set_data_affinity( FLASH_Data_aff data_affinity );
FLASH_Data_aff FLASH_Queue_get_data_affinity( void );
double         FLASH_Queue_get_total_time( void );
double         FLASH_Queue_get_parallel_time( void );

void           FLASH_Queue_exec( void );


// --- helper functions -------------------------------------------------------

void           FLASH_Queue_set_parallel_time( double dtime );
void           FLASH_Queue_set_block_size( dim_t size );
dim_t          FLASH_Queue_get_block_size( void );
void           FLASH_Queue_set_cache_size( dim_t size );
dim_t          FLASH_Queue_get_cache_size( void );
void           FLASH_Queue_set_cache_line_size( dim_t size );
dim_t          FLASH_Queue_get_cache_line_size( void );
void           FLASH_Queue_set_cores_per_cache( int cores );
int            FLASH_Queue_get_cores_per_cache( void );
void           FLASH_Queue_set_cores_per_queue( int cores );
int            FLASH_Queue_get_cores_per_queue( void );
void           FLASH_Queue_reset( void );
FLASH_Task*    FLASH_Queue_get_head_task( void );
FLASH_Task*    FLASH_Queue_get_tail_task( void );
void           FLASH_Queue_push( void *func, void *cntl, char *name,
                                 FLA_Bool enabled_gpu,
                                 int n_int_args, int n_fla_args,
                                 int n_input_args, int n_output_args, ... );
void           FLASH_Queue_push_input( FLA_Obj obj, FLASH_Task* t );
void           FLASH_Queue_push_output( FLA_Obj obj, FLASH_Task* t );
FLASH_Task*    FLASH_Task_alloc( void *func, void *cntl, char *name,
                                 FLA_Bool enabled_gpu,
                                 int n_int_args, int n_fla_args,
                                 int n_input_args, int n_output_args );
void           FLASH_Task_free( FLASH_Task *t );
void           FLASH_Queue_exec_task( FLASH_Task *t );
void           FLASH_Queue_verbose_output( void );

void           FLASH_Queue_init_tasks( void *arg );
void           FLASH_Queue_wait_enqueue( FLASH_Task *t, void *arg );
FLASH_Task*    FLASH_Queue_wait_dequeue( int queue, int cache, void *arg );
FLASH_Task*    FLASH_Queue_wait_dequeue_block( int queue, int cache, void *arg );
void           FLASH_Queue_update_cache( FLASH_Task *t, void *arg );
void           FLASH_Queue_update_cache_block( FLA_Obj obj, int cache, FLA_Bool output, void *arg );
void           FLASH_Queue_prefetch( int cache, void *arg );
void           FLASH_Queue_prefetch_block( FLA_Obj obj );
FLASH_Task*    FLASH_Queue_work_stealing( int queue, void *arg );
#ifdef FLA_ENABLE_GPU
void           FLASH_Queue_create_gpu( int thread, void *arg );
void           FLASH_Queue_destroy_gpu( int thread, void *arg );
FLA_Bool       FLASH_Queue_exec_gpu( FLASH_Task *t, void *arg );
FLA_Bool       FLASH_Queue_check_gpu( FLASH_Task *t, void *arg );
FLA_Bool       FLASH_Queue_check_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_update_gpu( FLASH_Task *t, void **input_arg, void **output_arg, void *arg );
void           FLASH_Queue_update_block_gpu( FLA_Obj obj, void **buffer_gpu, int thread, void *arg );
void           FLASH_Queue_mark_gpu( FLASH_Task *t, void *arg );
void           FLASH_Queue_invalidate_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_flush_block_gpu( FLA_Obj obj, int thread, void *arg );
void           FLASH_Queue_flush_gpu( int thread, void *arg );
#endif
void           FLASH_Queue_exec_parallel( void *arg );
void*          FLASH_Queue_exec_parallel_function( void *arg );
FLASH_Task*    FLASH_Task_update_dependencies( FLASH_Task *t, void *arg );
FLASH_Task*    FLASH_Task_update_binding( FLASH_Task *t, FLASH_Task *r, void *arg );
void           FLASH_Task_free_parallel( FLASH_Task *t, void *arg );

void           FLASH_Queue_exec_simulation( void *arg );


#endif // FLA_ENABLE_SUPERMATRIX


#endif // FLASH_QUEUE_MAIN_PROTOTYPES_H
