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


typedef struct FLA_Queue_s      FLA_Queue;
typedef struct FLA_Task_s       FLA_Task;


struct FLA_Queue_s
{
	// Number of tasks currently in queue.
	int            n_tasks;

	// Pointers to head (front) and tail (back) of queue.
	FLA_Task*      head;
	FLA_Task*      tail;
	
	// An array that provides random-access to tasks.
	FLA_Task**     task_array;
};

struct FLA_Task_s
{
	// FLOP cost of task
	double         cost;
	
	// Function pointer
	void*          func;
	
	// Function integer arguments
	int            n_int_args;
	int*           int_arg;
	
	// Function FLA_Obj arguments
	int            n_fla_args;
	FLA_Obj*       fla_arg;
	
	// Control tree pointer argument
	void*          cntl_arg;
	
	// Support for linked list of FLA_Tasks
	FLA_Task*      next_task;
};


// -----------------------------------------------------------------------------


void      FLA_Queue_init( void );
void      FLA_Queue_exec( void );
void      FLA_Queue_finalize( void );

void      FLA_Queue_push( void* func, double cost, int n_int_params, int n_fla_params, ... );

void      FLA_Queue_set_num_threads( int n_threads );
int       FLA_Queue_get_num_threads( void );

void      FLA_queue_exec_sync( void );
void      FLA_queue_exec_task( FLA_Task* t );
void      FLA_queue_push_unsorted( FLA_Task* t );
void      FLA_queue_insert_sorted( FLA_Task* t );
void      FLA_queue_sort_task_array( void );
void      FLA_queue_flush( void );
void      FLA_queue_create_task_array( void );
void      FLA_queue_free_task_array( void );
int       FLA_queue_task_cost_compare( const void* t0, const void* t1 );
void      FLA_queue_print_costs( void );

FLA_Task* FLA_task_alloc_init( void* func, double cost, int n_int_args, int n_fla_args );
void      FLA_task_free( FLA_Task* t );


// -----------------------------------------------------------------------------



#define ENQUEUE_FLA_Gemm( transA, transB, alpha, A, B, beta, C ) \
        FLA_Queue_push( (void*)FLA_Gemm_external, \
                        2.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            ( transA == FLA_TRANSPOSE ? FLA_Obj_length( A ) \
                                                      :  FLA_Obj_width( A ) ), \
                        2, 5, \
                        transA, transB, \
                        alpha, A, B, beta, C )

#define ENQUEUE_FLA_Symm( side, uplo, alpha, A, B, beta, C ) \
        FLA_Queue_push( (void*)FLA_Symm_external, \
                        2.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            FLA_Obj_length( A ), \
                        2, 5, \
                        side, uplo, \
                        alpha, A, B, beta, C )

#define ENQUEUE_FLA_Syrk( uplo, transA, alpha, A, beta, C ) \
        FLA_Queue_push( (void*)FLA_Syrk_external, \
                        1.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            ( transA == FLA_TRANSPOSE ? FLA_Obj_length( A ) \
                                                      :  FLA_Obj_width( A ) ), \
                        2, 4, \
                        uplo, transA, \
                        alpha, A, beta, C )

#define ENQUEUE_FLA_Syr2k( uplo, transA, alpha, A, B, beta, C ) \
        FLA_Queue_push( (void*)FLA_Syr2k_external, \
                        2.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            ( transA == FLA_TRANSPOSE ? FLA_Obj_length( A ) \
                                                      :  FLA_Obj_width( A ) ), \
                        2, 5, \
                        uplo, transA, \
                        alpha, A, B, beta, C )

#define ENQUEUE_FLA_Trmm( side, uplo, trans, diag, alpha, A, C ) \
        FLA_Queue_push( (void*)FLA_Trmm_external, \
                        1.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            FLA_Obj_length( A ), \
                        4, 3, \
                        side, uplo, trans, diag, \
                        alpha, A, C )

#define ENQUEUE_FLA_Trsm( side, uplo, trans, diag, alpha, A, C ) \
        FLA_Queue_push( (void*)FLA_Trsm_external, \
                        1.0 * FLA_Obj_length( C ) * FLA_Obj_width( C ) * \
                            FLA_Obj_length( A ), \
                        4, 3, \
                        side, uplo, trans, diag, \
                        alpha, A, C )

