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


typedef struct FLA_Queue_s FLA_Queue;
typedef struct FLA_Task_s FLA_Task;

struct FLA_Queue_s
{
	// Number of tasks currently in queue
	int             n_tasks;

	// Pointers to head (front) and tail (back) of queue
	FLA_Task*       head;
	FLA_Task*       tail;
	
	// pointer to pthreads that will execute on the queue.
	pthread_t*      thread;
	
	// pthread attribute stores information about the pthread objects
	// that is required when pthread_create() is invoked.
	pthread_attr_t  thread_attr;
	
	// pthread mutex provides mutual exclusion removing items from the queue.
	pthread_mutex_t lock;
	
};

struct FLA_Task_s
{
	// FLOP cost of task
	int       cost;
	
	// Function pointer
	void*     func;
	
	// Function integer arguments
	int       n_int_args;
	int*      int_arg;
	
	// Function FLA_Obj arguments
	int       n_fla_args;
	FLA_Obj*  fla_arg;
	
	// Support for linked list of FLA_Tasks
	FLA_Task* next_task;

};


// -----------------------------------------------------------------------------


void      FLA_Queue_init();
void      FLA_Queue_exec();
void      FLA_Queue_finalize();

void      FLA_Queue_set_num_threads( int n_threads );
int       FLA_Queue_get_num_threads();

void      FLA_Queue_push( void* func, int cost, int n_int_params, int n_fla_params, 
                          int param0, int param1,
                          FLA_Obj param2, FLA_Obj param3, FLA_Obj param4, 
                          FLA_Obj param5, FLA_Obj param6 );

void      FLA_queue_exec_sync();

void      FLA_queue_push( FLA_Task* t );
FLA_Task* FLA_queue_pop();
void      FLA_queue_flush();

void*     FLA_queue_exec_thread( void* p );
void      FLA_queue_engage_tasks();

void      FLA_queue_threads_awaken();
void      FLA_queue_threads_sleep_init();
void      FLA_queue_threads_sleep();
void      FLA_queue_threads_sleep_barrier();

void      FLA_queue_exec_task( FLA_Task* t );

FLA_Task* FLA_task_alloc_init( void* func, int cost, int n_int_args, int n_fla_args );
void      FLA_task_free( FLA_Task* t );


// -----------------------------------------------------------------------------


#define ENQUEUE_FLA_Gemm( transA, transB, alpha, A, B, beta, C ) \
        FLA_Queue_push( (void*)FLA_Gemm, \
                        2 * FLA_Obj_length( A ) * FLA_Obj_width( A ) * FLA_Obj_width( B ), \
                        2, 5, \
                        transA, transB, \
                        alpha, A, B, beta, C )

