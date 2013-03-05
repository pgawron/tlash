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


#ifdef FLA_ENABLE_MULTITHREADING


#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
#ifdef FLA_ENABLE_TIDSP
#include <ti/omp/omp.h>
#else
#include <omp.h>
#endif
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
#include <pthread.h>
#endif


void FLA_Lock_init( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_init

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_init_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_init( &(fla_lock_ptr->lock), NULL );
#endif
}


void FLA_Lock_acquire( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_acquire

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_set_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_lock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_Lock_release( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_release

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_unset_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_unlock( &(fla_lock_ptr->lock) );
#endif
}


void FLA_Lock_destroy( FLA_Lock* fla_lock_ptr )
/*----------------------------------------------------------------------------

   FLA_Lock_destroy

----------------------------------------------------------------------------*/
{
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_destroy_lock( &(fla_lock_ptr->lock) );
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_destroy( &(fla_lock_ptr->lock) );
#endif
}


#endif // FLA_ENABLE_MULTITHREADING

