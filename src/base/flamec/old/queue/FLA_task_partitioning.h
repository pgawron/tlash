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

void FLA_Task_partitioning_init( void );
void FLA_Task_partitioning_set( int n_threads, int tag0_val, int tag1_val, int tag2_val, int tag3_val, int tag4_val );

int  FLA_Task_compute_blocksize( int tag, FLA_Obj A, FLA_Obj A_proc, FLA_Quadrant from );

int  FLA_task_get_num_partitions( int n_threads, int tag );

int  FLA_task_determine_matrix_size( FLA_Obj A, FLA_Quadrant from );
int  FLA_task_determine_relative_blocksize( int A_size, int A_proc_size, int n_part );
int  FLA_task_determine_absolute_blocksize( int A_size, int A_proc_size, int nb_alg );

