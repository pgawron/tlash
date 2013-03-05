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

#ifdef OPENMP
FLA_Error REF_Axpy_sync_pipeline( FLA_Obj alpha, FLA_Obj X, FLA_Obj B )
{
  double* x_buf, *b_buf;
  int     x_m, x_n, x_ldim, b_ldim;
  int     b, j, j2, j_part, nb_alg;
  int     i_one = 1;
  double  alpha_value;

  // Compute the width of one lockable partition.
  nb_alg    = FLA_omp_compute_stage_width( X );
  x_n       = FLA_Obj_width( X );
  x_m       = FLA_Obj_length( X );
  x_buf     = FLA_Obj_buffer_at_view( X );
  b_buf     = FLA_Obj_buffer_at_view( B );
  x_ldim    = FLA_Obj_ldim( X );
  b_ldim    = FLA_Obj_ldim( B );
  alpha_value = FLA_DOUBLE_VALUE( alpha );

  for( j = 0; j < x_n; j += nb_alg )
  {
    b = min( x_n-j, nb_alg );

    /*------------------------------------------------------------*/

    // Get the index of the current partition.
    j_part = j/nb_alg;

    // Acquire lock[j_part] (the lock for X1 and B1).
    omp_set_lock( &fla_omp_lock[j_part] );

    // B1 := alpha * X1 + B1
    for( j2 = 0; j2 < b; ++j2 )
      FLA_C2F( daxpy) ( &x_m,
                        &alpha_value,
                        x_buf + (j+j2)*x_ldim, &i_one,
                        b_buf + (j+j2)*b_ldim, &i_one );

    // Release lock[j_part] (the lock for X1 and B1).
    omp_unset_lock( &fla_omp_lock[j_part] );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
#endif
