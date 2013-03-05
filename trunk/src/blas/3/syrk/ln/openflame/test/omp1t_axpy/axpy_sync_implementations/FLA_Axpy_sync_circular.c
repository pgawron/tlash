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
FLA_Error FLA_Axpy_sync_circular( FLA_Obj alpha, FLA_Obj X, FLA_Obj B )
{
  FLA_Obj XL,    XR,       X0,  X1,  X2;
  FLA_Obj BL,    BR,       B0,  B1,  B2;

  int n_stages    = FLA_omp_get_num_stages();
  int stage_width = FLA_omp_compute_stage_width( X );
  int thread_num  = omp_get_thread_num();
  int n_done      = 0;
  int b, i;

  // Start thread i on the ith panel partition of B.
  FLA_Part_1x2( X,    &XL,  &XR,    stage_width*thread_num, FLA_LEFT );
  FLA_Part_1x2( B,    &BL,  &BR,    stage_width*thread_num, FLA_LEFT );

  while ( n_done++ < n_stages ){

    // The last lockable partition may be smaller than the others.
    b = min( FLA_Obj_width( XR ), stage_width );
    
    FLA_Repart_1x2_to_1x3( XL,  /**/ XR,        &X0, /**/ &X1, &X2,
                           b, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    // Get the index of the current partition.
    i = FLA_Obj_width(XL)/stage_width;

    // Acquire lock[i] (the lock for X1 and B1).
    omp_set_lock( &fla_omp_lock[i] );

    // B1 := alpha * X1 + B1
    FLA_Axpy_external( alpha, X1, B1 );

    // Release lock[i] (the lock for X1 and B1).
    omp_unset_lock( &fla_omp_lock[i] );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &XL,  /**/ &XR,        X0, X1, /**/ X2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );

    // If this thread reaches the last partition, wrap back around to 
    // the first partition for the next iteration.
    if( FLA_Obj_width( XL ) == FLA_Obj_width( X ) )
    {
      FLA_Part_1x2( X,    &XL,  &XR,      0, FLA_LEFT );
      FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );
    }

  }

  return FLA_SUCCESS;
}
#endif
