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

fla_swap_t*      fla_swap_cntl_panel;
fla_swap_t*      fla_swap_cntl_blas;

fla_tpose_t*     fla_tpose_cntl;
fla_tpose_t*     fla_tpose_cntl_unb;
fla_blocksize_t* fla_tpose_bsize;
fla_blocksize_t* fla_tpose_swap_bsize;

void FLA_Transpose_cntl_init()
{
	// Set blocksizes based on libgoto query.
	fla_tpose_bsize      = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_tpose_swap_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that performs unblocked variant 2 transposition.
	fla_tpose_cntl_unb   = FLA_Cntl_tpose_obj_create( FLA_FLAT, 
	                                                  FLA_UNBLOCKED_VARIANT2,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that invokes an external implementation of swap.
	fla_swap_cntl_blas   = FLA_Cntl_swap_obj_create( FLA_FLAT,
	                                                 FLA_SUBPROBLEM,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree that invokes unblocked variant 2 of swap.
	fla_swap_cntl_panel  = FLA_Cntl_swap_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT2, 
	                                                 fla_tpose_swap_bsize,
	                                                 fla_swap_cntl_blas );

	// Create a control tree that assumes a large matrix argument.
	fla_tpose_cntl       = FLA_Cntl_tpose_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT2, 
	                                                  fla_tpose_bsize,
	                                                  fla_tpose_cntl_unb,
	                                                  fla_swap_cntl_panel );
}

void FLA_Transpose_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_tpose_cntl );
	FLA_Cntl_obj_free( fla_tpose_cntl_unb );
	FLA_Cntl_obj_free( fla_swap_cntl_panel );
	FLA_Cntl_obj_free( fla_swap_cntl_blas );

	FLA_Blocksize_free( fla_tpose_bsize );
	FLA_Blocksize_free( fla_tpose_swap_bsize );
}

