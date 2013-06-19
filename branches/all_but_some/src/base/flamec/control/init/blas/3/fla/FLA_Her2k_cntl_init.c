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

extern fla_scalr_t* fla_scalr_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;

fla_her2k_t*        fla_her2k_cntl_blas;
fla_her2k_t*        fla_her2k_cntl_ip;
fla_her2k_t*        fla_her2k_cntl_op;
fla_her2k_t*        fla_her2k_cntl_mm;
fla_blocksize_t*    fla_her2k_var3_bsize;
fla_blocksize_t*    fla_her2k_var9_bsize;

void FLA_Her2k_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_her2k_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_her2k_var9_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A and B are b x b blocks.
	fla_her2k_cntl_blas  = FLA_Cntl_her2k_obj_create( FLA_FLAT,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );
 
	// Create a control tree that assumes A and B form an inner panel product.
	fla_her2k_cntl_ip    = FLA_Cntl_her2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  fla_her2k_var9_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_her2k_cntl_blas,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A and B form an outer panel product.
	fla_her2k_cntl_op    = FLA_Cntl_her2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  fla_her2k_var3_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_her2k_cntl_blas,
	                                                  fla_gemm_cntl_blas,
	                                                  fla_gemm_cntl_blas );

	// Create a control tree that assumes A and B are both large.
	fla_her2k_cntl_mm    = FLA_Cntl_her2k_obj_create( FLA_FLAT,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  fla_her2k_var9_bsize,
	                                                  fla_scalr_cntl_blas,
	                                                  fla_her2k_cntl_op,
	                                                  NULL,
	                                                  NULL );
}

void FLA_Her2k_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_her2k_cntl_blas );

	FLA_Cntl_obj_free( fla_her2k_cntl_ip );
	FLA_Cntl_obj_free( fla_her2k_cntl_op );
	FLA_Cntl_obj_free( fla_her2k_cntl_mm );

	FLA_Blocksize_free( fla_her2k_var3_bsize );
	FLA_Blocksize_free( fla_her2k_var9_bsize );
}

