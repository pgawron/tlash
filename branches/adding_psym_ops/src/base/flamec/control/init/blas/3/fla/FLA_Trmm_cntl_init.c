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

extern fla_scal_t* fla_scal_cntl_blas;
extern fla_gemm_t* fla_gemm_cntl_blas;

fla_trmm_t*        fla_trmm_cntl_blas;
fla_trmm_t*        fla_trmm_cntl_bp;
fla_trmm_t*        fla_trmm_cntl_mp;
fla_trmm_t*        fla_trmm_cntl_mm;
fla_blocksize_t*   fla_trmm_var1_bsize;
fla_blocksize_t*   fla_trmm_var3_bsize;

void FLA_Trmm_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_trmm_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_trmm_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree that assumes A and B are b x b blocks.
	fla_trmm_cntl_blas  = FLA_Cntl_trmm_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	fla_trmm_cntl_bp    = FLA_Cntl_trmm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_trmm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trmm_cntl_blas,
	                                                NULL );

	// Create a control tree that assumes A is large and B is a panel.
	fla_trmm_cntl_mp    = FLA_Cntl_trmm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT1,
	                                                fla_trmm_var1_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trmm_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create a control tree that assumes A and B are both large.
	fla_trmm_cntl_mm    = FLA_Cntl_trmm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_trmm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_trmm_cntl_mp,
	                                                NULL );
}

void FLA_Trmm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_trmm_cntl_blas );

	FLA_Cntl_obj_free( fla_trmm_cntl_bp );
	FLA_Cntl_obj_free( fla_trmm_cntl_mp );
	FLA_Cntl_obj_free( fla_trmm_cntl_mm );

	FLA_Blocksize_free( fla_trmm_var1_bsize );
	FLA_Blocksize_free( fla_trmm_var3_bsize );
}

