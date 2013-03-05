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

extern fla_scal_t* flash_scal_cntl;

fla_gemm_t*      flash_gemm_cntl_blas;
fla_gemm_t*      flash_gemm_cntl_mm_mp;
fla_gemm_t*      flash_gemm_cntl_mm_pm;
fla_gemm_t*      flash_gemm_cntl_mm_op;
fla_gemm_t*      flash_gemm_cntl_mp_pb;
fla_gemm_t*      flash_gemm_cntl_mp_ip;
fla_gemm_t*      flash_gemm_cntl_pm_bp;
fla_gemm_t*      flash_gemm_cntl_pm_ip;
fla_gemm_t*      flash_gemm_cntl_op_bp;
fla_gemm_t*      flash_gemm_cntl_op_pb;
fla_gemm_t*      flash_gemm_cntl_pb_bb;
fla_gemm_t*      flash_gemm_cntl_bp_bb;
fla_gemm_t*      flash_gemm_cntl_ip_bb;

fla_gemm_t*      flash_gemm_cntl_mm;
fla_gemm_t*      flash_gemm_cntl_mp;
fla_gemm_t*      flash_gemm_cntl_pm;
fla_gemm_t*      flash_gemm_cntl_op;
fla_gemm_t*      flash_gemm_cntl_pb;
fla_gemm_t*      flash_gemm_cntl_bp;
fla_gemm_t*      flash_gemm_cntl_ip;

fla_blocksize_t* flash_gemm_bsize;

void FLASH_Gemm_cntl_init()
{
	// Set gemm blocksize for hierarchical storage.
	flash_gemm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree node that executes a gemm subproblem.
	flash_gemm_cntl_blas  = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create control trees for situations where one dimension is large.
	flash_gemm_cntl_pb_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );
	flash_gemm_cntl_bp_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );
	flash_gemm_cntl_ip_bb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_blas );

	// Create control trees for situations where two dimensions are large.
	flash_gemm_cntl_mp_ip = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_ip_bb );
	flash_gemm_cntl_op_bp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_bp_bb );
	flash_gemm_cntl_pm_ip = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_ip_bb );
	flash_gemm_cntl_op_pb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pb_bb );
	flash_gemm_cntl_mp_pb = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pb_bb );
	flash_gemm_cntl_pm_bp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_bp_bb );

	// Create control trees for situations where all dimensions are large.
	flash_gemm_cntl_mm_pm = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_pm_ip );
	flash_gemm_cntl_mm_mp = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT3,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_mp_ip );
	flash_gemm_cntl_mm_op = FLA_Cntl_gemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemm_cntl_op_bp );

	// Alias select control trees for convenience, when the caller doesn't
	// care (as he usually doesn't when partitioning hierarchical matricies)
	// which order the matrix is partitioned into blocks
	flash_gemm_cntl_mm = flash_gemm_cntl_mm_op;
	flash_gemm_cntl_mp = flash_gemm_cntl_mp_pb;
	flash_gemm_cntl_pm = flash_gemm_cntl_pm_bp;
	flash_gemm_cntl_op = flash_gemm_cntl_op_pb;
	flash_gemm_cntl_pb = flash_gemm_cntl_pb_bb;
	flash_gemm_cntl_bp = flash_gemm_cntl_bp_bb;
	flash_gemm_cntl_ip = flash_gemm_cntl_ip_bb;
	
}

void FLASH_Gemm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_gemm_cntl_blas );

	FLA_Cntl_obj_free( flash_gemm_cntl_pb_bb );
	FLA_Cntl_obj_free( flash_gemm_cntl_bp_bb );
	FLA_Cntl_obj_free( flash_gemm_cntl_ip_bb );

	FLA_Cntl_obj_free( flash_gemm_cntl_mp_ip );
	FLA_Cntl_obj_free( flash_gemm_cntl_op_bp );
	FLA_Cntl_obj_free( flash_gemm_cntl_pm_ip );
	FLA_Cntl_obj_free( flash_gemm_cntl_op_pb );
	FLA_Cntl_obj_free( flash_gemm_cntl_mp_pb );
	FLA_Cntl_obj_free( flash_gemm_cntl_pm_bp );

	FLA_Cntl_obj_free( flash_gemm_cntl_mm_pm );
	FLA_Cntl_obj_free( flash_gemm_cntl_mm_mp );
	FLA_Cntl_obj_free( flash_gemm_cntl_mm_op );

	FLA_Blocksize_free( flash_gemm_bsize );
}

