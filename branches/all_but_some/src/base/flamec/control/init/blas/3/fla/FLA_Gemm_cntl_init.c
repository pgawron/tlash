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

fla_gemm_t*      fla_gemm_cntl_blas;

fla_gemm_t*      fla_gemm_cntl_pb_bb;
fla_gemm_t*      fla_gemm_cntl_bp_bb;
fla_gemm_t*      fla_gemm_cntl_ip_bb;

fla_gemm_t*      fla_gemm_cntl_mp_ip;
fla_gemm_t*      fla_gemm_cntl_mp_ip_bb;
fla_gemm_t*      fla_gemm_cntl_op_bp;
fla_gemm_t*      fla_gemm_cntl_op_bp_bb;
fla_gemm_t*      fla_gemm_cntl_pm_ip;
fla_gemm_t*      fla_gemm_cntl_pm_ip_bb;
fla_gemm_t*      fla_gemm_cntl_op_pb;
fla_gemm_t*      fla_gemm_cntl_op_pb_bb;
fla_gemm_t*      fla_gemm_cntl_mp_pb;
fla_gemm_t*      fla_gemm_cntl_mp_pb_bb;
fla_gemm_t*      fla_gemm_cntl_pm_bp;
fla_gemm_t*      fla_gemm_cntl_pm_bp_bb;

fla_gemm_t*      fla_gemm_cntl_mm_pm;
fla_gemm_t*      fla_gemm_cntl_mm_pm_ip;
fla_gemm_t*      fla_gemm_cntl_mm_pm_ip_bb;
fla_gemm_t*      fla_gemm_cntl_mm_mp;
fla_gemm_t*      fla_gemm_cntl_mm_mp_ip;
fla_gemm_t*      fla_gemm_cntl_mm_mp_ip_bb;
fla_gemm_t*      fla_gemm_cntl_mm_op;
fla_gemm_t*      fla_gemm_cntl_mm_op_bp;
fla_gemm_t*      fla_gemm_cntl_mm_op_bp_bb;

fla_blocksize_t* fla_gemm_var1_bsize;
fla_blocksize_t* fla_gemm_var3_bsize;
fla_blocksize_t* fla_gemm_var5_bsize;

void FLA_Gemm_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_gemm_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_gemm_var3_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_gemm_var5_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree node that executes a gemm subproblem.
	fla_gemm_cntl_blas  = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_SUBPROBLEM,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create control trees for situations where one dimension is large.
	fla_gemm_cntl_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT1,
	                                                fla_gemm_var1_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );
	fla_gemm_cntl_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT3,
	                                                fla_gemm_var3_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );
	fla_gemm_cntl_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                FLA_BLOCKED_VARIANT5,
	                                                fla_gemm_var5_bsize,
	                                                fla_scal_cntl_blas,
	                                                fla_gemm_cntl_blas );

	// Create control trees for situations where two dimensions are large.
	fla_gemm_cntl_mp_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_mp_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_ip_bb );
	fla_gemm_cntl_op_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_op_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT1,
	                                                   fla_gemm_var1_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_bp_bb );
	fla_gemm_cntl_pm_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_pm_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_ip_bb );
	fla_gemm_cntl_op_pb    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_op_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   fla_gemm_var3_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_pb_bb );
	fla_gemm_cntl_mp_pb    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_mp_pb_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_pb_bb );
	fla_gemm_cntl_pm_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_blas );
	fla_gemm_cntl_pm_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_gemm_var5_bsize,
	                                                   fla_scal_cntl_blas,
	                                                   fla_gemm_cntl_bp_bb );

	// Create control trees for situations where all dimensions are large.
	fla_gemm_cntl_mm_pm       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_pm_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_pm_ip );
	fla_gemm_cntl_mm_pm_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      fla_gemm_var1_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_pm_ip_bb );
	fla_gemm_cntl_mm_mp       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_mp_ip    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_mp_ip );
	fla_gemm_cntl_mm_mp_ip_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT3,
	                                                      fla_gemm_var3_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_mp_ip_bb );
	fla_gemm_cntl_mm_op       = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_blas );
	fla_gemm_cntl_mm_op_bp    = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_op_bp );
	fla_gemm_cntl_mm_op_bp_bb = FLA_Cntl_gemm_obj_create( FLA_FLAT,
	                                                      FLA_BLOCKED_VARIANT5,
	                                                      fla_gemm_var5_bsize,
	                                                      fla_scal_cntl_blas,
	                                                      fla_gemm_cntl_op_bp_bb );
}

void FLA_Gemm_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_gemm_cntl_blas );

	FLA_Cntl_obj_free( fla_gemm_cntl_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_bp_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_ip_bb );

	FLA_Cntl_obj_free( fla_gemm_cntl_mp_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_bp_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_pb );
	FLA_Cntl_obj_free( fla_gemm_cntl_op_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_pb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mp_pb_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_pm_bp_bb );

	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_pm_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp_ip );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_mp_ip_bb );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op_bp );
	FLA_Cntl_obj_free( fla_gemm_cntl_mm_op_bp_bb );

	FLA_Blocksize_free( fla_gemm_var1_bsize );
	FLA_Blocksize_free( fla_gemm_var3_bsize );
	FLA_Blocksize_free( fla_gemm_var5_bsize );
}

