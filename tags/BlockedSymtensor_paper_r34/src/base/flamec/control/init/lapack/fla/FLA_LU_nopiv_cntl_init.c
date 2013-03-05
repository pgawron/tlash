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

extern fla_gemm_t* fla_gemm_cntl_blas;
extern fla_trsm_t* fla_trsm_cntl_blas;

fla_lu_t*          fla_lu_nopiv_cntl;
fla_lu_t*          fla_lu_nopiv_cntl2;

fla_lu_t*          fla_lu_nopiv_cntl_in;
fla_lu_t*          fla_lu_nopiv_cntl_leaf;
fla_blocksize_t*   fla_lu_nopiv_var5_bsize;
fla_blocksize_t*   fla_lu_nopiv_var5_bsize_in;
double             fla_lu_nopiv_var5_in_to_ou_bsize_ratio = 0.25;

void FLA_LU_nopiv_cntl_init()
{
	// Set blocksizes with default values for conventional storage.
	fla_lu_nopiv_var5_bsize    = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	fla_lu_nopiv_var5_bsize_in = FLA_Blocksize_create_copy( fla_lu_nopiv_var5_bsize );
	FLA_Blocksize_scale( fla_lu_nopiv_var5_bsize_in, fla_lu_nopiv_var5_in_to_ou_bsize_ratio );

	// Create a control tree to invoke unblocked variant 1.
	fla_lu_nopiv_cntl_leaf   = FLA_Cntl_lu_obj_create( FLA_FLAT,
	                                                   FLA_UNB_OPT_VARIANT5,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree for small subproblems.
	fla_lu_nopiv_cntl_in     = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_lu_nopiv_var5_bsize_in,
	                                                   fla_lu_nopiv_cntl_leaf,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree for larger problems with one level of recursion.
	fla_lu_nopiv_cntl2       = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_lu_nopiv_var5_bsize,
	                                                   fla_lu_nopiv_cntl_in,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree for large problems with no extra recursion.
	fla_lu_nopiv_cntl        = FLA_Cntl_lu_obj_create( FLA_FLAT, 
	                                                   FLA_BLOCKED_VARIANT5,
	                                                   fla_lu_nopiv_var5_bsize,
	                                                   fla_lu_nopiv_cntl_leaf,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_gemm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   fla_trsm_cntl_blas,
	                                                   NULL,
	                                                   NULL );
}

void FLA_LU_nopiv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_lu_nopiv_cntl );
	FLA_Cntl_obj_free( fla_lu_nopiv_cntl2 );
	FLA_Cntl_obj_free( fla_lu_nopiv_cntl_in );
	FLA_Cntl_obj_free( fla_lu_nopiv_cntl_leaf );

	FLA_Blocksize_free( fla_lu_nopiv_var5_bsize );
	FLA_Blocksize_free( fla_lu_nopiv_var5_bsize_in );
}

