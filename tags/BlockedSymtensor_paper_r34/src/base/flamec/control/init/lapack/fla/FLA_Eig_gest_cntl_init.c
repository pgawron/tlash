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

extern fla_axpy_t*  fla_axpy_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;

fla_eig_gest_t*     fla_eig_gest_ix_cntl;
fla_eig_gest_t*     fla_eig_gest_nx_cntl;
fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf;
fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf;
fla_blocksize_t*    fla_eig_gest_var1_bsize;

void FLA_Eig_gest_cntl_init()
{
	// Set blocksize with default values for conventional storage.
	fla_eig_gest_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke the unblocked subproblem (inverse cases).
	fla_eig_gest_ix_cntl_leaf = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                          FLA_BLOCKED_EXTERN,
#else
	                                                          FLA_UNB_OPT_VARIANT3,
#endif
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL );

	// Create a control tree to invoke the unblocked subproblem (no inverse cases).
	fla_eig_gest_nx_cntl_leaf = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                          FLA_BLOCKED_EXTERN,
#else
	                                                          FLA_UNB_OPT_VARIANT2,
#endif
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL,
	                                                          NULL );

	// Create a control tree for large problems with no extra recursion
	// (inverse cases).
	fla_eig_gest_ix_cntl      = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
	                                                          FLA_BLOCKED_VARIANT4, 
	                                                          fla_eig_gest_var1_bsize,
	                                                          fla_eig_gest_ix_cntl_leaf,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_hemm_cntl_blas,
	                                                          fla_her2k_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trsm_cntl_blas,
	                                                          fla_trsm_cntl_blas );

	// Create a control tree for large problems with no extra recursion
	// (no inverse cases).
	fla_eig_gest_nx_cntl      = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
	                                                          FLA_BLOCKED_VARIANT2, 
	                                                          fla_eig_gest_var1_bsize,
	                                                          fla_eig_gest_nx_cntl_leaf,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_axpy_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_gemm_cntl_blas,
	                                                          fla_hemm_cntl_blas,
	                                                          fla_her2k_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trmm_cntl_blas,
	                                                          fla_trsm_cntl_blas,
	                                                          fla_trsm_cntl_blas );
}

void FLA_Eig_gest_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_eig_gest_ix_cntl );
	FLA_Cntl_obj_free( fla_eig_gest_nx_cntl );
	FLA_Cntl_obj_free( fla_eig_gest_ix_cntl_leaf );
	FLA_Cntl_obj_free( fla_eig_gest_nx_cntl_leaf );

	FLA_Blocksize_free( fla_eig_gest_var1_bsize );
}

