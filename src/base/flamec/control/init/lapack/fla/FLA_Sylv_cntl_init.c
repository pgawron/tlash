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

fla_sylv_t*        fla_sylv_cntl_leaf;
fla_sylv_t*        fla_sylv_cntl_mb;
fla_sylv_t*        fla_sylv_cntl;
fla_blocksize_t*   fla_sylv_bsize;

void FLA_Sylv_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_sylv_bsize       = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke LAPACK.
	fla_sylv_cntl_leaf   = FLA_Cntl_sylv_obj_create( FLA_FLAT,
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS
	                                                 FLA_BLOCKED_EXTERN, 
#else
	                                                 FLA_UNB_OPT_VARIANT1,
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
	                                                 NULL );

	// Create a control tree to invoke variant 15.
	fla_sylv_cntl_mb     = FLA_Cntl_sylv_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT15,
	                                                 fla_sylv_bsize,
	                                                 fla_sylv_cntl_leaf,
	                                                 NULL,
	                                                 NULL,
	                                                 fla_gemm_cntl_blas,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke variant 17.
	fla_sylv_cntl        = FLA_Cntl_sylv_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT17,
	                                                 fla_sylv_bsize,
	                                                 fla_sylv_cntl_mb,
	                                                 NULL,
	                                                 NULL,
	                                                 fla_gemm_cntl_blas,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );
}

void FLA_Sylv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_sylv_cntl_leaf );
	FLA_Cntl_obj_free( fla_sylv_cntl_mb );
	FLA_Cntl_obj_free( fla_sylv_cntl );

	FLA_Blocksize_free( fla_sylv_bsize );
}

