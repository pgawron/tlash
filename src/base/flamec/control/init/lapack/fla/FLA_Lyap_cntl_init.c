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

extern fla_scal_t*  fla_scal_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_sylv_t*  fla_sylv_cntl;

fla_lyap_t*         fla_lyap_cntl_leaf;
fla_lyap_t*         fla_lyap_cntl;
fla_blocksize_t*    fla_lyap_bsize;

void FLA_Lyap_cntl_init()
{
	// Set blocksize with default value for conventional storage.
	fla_lyap_bsize       = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Create a control tree to invoke an unblocked variant.
	fla_lyap_cntl_leaf   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                                 FLA_UNBLOCKED_VARIANT1,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );

	// Create a control tree to invoke a blocked variant.
	fla_lyap_cntl        = FLA_Cntl_lyap_obj_create( FLA_FLAT, 
	                                                 FLA_BLOCKED_VARIANT1,
	                                                 fla_lyap_bsize,
	                                                 fla_scal_cntl_blas,
	                                                 fla_lyap_cntl_leaf,
	                                                 fla_sylv_cntl,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_hemm_cntl_blas,
	                                                 fla_her2k_cntl_blas );
}

void FLA_Lyap_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_lyap_cntl_leaf );
	FLA_Cntl_obj_free( fla_lyap_cntl );

	FLA_Blocksize_free( fla_lyap_bsize );
}

