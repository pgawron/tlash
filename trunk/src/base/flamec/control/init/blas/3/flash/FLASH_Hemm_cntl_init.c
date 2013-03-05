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
extern fla_gemm_t* flash_gemm_cntl_op_bp;
extern fla_gemm_t* flash_gemm_cntl_mm_pm;
extern fla_gemm_t* flash_gemm_cntl_mm_mp;

fla_hemm_t*        flash_hemm_cntl_blas;
fla_hemm_t*        flash_hemm_cntl_bp;
fla_hemm_t*        flash_hemm_cntl_mp;
fla_hemm_t*        flash_hemm_cntl_mm;
fla_blocksize_t*   flash_hemm_bsize;

void FLASH_Hemm_cntl_init()
{
	// Set hemm blocksize for hierarchical storage.
	flash_hemm_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_hemm_cntl_blas  = FLA_Cntl_hemm_obj_create( FLA_HIER, 
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is a block and B is a panel.
	flash_hemm_cntl_bp    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_blas,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A is large and B is a panel.
	flash_hemm_cntl_mp    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_blas,
	                                                  flash_gemm_cntl_op_bp,
	                                                  flash_gemm_cntl_mm_mp );

	// Create a control tree that assumes A and B are both large.
	flash_hemm_cntl_mm    = FLA_Cntl_hemm_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT9,
	                                                  flash_hemm_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_hemm_cntl_mp,
	                                                  NULL,
	                                                  NULL );
}

void FLASH_Hemm_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_hemm_cntl_blas );

	FLA_Cntl_obj_free( flash_hemm_cntl_bp );
	FLA_Cntl_obj_free( flash_hemm_cntl_mp );
	FLA_Cntl_obj_free( flash_hemm_cntl_mm );

	FLA_Blocksize_free( flash_hemm_bsize );
}

