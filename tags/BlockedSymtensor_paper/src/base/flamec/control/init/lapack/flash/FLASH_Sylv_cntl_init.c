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

extern fla_gemm_t* flash_gemm_cntl_pm_bp;
extern fla_gemm_t* flash_gemm_cntl_ip_bb;

fla_sylv_t*        flash_sylv_cntl_leaf;
fla_sylv_t*        flash_sylv_cntl_mb;
fla_sylv_t*        flash_sylv_cntl;
fla_blocksize_t*   flash_sylv_bsize;

void FLASH_Sylv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_sylv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are b x b blocks.
	flash_sylv_cntl_leaf   = FLA_Cntl_sylv_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
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

	// Create a control tree that assumes A is a matrix and B is a block.
	flash_sylv_cntl_mb     = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT17,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_leaf,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_ip_bb,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that assumes A is a matrix and B is a matrix.
	flash_sylv_cntl        = FLA_Cntl_sylv_obj_create( FLA_HIER, 
	                                                   FLA_BLOCKED_VARIANT15,
	                                                   flash_sylv_bsize,
	                                                   flash_sylv_cntl_mb,
	                                                   NULL,
	                                                   NULL,
	                                                   flash_gemm_cntl_pm_bp,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );
}

void FLASH_Sylv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_sylv_cntl_leaf );
	FLA_Cntl_obj_free( flash_sylv_cntl_mb );
	FLA_Cntl_obj_free( flash_sylv_cntl );

	FLA_Blocksize_free( flash_sylv_bsize );
}

