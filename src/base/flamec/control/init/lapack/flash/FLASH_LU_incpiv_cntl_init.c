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

extern fla_gemm_t*  flash_gemm_cntl_bp_bb;
extern fla_trsm_t*  flash_trsm_cntl_bp;
extern fla_appiv_t* flash_appiv_cntl_bp;

fla_lu_t*           flash_lu_incpiv_cntl_leaf;
fla_lu_t*           flash_lu_incpiv_cntl;
fla_blocksize_t*    flash_lu_incpiv_bsize;

void FLASH_LU_incpiv_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_lu_incpiv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_lu_incpiv_cntl_leaf   = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                      FLA_SUBPROBLEM,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL,
	                                                      NULL );

	// Create a control tree that assumes A is large.
	flash_lu_incpiv_cntl        = FLA_Cntl_lu_obj_create( FLA_HIER,
	                                                      FLA_BLOCKED_VARIANT1,
	                                                      flash_lu_incpiv_bsize,
	                                                      flash_lu_incpiv_cntl_leaf,
	                                                      flash_gemm_cntl_bp_bb,
	                                                      NULL,
	                                                      NULL,
	                                                      flash_trsm_cntl_bp,
	                                                      flash_trsm_cntl_bp,
	                                                      flash_appiv_cntl_bp,
	                                                      NULL );
}

void FLASH_LU_incpiv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_lu_incpiv_cntl_leaf );
	FLA_Cntl_obj_free( flash_lu_incpiv_cntl );

	FLA_Blocksize_free( flash_lu_incpiv_bsize );
}

