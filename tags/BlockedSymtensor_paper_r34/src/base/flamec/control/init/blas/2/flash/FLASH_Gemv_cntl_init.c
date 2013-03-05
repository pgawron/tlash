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

fla_gemv_t*      flash_gemv_cntl_blas;
fla_gemv_t*      flash_gemv_cntl_fm_rp;
fla_gemv_t*      flash_gemv_cntl_fm_cp;
fla_gemv_t*      flash_gemv_cntl_rp_bv;
fla_gemv_t*      flash_gemv_cntl_cp_bv;
fla_blocksize_t* flash_gemv_bsize;

void FLASH_Gemv_cntl_init()
{
	// Set gemv blocksize for hierarchical storage.
	flash_gemv_bsize       = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree node that executes a gemv subproblem.
	flash_gemv_cntl_blas   = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create control trees for situations where one dimension is large.
	flash_gemv_cntl_cp_bv = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_blas );
	flash_gemv_cntl_rp_bv = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_blas );

	// Create control trees for situations where both dimensions are large.
	flash_gemv_cntl_fm_rp = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_rp_bv );
	flash_gemv_cntl_fm_cp = FLA_Cntl_gemv_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_gemv_bsize,
	                                                  flash_scal_cntl,
	                                                  flash_gemv_cntl_cp_bv );
}

void FLASH_Gemv_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_gemv_cntl_blas );

	FLA_Cntl_obj_free( flash_gemv_cntl_cp_bv );
	FLA_Cntl_obj_free( flash_gemv_cntl_rp_bv );

	FLA_Cntl_obj_free( flash_gemv_cntl_fm_rp );
	FLA_Cntl_obj_free( flash_gemv_cntl_fm_cp );

	FLA_Blocksize_free( flash_gemv_bsize );
}

