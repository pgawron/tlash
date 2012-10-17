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

extern fla_scalr_t* flash_scalr_cntl;
extern fla_gemm_t*  flash_gemm_cntl_pb_bb;

fla_herk_t*         flash_herk_cntl_blas;
fla_herk_t*         flash_herk_cntl_ip;
fla_herk_t*         flash_herk_cntl_op;
fla_herk_t*         flash_herk_cntl_mm;
fla_blocksize_t*    flash_herk_bsize;

void FLASH_Herk_cntl_init()
{
	// Set herk blocksize for hierarchical storage.
	flash_herk_bsize      = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A is a b x b block.
	flash_herk_cntl_blas  = FLA_Cntl_herk_obj_create( FLA_HIER,
	                                                  FLA_SUBPROBLEM,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL,
	                                                  NULL );

	// Create a control tree that assumes A * A' forms an inner panel product.
	flash_herk_cntl_ip    = FLA_Cntl_herk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_herk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_herk_cntl_blas,
	                                                  NULL );

	// Create a control tree that assumes A * A' forms an outer panel product.
	flash_herk_cntl_op    = FLA_Cntl_herk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT2,
	                                                  flash_herk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_herk_cntl_blas,
	                                                  flash_gemm_cntl_pb_bb );

	// Create a control tree that assumes A is large.
	flash_herk_cntl_mm    = FLA_Cntl_herk_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT5,
	                                                  flash_herk_bsize,
	                                                  flash_scalr_cntl,
	                                                  flash_herk_cntl_op,
	                                                  NULL );
}

void FLASH_Herk_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_herk_cntl_blas );

	FLA_Cntl_obj_free( flash_herk_cntl_ip );
	FLA_Cntl_obj_free( flash_herk_cntl_op );
	FLA_Cntl_obj_free( flash_herk_cntl_mm );

	FLA_Blocksize_free( flash_herk_bsize );
}

