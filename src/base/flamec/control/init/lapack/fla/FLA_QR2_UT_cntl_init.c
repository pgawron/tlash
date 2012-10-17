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
extern fla_copy_t* fla_copy_cntl_blas;
extern fla_axpy_t* fla_axpy_cntl_blas;

fla_qr2ut_t*       fla_qr2ut_cntl_unb;
fla_qr2ut_t*       fla_qr2ut_cntl_leaf;
fla_blocksize_t*   fla_qr2ut_var1_bsize;

void FLA_QR2_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_qr2ut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_qr2ut_var1_bsize, FLA_QR_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke unblocked variant 1.
	fla_qr2ut_cntl_unb = FLA_Cntl_qr2ut_obj_create( FLA_FLAT,
	                                                FLA_UNB_OPT_VARIANT1, 
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL,
	                                                NULL );

	// Create a control tree for small-to-medium sequential problems and
	// as the means to compute on FLASH blocks.
	fla_qr2ut_cntl_leaf = FLA_Cntl_qr2ut_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT1, 
	                                                 fla_qr2ut_var1_bsize,
	                                                 fla_qr2ut_cntl_unb,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_copy_cntl_blas,
	                                                 fla_axpy_cntl_blas );

}

void FLA_QR2_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_qr2ut_cntl_unb );
	FLA_Cntl_obj_free( fla_qr2ut_cntl_leaf );

	FLA_Blocksize_free( fla_qr2ut_var1_bsize );
}

