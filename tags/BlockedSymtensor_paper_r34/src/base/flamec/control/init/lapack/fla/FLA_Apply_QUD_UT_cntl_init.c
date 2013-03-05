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

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

fla_apqudut_t*      fla_apqudut_cntl_leaf;
fla_blocksize_t*    fla_apqudut_var1_bsize;

void FLA_Apply_QUD_UT_cntl_init()
{
	// Set the blocksize to the default value for conventional storage,
	// but scaled down.
	fla_apqudut_var1_bsize = FLA_Query_blocksizes( FLA_DIMENSION_MIN );
	FLA_Blocksize_scale( fla_apqudut_var1_bsize, FLA_UDDATE_INNER_TO_OUTER_B_RATIO );

	// Create a control tree to invoke variant 1.
	fla_apqudut_cntl_leaf = FLA_Cntl_apqudut_obj_create( FLA_FLAT,
	                                                     FLA_BLOCKED_VARIANT1, 
	                                                     fla_apqudut_var1_bsize,
	                                                     NULL,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_gemm_cntl_blas,
	                                                     fla_trsm_cntl_blas,
	                                                     fla_copyt_cntl_blas,
	                                                     fla_axpyt_cntl_blas );
}

void FLA_Apply_QUD_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_apqudut_cntl_leaf );

	FLA_Blocksize_free( fla_apqudut_var1_bsize );
}

