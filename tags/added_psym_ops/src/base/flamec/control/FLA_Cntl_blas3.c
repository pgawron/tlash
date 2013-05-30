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

//
// Level-3 BLAS
//

fla_gemm_t* FLA_Cntl_gemm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_gemm_t*      sub_gemm )
{
	fla_gemm_t* cntl;
	
	cntl = ( fla_gemm_t* ) FLA_malloc( sizeof(fla_gemm_t) );
	
	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_hemm_t* FLA_Cntl_hemm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_hemm_t*      sub_hemm,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2 )
{
	fla_hemm_t* cntl;
	
	cntl = ( fla_hemm_t* ) FLA_malloc( sizeof(fla_hemm_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_hemm    = sub_hemm;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;

	return cntl;
}

fla_herk_t* FLA_Cntl_herk_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scalr_t*     sub_scalr,
                                      fla_herk_t*      sub_herk,
                                      fla_gemm_t*      sub_gemm )
{
	fla_herk_t* cntl;
	
	cntl = ( fla_herk_t* ) FLA_malloc( sizeof(fla_herk_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_herk    = sub_herk;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_her2k_t* FLA_Cntl_her2k_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize, 
                                        fla_scalr_t*     sub_scalr,
                                        fla_her2k_t*     sub_her2k,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2 )
{
	fla_her2k_t* cntl;
	
	cntl = ( fla_her2k_t* ) FLA_malloc( sizeof(fla_her2k_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_her2k   = sub_her2k;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;

	return cntl;
}

fla_symm_t* FLA_Cntl_symm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_symm_t*      sub_symm,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2 )
{
	fla_symm_t* cntl;
	
	cntl = ( fla_symm_t* ) FLA_malloc( sizeof(fla_symm_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_symm    = sub_symm;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;

	return cntl;
}

fla_syrk_t* FLA_Cntl_syrk_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scalr_t*     sub_scalr,
                                      fla_syrk_t*      sub_syrk,
                                      fla_gemm_t*      sub_gemm )
{
	fla_syrk_t* cntl;
	
	cntl = ( fla_syrk_t* ) FLA_malloc( sizeof(fla_syrk_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_syrk    = sub_syrk;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_syr2k_t* FLA_Cntl_syr2k_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_syr2k_t*     sub_syr2k,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2 )
{
	fla_syr2k_t* cntl;
	
	cntl = ( fla_syr2k_t* ) FLA_malloc( sizeof(fla_syr2k_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_syr2k   = sub_syr2k;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;

	return cntl;
}

fla_trmm_t* FLA_Cntl_trmm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_trmm_t*      sub_trmm,
                                      fla_gemm_t*      sub_gemm )
{
	fla_trmm_t* cntl;
	
	cntl = ( fla_trmm_t* ) FLA_malloc( sizeof(fla_trmm_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_trmm    = sub_trmm;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_trsm_t* FLA_Cntl_trsm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_trsm_t*      sub_trsm,
                                      fla_gemm_t*      sub_gemm )
{
	fla_trsm_t* cntl;
	
	cntl = ( fla_trsm_t* ) FLA_malloc( sizeof(fla_trsm_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_trsm    = sub_trsm;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

