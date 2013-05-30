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
// Level-1 BLAS
//

fla_axpy_t* FLA_Cntl_axpy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_axpy_t*      sub_axpy )
{
	fla_axpy_t* cntl;
	
	cntl = ( fla_axpy_t* ) FLA_malloc( sizeof(fla_axpy_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_axpy    = sub_axpy;

	return cntl;
}

fla_axpyt_t* FLA_Cntl_axpyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_axpyt_t*     sub_axpyt )
{
	fla_axpyt_t* cntl;
	
	cntl = ( fla_axpyt_t* ) FLA_malloc( sizeof(fla_axpyt_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_axpyt   = sub_axpyt;

	return cntl;
}

fla_copy_t* FLA_Cntl_copy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_copy_t*      sub_copy )
{
	fla_copy_t* cntl;
	
	cntl = ( fla_copy_t* ) FLA_malloc( sizeof(fla_copy_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copy    = sub_copy;

	return cntl;
}

fla_copyt_t* FLA_Cntl_copyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyt_t*     sub_copyt )
{
	fla_copyt_t* cntl;
	
	cntl = ( fla_copyt_t* ) FLA_malloc( sizeof(fla_copyt_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copyt   = sub_copyt;

	return cntl;
}

fla_copyr_t* FLA_Cntl_copyr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyr_t*     sub_copyr,
                                        fla_copy_t*      sub_copy )
{
	fla_copyr_t* cntl;
	
	cntl = ( fla_copyr_t* ) FLA_malloc( sizeof(fla_copyr_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copyr   = sub_copyr;
	cntl->sub_copy    = sub_copy;

	return cntl;
}

fla_scal_t* FLA_Cntl_scal_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal )
{
	fla_scal_t* cntl;
	
	cntl = ( fla_scal_t* ) FLA_malloc( sizeof(fla_scal_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;

	return cntl;
}

fla_scalr_t* FLA_Cntl_scalr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_scal_t*      sub_scal )
{
	fla_scalr_t* cntl;
	
	cntl = ( fla_scalr_t* ) FLA_malloc( sizeof(fla_scalr_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_scal    = sub_scal;

	return cntl;
}

fla_swap_t* FLA_Cntl_swap_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_swap_t*      sub_swap )
{
	fla_swap_t* cntl;
	
	cntl = ( fla_swap_t* ) FLA_malloc( sizeof(fla_swap_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_swap    = sub_swap;

	return cntl;
}

fla_tpose_t* FLA_Cntl_tpose_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_tpose_t*     sub_trans,
                                        fla_swap_t*      sub_swap )
{
	fla_tpose_t* cntl;
	
	cntl = ( fla_tpose_t* ) FLA_malloc( sizeof(fla_tpose_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_trans   = sub_trans;
	cntl->sub_swap    = sub_swap;

	return cntl;
}

