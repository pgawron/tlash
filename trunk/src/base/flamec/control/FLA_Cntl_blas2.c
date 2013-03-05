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
// Level-2 BLAS
//

fla_gemv_t* FLA_Cntl_gemv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize, 
                                      fla_scal_t*      sub_scal,
                                      fla_gemv_t*      sub_gemv )
{
	fla_gemv_t* cntl;
	
	cntl = ( fla_gemv_t* ) FLA_malloc( sizeof(fla_gemv_t) );
	
	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_gemv    = sub_gemv;

	return cntl;
}

fla_trsv_t* FLA_Cntl_trsv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize, 
                                      fla_trsv_t*      sub_trsv,
                                      fla_gemv_t*      sub_gemv )
{
	fla_trsv_t* cntl;
	
	cntl = ( fla_trsv_t* ) FLA_malloc( sizeof(fla_trsv_t) );
	
	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_trsv    = sub_trsv;
	cntl->sub_gemv    = sub_gemv;

	return cntl;
}

