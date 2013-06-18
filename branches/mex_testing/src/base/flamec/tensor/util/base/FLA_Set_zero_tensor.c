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

FLA_Error FLA_Set_zero_tensor( FLA_Obj A )
{
	dim_t i;
	dim_t nElem = FLA_Obj_num_elem_alloc(A);

	//If it is blocked, set each block
	if(FLA_Obj_elemtype(A) == FLA_MATRIX || FLA_Obj_elemtype(A) == FLA_TENSOR){
		FLA_Obj* buffer = (FLA_Obj*)FLA_Obj_base_buffer(A);
		for(i = 0; i < nElem; i++)
			FLA_Set_zero_tensor(buffer[i]);
	}
	//If it is a block, set it (only if stored)
	else{
		if(A.isStored){
			memset(&(((double*)FLA_Obj_base_buffer(A))[0]), 0, nElem * sizeof(double));
		}
	}

  return FLA_SUCCESS;
}

