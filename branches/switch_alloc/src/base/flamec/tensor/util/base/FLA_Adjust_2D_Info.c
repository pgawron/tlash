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

//Adjusts tensor object to fill in for 2D matrix object
FLA_Error FLA_Adjust_2D_info( FLA_Obj *A )
{
	dim_t i;
	dim_t  order = A->order;
    dim_t* size_inner = A->size_inner;
    dim_t* size = A->size;
    dim_t* offset = A->offset;
	dim_t* base_index = A->base->index;
    dim_t* base_stride = A->base->stride;
    dim_t* base_size_inner = A->base->size_inner;
    dim_t* base_size = A->base->size;
    
	A->base->m = base_size[0];
	A->base->n = base_size[1];
	A->base->rs = base_stride[0];
	A->base->cs = (A->base->rs == 1) ? (A->base->rs * A->base->m) : 1;
	A->base->m_inner = base_size_inner[0];
	A->base->n_inner = base_size_inner[1];
	A->base->m_index = base_index[0];
	A->base->n_index = base_index[1];

	A->offm = offset[0];
	A->offn = offset[1];
	A->m = size[0];
	A->n = size[1];
	A->m_inner = size_inner[0];
	A->n_inner = size_inner[1];

	for(i = 2; i < order; i++){
		A->base->n *= base_size[i];
		A->base->n_inner *= base_size_inner[i];
		//A->base->n_index *= base_index[i];
		//A->offn *= offset[i];
		A->n *= size[i];
		A->n_inner *= size_inner[i];
	}

  return FLA_SUCCESS;
}

