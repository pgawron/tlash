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

#include "FLA_Permute.h"

FLA_Error FLA_Permute_single( FLA_Obj A, dim_t permutation[], FLA_Obj* B){
 
    FLA_Datatype datatype = FLA_Obj_datatype( A );	
	dim_t order;
	dim_t* stride_A;
	dim_t* size_A;
	dim_t* offset_A;
//	dim_t* base_index_A;

	order = FLA_Obj_order(A);
	stride_A = FLA_Obj_stride(A);
	offset_A = FLA_Obj_offset(A);
	size_A = FLA_Obj_size(A);

    FLA_Obj_create_tensor_without_buffer(datatype, order, size_A, B);
	FLA_Obj_attach_buffer_to_tensor((A.base)->buffer, order, stride_A, B);

	//Permute tensor info
	for(dim_t i = 0; i < order; i++){
		(B->offset)[i] = offset_A[permutation[i]];
		(B->size)[i] = size_A[permutation[i]];
		(B->size_inner)[i] = size_A[permutation[i]];
	}
	B->base->order = order;
	for(dim_t i = 0; i < order; i++){
		(B->base->size)[i] = size_A[permutation[i]];
		(B->base->stride)[i] = stride_A[permutation[i]];
		(B->base->size_inner)[i] = size_A[permutation[i]];
	}

	FLA_Adjust_2D_info(B);

	return FLA_SUCCESS;
}

FLA_Error FLA_Permute_hier( FLA_Obj A, dim_t permutation[], FLA_Obj* B){

	if(FLA_Obj_elemtype(A) == FLA_MATRIX || FLA_Obj_elemtype(A) == FLA_TENSOR){
		//Permute the hier's then permute inner layer
		dim_t i;
		dim_t order;
		dim_t* stride_A;
		dim_t* size_A;
		dim_t* offset_A;
		dim_t nElem = 1;
		FLA_Elemtype elemtype;

		order = FLA_Obj_order(A);
		stride_A = FLA_Obj_stride(A);
		offset_A = FLA_Obj_offset(A);
		size_A = FLA_Obj_size(A);
	
		for(i = 0; i < order; i++)
			nElem *= FLA_Obj_dimsize(A, i);

	    FLA_Obj_create_tensor_without_buffer(FLA_Obj_datatype(A), order, size_A, B);
		B->base->elemtype = FLA_Obj_elemtype(A);
		void* buf = FLA_malloc(nElem * sizeof(FLA_Obj));
		memcpy(&(((FLA_Obj*)buf)[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0]), nElem * sizeof(FLA_Obj));

		FLA_Obj_attach_buffer_to_tensor(buf, order, stride_A, B);

		//Permute tensor info
		for(dim_t i = 0; i < order; i++){
			(B->offset)[i] = offset_A[permutation[i]];
			(B->size)[i] = size_A[permutation[i]];
			(B->size_inner)[i] = size_A[permutation[i]];
		}
		B->base->order = order;
		for(dim_t i = 0; i < order; i++){
			(B->base->size)[i] = size_A[permutation[i]];
			(B->base->stride)[i] = stride_A[permutation[i]];
			(B->base->size_inner)[i] = size_A[permutation[i]];
		}

		FLA_Adjust_2D_info(B);

		//Permute inner FLA_Obj info
		FLA_Obj* buf_A = (FLA_Obj*)FLA_Obj_base_buffer(A);
		FLA_Obj* buf_B = (FLA_Obj*)FLA_Obj_base_buffer(*B);
		for(i = 0; i < nElem; i++)
			FLA_Permute_single(buf_A[i], permutation, &(buf_B[i]));		
	}
	else{
		FLA_Permute_single(A, permutation, B);
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Permute( dim_t permutation[], dim_t order, dim_t size[], FLA_Obj A[], dim_t* B_order, dim_t (*B_size)[], FLA_Obj (*B)[] )
{
  dim_t nElem = 1;
  dim_t stride_A[FLA_MAX_ORDER];
  dim_t stride_B[FLA_MAX_ORDER];
  dim_t idx_A = 0;
  dim_t idx_B = 0;
  dim_t i;
  for(i = 0; i < order; i++)
    nElem *= size[i];

  stride_A[0] = 1;
  stride_B[permutation[0]] = stride_A[0];
  for(i = 1; i < order; i++){
	stride_A[i] = stride_A[i-1]*size[i-1];
	stride_B[permutation[i]] = stride_A[i];
  }

  dim_t curIndex[order];
  memset(curIndex, 0, order * sizeof(dim_t));
  dim_t count = 0;
  dim_t ptr = 0;
  while(count < nElem){
    (*B)[idx_B].order = A[idx_A].order;
	(*B)[idx_B].base = A[idx_A].base;
//	((*B)[idx_B].base)->buffer = (A[idx_A].base)->buffer;
    for(i = 0; i < order; i++){
	  (*B)[idx_B].size[permutation[i]] = A[idx_A].size[i];
	  (*B)[idx_B].size_inner[permutation[i]] = A[idx_A].size_inner[i];
	  (((*B)[idx_B].base)->stride)[permutation[i]] = ((A[idx_A].base)->stride)[i];
	}

    count++;
    curIndex[ptr]++;
	idx_A += stride_A[ptr];
	idx_B += stride_B[ptr];

    if(curIndex[ptr] == size[ptr]){
		ptr++;
		curIndex[ptr]++;
		idx_A += stride_A[ptr];
		idx_B += stride_B[ptr];
		while(ptr < order && curIndex[ptr] == size[ptr]){
			ptr++;
			curIndex[ptr]++;
			idx_A += stride_A[ptr];
			idx_B += stride_B[ptr];
		}
		ptr--;
		memset(curIndex, 0, (ptr+1) * sizeof(dim_t));
		for(; ptr > 0; ptr--){
			idx_A -= stride_A[ptr] * size[ptr];
			idx_B -= stride_B[ptr] * size[ptr];
		}
		idx_A -= stride_A[ptr] * size[ptr];
		idx_B -= stride_B[ptr] * size[ptr];
	}
  }

  *B_order = order;
  return FLA_SUCCESS;  
}

