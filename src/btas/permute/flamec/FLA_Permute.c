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
FLA_Error FLA_Permute_single_inplace( FLA_Obj* A, dim_t permutation[]){
	
    FLA_Datatype datatype = FLA_Obj_datatype( *A );	
	dim_t order;
	dim_t* stride_A_old;
	dim_t* size_A_old;
	
	order = FLA_Obj_order(*A);
	stride_A_old = FLA_Obj_stride(*A);
	size_A_old = FLA_Obj_size(*A);
	
	//Permute tensor info
	for(dim_t i = 0; i < order; i++){
		(A->size)[i] = size_A_old[permutation[i]];
		(A->size_inner)[i] = size_A_old[permutation[i]];
	}
	for(dim_t i = 0; i < order; i++){
		(A->base->size)[i] = size_A_old[permutation[i]];
		(A->base->stride)[i] = stride_A_old[permutation[i]];
		(A->base->size_inner)[i] = size_A_old[permutation[i]];
	}
	(A->base->stride)[0] = 1;
	for(dim_t i = 1; i < order; i++)
		(A->base->stride)[i] = (A->base->stride)[i-1] * (A->base->size)[i-1];
	
	dim_t nElem = 1;
	for(dim_t i = 0; i < order; i++)
		nElem *= (A->base->size)[i];
	
	void* buffer;
	if(FLA_Obj_elemtype(*A) == FLA_SCALAR)
		buffer = FLA_malloc(nElem * sizeof(double));
	else
		buffer = FLA_malloc(nElem * sizeof(FLA_Obj));	
	
    //Explicitely permute the data as well
	void* buf_A = FLA_Obj_base_buffer(*A);
	dim_t* newStride = FLA_Obj_stride(*A);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	dim_t* endIndex = FLA_Obj_size(*A);
	dim_t updatePtr = order - 1;
    while(TRUE){
		//Calculate linear index fro and to
		dim_t linIndexFro = 0;
		dim_t linIndexTo = 0;
		
		dim_t permutedIndex[order];
		FLA_Permute_array(order, curIndex, permutation, &(permutedIndex[0]));
		
		FLA_TIndex_to_LinIndex(order, curIndex, stride_A_old, &linIndexFro);
		FLA_TIndex_to_LinIndex(order, permutedIndex, newStride, &linIndexTo);
		
		if(FLA_Obj_elemtype(*A) == FLA_SCALAR)
			((double*)buffer)[linIndexTo] = ((double*)buf_A)[linIndexFro];
		else
			((FLA_Obj*)buffer)[linIndexTo] = ((FLA_Obj*)buf_A)[linIndexFro];
		
		//Update
		curIndex[updatePtr]++;
		while(updatePtr < order && curIndex[updatePtr] == endIndex[updatePtr]){
			updatePtr--;
			if(updatePtr < order)
				curIndex[updatePtr]++;
		}
		if(updatePtr >= order)
			break;
		for(dim_t i = updatePtr+1; i < order; i++)
			curIndex[i] = 0;
		updatePtr = order - 1;
	}

	//Copy the data from temporary buffer to A's buffer
	dim_t numElem = FLA_Obj_num_elem_alloc(*A);
	for (dim_t i = 0; i < numElem; i++) {
		if(FLA_Obj_elemtype(*A) == FLA_SCALAR)
			((double*)FLA_Obj_base_buffer(*A))[i] = ((double*)buffer)[i];
		else
			((FLA_Obj*)FLA_Obj_base_buffer(*A))[i] = ((FLA_Obj*)buffer)[i];
	}
	
//	FLA_Obj_attach_buffer_to_tensor(buffer, order, A->base->stride, A);
	
	FLA_Adjust_2D_info(A);
	
	return FLA_SUCCESS;	
}


/***
 *
 *	TODO: Need consistency between what a View is.  Right now, I need to permute
 *  data first based on the permutation array in the View, then based on the
 *  permutation given in the parameter.
 *  Does this mean that the FLA_View size value should reflect the already
 *  permuted data? (same with other associated fields.
 *
 ***/
FLA_Error FLA_Permute_single( FLA_Obj A, dim_t permutation[], FLA_Obj* B){
	
	dim_t i;
	FLA_Datatype datatype = FLA_Obj_datatype(A);
	dim_t order;
	dim_t* stride_A;
	dim_t* size_A;
	dim_t* offset_A;
	dim_t* size_B;
	dim_t nElem_A;

	order = FLA_Obj_order(A);
	stride_A = FLA_Obj_stride(A);
	offset_A = FLA_Obj_offset(A);
	size_A = FLA_Obj_size(A);
	nElem_A = FLA_Obj_num_elem_alloc(A);
	size_B = FLA_Obj_size(*B);
	dim_t stride_B[order];
	dim_t offset_B[order];
	
	dim_t viewPermuted = FALSE;
	for(i = 0; i < order; i++){
		if(A.permutation[i] != i){
			viewPermuted = TRUE;
			break;
		}
	}

	if(viewPermuted){
		FLA_Obj unpermA, tmp;
		dim_t sizeTmp[order];
		dim_t strideTmp[order];
		FLA_Permute_array(order, size_A, A.permutation, &(sizeTmp[0]));
		FLA_Set_tensor_stride(order, sizeTmp, &(strideTmp[0]));

		FLA_Obj_create_tensor_without_buffer(datatype, order, size_A, &unpermA);
		FLA_Obj_create_tensor_without_buffer(datatype, order, sizeTmp, &tmp);

		dim_t numElemTmp = 1;
		for(i = 0; i < order; i++)
			numElemTmp *= sizeTmp[i];

		void* tmpBuf = FLA_malloc(numElemTmp * sizeof(double));

		unpermA.base = A.base;
		memcpy(&(unpermA.offset[0]), &(A.offset[0]), order * sizeof(dim_t));
		FLA_Adjust_2D_info(&unpermA);

		FLA_Obj_attach_buffer_to_tensor(tmpBuf, order, strideTmp, &tmp);
		FLA_Permute_single(unpermA, A.permutation, &tmp);
		FLA_Permute_single(tmp, permutation, B);
		return FLA_SUCCESS;
	}
	
	FLA_Set_tensor_stride(order, size_B, &(stride_B[0]));
	
	void* buffer;
	if(FLA_Obj_elemtype(A) == FLA_SCALAR)
		buffer = FLA_malloc(nElem_A * sizeof(double));
	else
		buffer = FLA_malloc(nElem_A * sizeof(FLA_Obj));	

    //Explicitely permute the data
	void* buf_A = FLA_Obj_base_buffer(A);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	
	dim_t updatePtr = order - 1;
    while(TRUE){
		//Calculate linear index fro and to
		dim_t linIndexFro = 0;
		dim_t linIndexTo = 0;

		FLA_TIndex_to_LinIndex(order, stride_A, curIndex, &(linIndexFro));

		dim_t permutedIndex[order];
		FLA_Permute_array(order, curIndex, permutation, &(permutedIndex[0]));
		//Check if this math is right.
		FLA_TIndex_to_LinIndex(order, stride_B, permutedIndex, &(linIndexTo));

		if(FLA_Obj_elemtype(A) == FLA_SCALAR)
			((double*)buffer)[linIndexTo] = ((double*)buf_A)[linIndexFro];
		else
			((FLA_Obj*)buffer)[linIndexTo] = ((FLA_Obj*)buf_A)[linIndexFro];

		//Update
		curIndex[updatePtr]++;
		while(updatePtr < order && curIndex[updatePtr] == size_A[updatePtr]){
			updatePtr--;
			if(updatePtr < order)
				curIndex[updatePtr]++;
		}
		if(updatePtr >= order)
			break;
		for(dim_t i = updatePtr+1; i < order; i++)
			curIndex[i] = 0;
		updatePtr = order - 1;
	}

	dim_t nElemB = 1;
	for(i = 0; i < order; i++)
		nElemB *= FLA_Obj_dimsize(*B, i);
	
	if(FLA_Obj_elemtype(A) == FLA_SCALAR){
		double* buf_B = FLA_Obj_base_buffer(*B);
		memcpy(&(buf_B[0]), &(((double*)buffer)[0]), nElemB * sizeof(double));
	}else{
		FLA_Obj* buf_B = FLA_Obj_base_buffer(*B);
		memcpy(&(buf_B[0]), &(((double*)buffer)[0]), nElemB * sizeof(FLA_Obj*));
	}
	memcpy(&((B->permutation)[0]), &(A.permutation[0]), order * sizeof(dim_t));
	
	FLA_Adjust_2D_info(B);

	return FLA_SUCCESS;
}

FLA_Error FLA_Permute_hier( FLA_Obj A, dim_t permutation[], FLA_Obj* B){

	printf("A pre permute_single in hier:\n");
	FLA_Obj_print_tensor(A);

//	printf("B pre permute_single in hier:\n");
//	FLA_Obj_print_tensor(*B);
	
	FLA_Permute_single(A, permutation, B);

	if(FLA_Obj_elemtype(A) == FLA_MATRIX || FLA_Obj_elemtype(A) == FLA_TENSOR){
		//Permute inner FLA_Obj info
		dim_t nElem = 1;
		for(dim_t i = 0; i < FLA_Obj_order(A); i++)
			nElem *= FLA_Obj_dimsize(A, i);
		FLA_Obj* buf_B = (FLA_Obj*)FLA_Obj_base_buffer(*B);
		for(dim_t i = 0; i < nElem; i++){
			//Only permute the unique entries
			//Othersize just adjust the permutation field
			dim_t order = FLA_Obj_order(A);
			dim_t isUnique = TRUE;
			for(dim_t j = 0; j < order; j++)
				if(A.permutation[j] != j){
					isUnique = FALSE;
					break;
				}
			if(isUnique)
				FLA_Permute_single_inplace(&(buf_B[i]), permutation);
			else{
				dim_t newPerm[order];
				FLA_Permute_array(order, &(A.permutation[0]), permutation, &(newPerm[0]));
				for(dim_t j = 0; j < order; j++)
					A.permutation[j] = newPerm[j];
			}
		}
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

