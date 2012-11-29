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

/***
 *
 *  PERMUTES A SCALAR TENSOR ONLY!!!!!
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

	dim_t ipermutation[order];
	for(i = 0; i < order; i++)
		ipermutation[permutation[i]] = i;
	
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

		//WARNING: HACK NEED TO FIX THIS
		FLA_free(unpermA.base);
		unpermA.base = A.base;
		memcpy(&(unpermA.offset[0]), &(A.offset[0]), order * sizeof(dim_t));
		FLA_Adjust_2D_info(&unpermA);

		FLA_Obj_attach_buffer_to_tensor(tmpBuf, order, strideTmp, &tmp);
		FLA_Permute_single(unpermA, A.permutation, &tmp);
		FLA_Permute_single(tmp, permutation, B);

		FLA_Obj_free_buffer(&tmp);
		FLA_Obj_free_without_buffer(&tmp);
		FLA_free(stride_A);
		FLA_free(size_A);
		FLA_free(offset_A);
		FLA_free(size_B);
		return FLA_SUCCESS;
	}
	
	FLA_Set_tensor_stride(order, size_B, &(stride_B[0]));
	
	void* buffer;
	buffer = FLA_malloc(nElem_A * sizeof(double));

    //Explicitely permute the data
	void* buf_A = FLA_Obj_base_buffer(A);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));

	dim_t updatePtr = 0;
	dim_t linIndexFro = 0;
	dim_t linIndexTo = 0;
    while(TRUE){
		//Calculate linear index fro and to
//		dim_t linIndexFro = 0;
//		dim_t linIndexTo = 0;

//		FLA_TIndex_to_LinIndex(order, stride_A, curIndex, &(linIndexFro));

//		dim_t permutedIndex[order];
//		FLA_Permute_array(order, curIndex, permutation, &(permutedIndex[0]));
		//Check if this math is right.
//		FLA_TIndex_to_LinIndex(order, stride_B, permutedIndex, &(linIndexTo));


		((double*)buffer)[linIndexTo] = ((double*)buf_A)[linIndexFro];

		//Update
		curIndex[updatePtr]++;
		linIndexFro += stride_A[updatePtr];
		linIndexTo += stride_B[ipermutation[updatePtr]];
		while(updatePtr < order && curIndex[updatePtr] == size_A[updatePtr]){
			updatePtr++;
			if(updatePtr < order){
				curIndex[updatePtr]++;
				linIndexFro += stride_A[updatePtr];
				linIndexTo += stride_B[ipermutation[updatePtr]];
			}
		}
		if(updatePtr >= order)
			break;
		for(dim_t i = updatePtr-1; i < order; i--){
			curIndex[i] = 0;
			linIndexFro -= stride_A[i]*size_A[i];
			linIndexTo -= stride_B[ipermutation[i]]*size_B[ipermutation[i]];
		}
		updatePtr = 0;
	}

	dim_t nElemB = 1;
	for(i = 0; i < order; i++)
		nElemB *= FLA_Obj_dimsize(*B, i);
	
	double* buf_B = FLA_Obj_base_buffer(*B);
	memcpy(&(buf_B[0]), &(((double*)buffer)[0]), nElemB * sizeof(double));

	memcpy(&((B->permutation)[0]), &(A.permutation[0]), order * sizeof(dim_t));

	FLA_Adjust_2D_info(B);

	FLA_free(buffer);	
	FLA_free(stride_A);
	FLA_free(size_A);
	FLA_free(offset_A);
	FLA_free(size_B);

	return FLA_SUCCESS;
}
