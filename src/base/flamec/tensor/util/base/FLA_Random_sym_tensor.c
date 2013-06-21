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

//Used to generate all permutations
//From Wikipedia article
FLA_Bool next_permutation(dim_t nElem, dim_t* perm){
	dim_t k = -1;
	dim_t i,l;
	dim_t tmp;

	for(i = 0; i < nElem - 1; i++)
		if(perm[i] < perm[i+1])
			k = i;

	if(k == -1)
		return FALSE;
	
	l = k+1;
	for(i = k+1; i < nElem; i++)
		if(perm[k] < perm[i])
			l = i;

	tmp = perm[k];
	perm[k] = perm[l];
	perm[l] = tmp;

	qsort(&(perm[k+1]), (nElem - (k+1)), sizeof(dim_t), compare_dim_t);	
	return TRUE;
}

//Creates dense psym tensor.
//Eventually should use psttm routine (must implement first)
//For now, algorithm is as follows:
//Start with scalar 1.  For each sym group:
//Create random vector and multiply in all modes of symGroup.
//Repeat for all symGroups
FLA_Error FLA_Random_scalar_psym_tensor(FLA_Obj obj){
	dim_t i,j;
	dim_t order = FLA_Obj_order(obj);
	FLA_Obj tmp;
	dim_t tmpSize[FLA_MAX_ORDER];
	dim_t tmpStride[FLA_MAX_ORDER];
	TLA_sym objSym;

	void* tmpBuffer;
	void* objBuffer;

	for(i = 0; i < order; i++){
		tmpSize[i] = 1;
		tmpStride[i] = 1;
	}
	FLA_Obj_create_tensor(FLA_DOUBLE, order, tmpSize, tmpStride, &tmp);
	((double*)((tmp.base)->buffer))[0] = 1;
	
	objSym = obj.sym;

	//Loop over all symGroups and multiply by random vector in all modes of symGroup
	for(i = 0; i < objSym.nSymGroups; i++){
		//Create vector to multiply by
		FLA_Obj vec;
		dim_t symGroupOffset = TLA_sym_group_mode_offset(objSym, i);
		dim_t symGroupMode = objSym.symModes[symGroupOffset];
		dim_t out_mode_size = FLA_Obj_dimsize(obj, symGroupMode);
		dim_t vecSize[] = {out_mode_size, 1};
		dim_t vecStride[] = {1, out_mode_size};
		FLA_Obj_create_tensor(FLA_DOUBLE, 2, vecSize, vecStride, &vec);
		FLA_Random_tensor(vec);

		//Multiply in all modes of symGroup
		for(j = 0; j < objSym.symGroupLens[i]; j++){
			//First determine which mode being multiplied in now
			dim_t mode_mult = objSym.symModes[symGroupOffset + j];
			//Set up temporary
			FLA_Obj tmpRes;
			dim_t tmpResSize[FLA_MAX_ORDER];
			dim_t tmpResStride[FLA_MAX_ORDER];
			memcpy(&(tmpResSize[0]), &(tmp.size[0]), order * sizeof(dim_t));
			tmpResSize[mode_mult] = out_mode_size;
			FLA_Set_tensor_stride(order, tmpResSize, tmpResStride);
			FLA_Obj_create_tensor(FLA_DOUBLE, order, tmpResSize, tmpResStride, &tmpRes);
			FLA_Set_zero_tensor(tmpRes);

			//Perform multiply
			FLA_Ttm_scalar_permC(FLA_ONE, tmp, mode_mult, FLA_ONE, vec, tmpRes);

			//Output is input for next iteration
			FLA_Obj_free_buffer(&tmp);
			FLA_Obj_free_without_buffer(&tmp);
			tmp = tmpRes;
		}

		//Clear vector for next symGroup
		FLA_Obj_free_buffer(&vec);
		FLA_Obj_free_without_buffer(&vec);
	}

	//Final result created, copy over to obj
	tmpBuffer = FLA_Obj_base_buffer(tmp);
	objBuffer = FLA_Obj_base_buffer(obj);
	memcpy(objBuffer, tmpBuffer, FLA_Obj_num_elem_alloc(obj) * sizeof(double));

	//Free locally alloc'd data
	FLA_Obj_free_buffer(&tmp);
	FLA_Obj_free_without_buffer(&tmp);
	return FLA_SUCCESS;
}


FLA_Error FLA_Random_psym_tensor(FLA_Obj obj){
	dim_t i;

	//Obj data
	dim_t order = FLA_Obj_order(obj);
	dim_t size[FLA_MAX_ORDER];
	dim_t stride[FLA_MAX_ORDER];
	dim_t blk_size[FLA_MAX_ORDER];
	dim_t blk_stride[FLA_MAX_ORDER];
	dim_t n_elem_alloc;
	FLA_Obj* obj_buffer = (FLA_Obj*)FLA_Obj_base_buffer(obj);

	//Loop data
	dim_t curIndex[FLA_MAX_ORDER];
	dim_t linIndex;
	dim_t isUnique;

	//If scalar, call appropriate subroutine
	if(FLA_Obj_elemtype(obj) == FLA_SCALAR){
		return FLA_Random_scalar_psym_tensor(obj);
	}

	//Initialize Obj data
	memcpy(&(size[0]), &(obj.size[0]), order * sizeof(dim_t));
	memcpy(&(stride[0]), &(((obj.base)->stride)[0]), order * sizeof(dim_t));
	memcpy(&(blk_size[0]), &(obj_buffer[0].size[0]), order * sizeof(dim_t));
	FLA_Set_tensor_stride(order, blk_size, blk_stride);
	n_elem_alloc = FLA_Obj_num_elem_alloc(obj);


	//Loop over blocks, if they're unique, fill the blocks
	for(i = 0; i < n_elem_alloc; i++){
		linIndex = i;
		FLA_LinIndex_to_TIndex(order, stride, linIndex, curIndex);
		isUnique = obj_buffer[linIndex].isStored;

		if(isUnique){
			//Determine sym groups
			FLA_Obj* curObj = &(obj_buffer[linIndex]);
			memcpy(&((curObj->offset)[0]), &(curIndex[0]), order * sizeof(dim_t));
			TLA_update_sym_based_offset(obj.sym, curObj);

			//Fill block with data
			FLA_Random_scalar_psym_tensor(*curObj);
		}
	}

	return FLA_SUCCESS;
}
