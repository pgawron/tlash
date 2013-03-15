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

FLA_Bool next_permutation(dim_t nElem, dim_t perm[]){
	dim_t k = -1;
	dim_t i;

	for(i = 0; i < nElem - 1; i++)
		if(perm[i] < perm[i+1])
			k = i;

	if(k == -1)
		return FALSE;
	
	dim_t l = k+1;
	for(i = k+1; i < nElem; i++)
		if(perm[k] < perm[i])
			l = i;

	dim_t tmp = perm[k];
	perm[k] = perm[l];
	perm[l] = tmp;

	qsort(&(perm[k+1]), (nElem - (k+1)), sizeof(dim_t), compare_dim_t);	
	return TRUE;
}

FLA_Error FLA_Random_scalar_psym_tensor(FLA_Obj obj){
	dim_t order = FLA_Obj_order(obj);
	dim_t* size = FLA_Obj_size(obj);
	dim_t* stride = FLA_Obj_stride(obj);
	
	FLA_Random_tensor(obj);

	FLA_Obj tmpBlk;
	FLA_Obj_create_tensor(FLA_DOUBLE, FLA_Obj_order(obj), size, stride, &tmpBlk);
	FLA_Set_zero_tensor(tmpBlk);
	
	FLA_Obj permC;
	FLA_Obj_create_tensor(FLA_DOUBLE, FLA_Obj_order(obj), size, stride, &permC);

	dim_t i,j;
	dim_t modeOffset = 0;
	for(i = 0; i < obj.nSymGroups; i++){
		dim_t lenGroup = obj.symGroupLens[i];
		if(lenGroup == 1){
			modeOffset++;
			continue;
		}

		dim_t perm[lenGroup];
		for(j = 0; j < lenGroup; j++)
			perm[j] = j;
		
		//Wrap the sum into a function that blocks, blah blah...
		double* bufpermC = (double*)FLA_Obj_base_buffer(permC);
		double* buftmpBlk = (double*)FLA_Obj_base_buffer(tmpBlk);
		for(j = 0; j < lenGroup; j++)
			permC.permutation[modeOffset + j] = obj.symModes[modeOffset + perm[j]];
		FLA_Permute(obj, &(permC.permutation[0]), permC);

		//Wrap the sum into a function that blocks, blah blah...
		for(j = 0; j < FLA_Obj_num_elem_alloc(permC); j++)
			buftmpBlk[j] += bufpermC[j];


		while(next_permutation(lenGroup, perm) == TRUE){
			for(j = 0; j < lenGroup; j++)
				permC.permutation[modeOffset + j] = obj.symModes[modeOffset + perm[j]];
			
			FLA_Permute(obj, &(permC.permutation[0]), permC);
			//Wrap the sum into a function that blocks, blah blah...
			for(j = 0; j < FLA_Obj_num_elem_alloc(permC); j++)
				buftmpBlk[j] += bufpermC[j];
			
		}
		
		modeOffset += lenGroup;
		for(j = 0; j < order; j++)
			permC.permutation[j] = j;
	}
	memcpy(&(((double*)FLA_Obj_base_buffer(obj))[0]), &(((double*)FLA_Obj_base_buffer(tmpBlk))[0]), FLA_Obj_num_elem_alloc(obj) * sizeof(double));
	FLA_Obj_free_buffer(&tmpBlk);
	FLA_Obj_free_without_buffer(&tmpBlk);

	FLA_Obj_free_buffer(&permC);
	FLA_Obj_free_without_buffer(&permC);
	
	FLA_free(size);
	FLA_free(stride);
	return FLA_SUCCESS;
}

/*
FLA_Error FLA_Random_dense_sym_tensor(dim_t b, FLA_Obj* obj)
	dim_t order = FLA_Obj_order(obj);
	dim_t* size = FLA_Obj_size(obj);
	dim_t* stride = FLA_Obj_stride(obj);

	dim_t i, thisSymmGroup;
	dim_t order = FLA_Obj_order(*obj);
	dim_t modeSize = FLA_Obj_dimsize(*obj, 0);
	dim_t* tmp;

	FLA_Obj tmpA, tmpB, tmpC;
	dim_t size_A[order];
	for(i = 0; i < order; i++)
		size_A[i] = 1;
	FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size_A, &tmpA);
	double* data_A = (double*)FLA_malloc(1 * sizeof(double));
	data_A[0] = 1;
	dim_t stride_A[order];
	for(i = 0; i < order; i++)
		stride_A[i] = 1;
	FLA_Obj_attach_buffer_to_tensor(data_A, order, stride_A, &tmpA);

	for(thisSymmGroup = 0; thisSymmGroup < nSymmGroups; thisSymmGroup++)
	{
		tmp = FLA_Obj_size( tmpA );
		memcpy(&(size_A[0]), &(tmp[0]), order * sizeof(dim_t));
		dim_t size_B[] = {modeSize, 1};
		dim_t stride_B[] = {1, modeSize};
		FLA_Obj_create_tensor(FLA_DOUBLE, 2, size_B, stride_B, &tmpB);
		FLA_Random_matrix(tmpB);

		dim_t size_C[order];
		memcpy(&(size_C[0]), &(size_A[0]), order * sizeof(dim_t));
		dim_t newdimsize = FLA_Obj_dimsize(tmpB, 0);
		for(i = 0; i < symmGroupLens[thisSymmGroup]; i++)
			size_C[symmetries[thisSymmGroup][i]] = newdimsize;
		dim_t stride_C[order];
		stride_C[0] = 1;
		for(i = 1; i < order; i++)
			stride_C[i] = stride_C[i-1] * size_C[i-1];

		FLA_Obj bMults[symmGroupLens[thisSymmGroup]];
		for(i = 0; i < symmGroupLens[thisSymmGroup]; i++)
			bMults[i] = tmpB;
		FLA_Obj_create_tensor(FLA_DOUBLE, order, size_C, stride_C, &tmpC);
		FLA_Ttm(FLA_ONE, tmpA, symmGroupLens[thisSymmGroup], symmetries[thisSymmGroup], FLA_ONE, bMults, tmpC);

        //Make tmpA = tmpC and clear tmpC for next iteration
        memcpy(&(tmpA.size[0]), &(tmpC.size[0]), order * sizeof(dim_t));
        memcpy(&((tmpA.base)->size[0]), &(tmpC.size[0]), order * sizeof(dim_t));
        memcpy(&(((tmpA.base)->stride)[0]), &(((tmpC.base)->stride)[0]), order * sizeof(dim_t));
        //Clear tmpA's buffer since we will point it elsewhere
        FLA_free((tmpA.base)->buffer);

        tmpA.base->buffer = FLA_malloc(tmpC.base->n_elem_alloc * sizeof(double));
        memcpy(&(((double*)tmpA.base->buffer)[0]), &(((double*)tmpC.base->buffer)[0]), tmpC.base->n_elem_alloc * sizeof(double));
        tmpA.base->n_elem_alloc = tmpC.base->n_elem_alloc;
        FLA_Adjust_2D_info(&tmpA);


		//Definitely no longer need tmpB
		FLA_Obj_free_buffer(&tmpB);
		FLA_Obj_free_without_buffer(&tmpB);
		FLA_Obj_free(&tmpC);
		FLA_free(tmp);
	}
	memcpy(&(((double*)obj->base->buffer)[0]), &(((double*)tmpA.base->buffer)[0]), obj->base->n_elem_alloc * sizeof(double));

	FLA_Obj_free(&tmpA);

	return FLA_SUCCESS;
}
*/

void create_sym_groups(dim_t order, dim_t index[order], dim_t* nGroups, dim_t* groupLens, dim_t* groups){
	dim_t i;

	(*nGroups) = 0;

	dim_t match = index[0];
	dim_t curGroupLen = 0;
	curGroupLen++;

	for(i = 0; i < order; i++)
		groups[i] = i;

	for(i = 1; i < order; i++){
		if(index[i] == match){
			curGroupLen++;
		}
		else{
			match = index[i];
			groupLens[*nGroups] = curGroupLen;
			//Update variables for next group
			(*nGroups)++;
			curGroupLen = 1;
		}
	}
	//Do final update
	groupLens[*nGroups] = curGroupLen;
	(*nGroups)++;
}

FLA_Error FLA_Random_sym_tensor(FLA_Obj obj){
	dim_t i;
	dim_t order = FLA_Obj_order(obj);

	if(FLA_Obj_elemtype(obj) == FLA_SCALAR){
		create_sym_groups(order, &(obj.offset[0]), &(obj.nSymGroups), &(obj.symGroupLens[0]), &(obj.symModes[0]));
		FLA_Random_scalar_psym_tensor(obj);
		return FLA_SUCCESS;
	}

	dim_t* size = FLA_Obj_size(obj);
	dim_t* stride = FLA_Obj_stride(obj);
	dim_t blkSize[order];
	dim_t blkStride[order];
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t linIndex;

	dim_t b = FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(obj))[0],0);
	blkSize[0] = b;
	blkStride[0] = 1;
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	endIndex[0] = size[0];
	for(i = 1; i < order; i++){
		blkSize[i] = b;
		endIndex[i] = size[i];
		blkStride[i] = blkStride[i-1]*blkSize[i-1];
	}

	FLA_Obj tmpBlk;

	dim_t update_ptr = order - 1;
	while(TRUE){
		//Create blk
		FLA_Obj_create_tensor(FLA_DOUBLE, order, blkSize, blkStride, &tmpBlk);

		//Determine symm groups
		create_sym_groups(order, curIndex, &(tmpBlk.nSymGroups), &(tmpBlk.symGroupLens[0]), &(tmpBlk.symModes[0]));
		
		FLA_Random_scalar_psym_tensor(tmpBlk);
		//Fill data

		FLA_TIndex_to_LinIndex(order, stride, curIndex, &linIndex);
		FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(obj))[linIndex];
		double* curObjBuf = (double*)FLA_Obj_base_buffer(curObj);
		double* tmpBlkBuf = (double*)FLA_Obj_base_buffer(tmpBlk);
		memcpy(&(curObjBuf[0]), &(tmpBlkBuf[0]), FLA_Obj_num_elem_alloc(tmpBlk) * sizeof(double));

		FLA_Obj_free_buffer(&tmpBlk);
		FLA_Obj_free_without_buffer(&tmpBlk);		
		
		//Update
		curIndex[update_ptr]++;
		//Hit the end of this index
		if(curIndex[update_ptr] == endIndex[update_ptr]){
			//Keep updating previous index loops until we get a value that hasn't hit the end
			while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
				update_ptr--;
				if (update_ptr < order) {
					curIndex[update_ptr]++;
				}
			}
			//We are done if we run off the edge
			if(update_ptr >= order)
				break;
			//Otherwise, reset all upper indices to this value
			//and reset pointer
			for(i = update_ptr+1; i < order; i++)
				curIndex[i] = curIndex[update_ptr];
			update_ptr = order - 1;
		}
	}

	FLA_free(size);
	FLA_free(stride);
	return FLA_SUCCESS;
}

FLA_Error FLA_Random_psym_tensor(FLA_Obj obj){
	dim_t i,j;
	dim_t order = FLA_Obj_order(obj);

	if(FLA_Obj_elemtype(obj) == FLA_SCALAR){
		FLA_Random_scalar_psym_tensor(obj);
		return FLA_SUCCESS;
	}

	dim_t* size = FLA_Obj_size(obj);
	dim_t* stride = FLA_Obj_stride(obj);
	dim_t blkSize[order];
	dim_t blkStride[order];
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t linIndex;
	dim_t nSymGroups = obj.nSymGroups;
	dim_t symGroupLens[order];
	dim_t symModes[order];
	
	memcpy(&(symGroupLens[0]), &(obj.symGroupLens[0]), nSymGroups * sizeof(dim_t));
	memcpy(&(symModes[0]), &(obj.symModes[0]), order * sizeof(dim_t));


	dim_t b = FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(obj))[0],0);
	blkSize[0] = b;
	blkStride[0] = 1;
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	endIndex[0] = size[0];
	for(i = 1; i < order; i++){
		blkSize[i] = b;
		endIndex[i] = size[i];
		blkStride[i] = blkStride[i-1]*blkSize[i-1];
	}

	FLA_Obj tmpBlk;

	dim_t update_ptr = order - 1;
	while(TRUE){
	
		//Check if index is unique (otherwise no need to set up the random data
		dim_t isUnique = TRUE;
		dim_t count = 0;
		for(i = 0; i < nSymGroups; i++)
			if(symGroupLens[i] > 1){
				for(j = 1; j < symGroupLens[i]; j++){
					if(curIndex[symModes[count - 1]] > curIndex[symModes[count]]){
						isUnique = FALSE;
						break;
					}
					count++;
				}
				if(isUnique == FALSE)
					break;
			}else{
				count++;
			}

		if(isUnique){		
			//Create blk
			FLA_Obj_create_tensor(FLA_DOUBLE, order, blkSize, blkStride, &tmpBlk);

			//Determine symm groups
			create_sym_groups(order, curIndex, &(tmpBlk.nSymGroups), &(tmpBlk.symGroupLens[0]), &(tmpBlk.symModes[0]));
		
			FLA_Random_scalar_psym_tensor(tmpBlk);
			//Fill data

			FLA_TIndex_to_LinIndex(order, stride, curIndex, &linIndex);
			FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(obj))[linIndex];
			double* curObjBuf = (double*)FLA_Obj_base_buffer(curObj);
			double* tmpBlkBuf = (double*)FLA_Obj_base_buffer(tmpBlk);
			memcpy(&(curObjBuf[0]), &(tmpBlkBuf[0]), FLA_Obj_num_elem_alloc(tmpBlk) * sizeof(double));

			FLA_Obj_free_buffer(&tmpBlk);
			FLA_Obj_free_without_buffer(&tmpBlk);		
		}
		//Update
		curIndex[update_ptr]++;
		//Hit the end of this index
		if(curIndex[update_ptr] == endIndex[update_ptr]){
			//Keep updating previous index loops until we get a value that hasn't hit the end
			while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
				update_ptr--;
				if (update_ptr < order) {
					curIndex[update_ptr]++;
				}
			}
			//We are done if we run off the edge
			if(update_ptr >= order)
				break;
			//Otherwise, reset all upper indices to this value
			//and reset pointer
			for(i = update_ptr+1; i < order; i++)
				curIndex[i] = 0;
			update_ptr = order - 1;
		}
	}

	FLA_free(size);
	FLA_free(stride);
	return FLA_SUCCESS;
}