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

FLA_Error FLA_Random_dense_symm_tensor(dim_t nSymmGroups, dim_t symmGroupLens[nSymmGroups], dim_t** symmetries, FLA_Obj *obj){
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

		//Clear tmpA since we no longer need it
		//Will set it to tmpC right after this clear
		FLA_Obj_free_buffer(&tmpA);
		FLA_Obj_free_without_buffer(&tmpA);

		tmpA = tmpC;

		//Definitely no longer need tmpB
		FLA_Obj_free_buffer(&tmpB);
		FLA_Obj_free_without_buffer(&tmpB);

		FLA_free(tmp);
	}
	(obj->base)->buffer = (tmpC.base)->buffer;

	return FLA_SUCCESS;
}

void create_symm_groups(dim_t order, dim_t index[order], dim_t* nGroups, dim_t** groupLens, dim_t*** groups){
	dim_t i;

	(*nGroups) = 0;
	(*groupLens) = (dim_t*)FLA_malloc(order * sizeof(dim_t));
	(*groups) = (dim_t**)FLA_malloc(order * sizeof(dim_t*));

	dim_t match = index[0];
	dim_t curGroupLen = 0;
	dim_t* curGroup = (dim_t*)FLA_malloc(order * sizeof(dim_t));
	curGroup[curGroupLen] = 0;
	curGroupLen++;

	for(i = 1; i < order; i++){
		if(index[i] == match){
			curGroup[curGroupLen] = i;
			curGroupLen++;
		}
		else{
			match = index[i];
			(*groupLens)[*nGroups] = curGroupLen;
			(*groups)[*nGroups] = curGroup;
			//Update variables for next group
			(*nGroups)++;
			curGroupLen = 0;
			curGroup = (dim_t*)FLA_malloc(order * sizeof(dim_t));
			curGroup[curGroupLen] = i;
			curGroupLen++;
		}
	}
	//Do final update
	(*groupLens)[*nGroups] = curGroupLen;
	(*groups)[*nGroups] = curGroup;
	(*nGroups)++;
}

FLA_Error FLA_Obj_create_Random_symm_tensor_data(dim_t b, FLA_Obj obj){
	dim_t i;
	dim_t order = FLA_Obj_order(obj);
	dim_t* size = FLA_Obj_size(obj);
	dim_t blkSize[order];
	dim_t blkStride[order];
	dim_t curIndex[order];
	dim_t endIndex[order];
	blkSize[0] = b;
	blkStride[0] = 1;
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	endIndex[0] = size[0];
	for(i = 1; i < order; i++){
		blkSize[i] = b;
		endIndex[i] = size[i];
		blkStride[i] = blkStride[i-1]*blkSize[i-1];
	}

	dim_t nUniques = binomial(order + size[0] - 1, order);
	void** uniqueBuffers = (void**)FLA_malloc(nUniques * sizeof(void*));
	
	FLA_Obj tmpBlk;

	dim_t update_ptr = order - 1;
	dim_t count = 0;
	while(TRUE){
		//Determine symm groups
		dim_t nSymmGroups;
		dim_t* symmGroupLens;
		dim_t** symmGroups;
		create_symm_groups(order, curIndex, &nSymmGroups, &symmGroupLens, &symmGroups);
		
		//Create blk
		FLA_Obj_create_tensor(FLA_DOUBLE, order, blkSize, blkStride, &tmpBlk);
		//Freeing block because it gets created in next call
		//WARNING: NEED TO FIX THIS TO NOT BE A HACK
		FLA_Obj_free_buffer(&tmpBlk);

		FLA_Random_dense_symm_tensor(nSymmGroups, symmGroupLens, symmGroups, &tmpBlk);
		//Fill data
		uniqueBuffers[count] = tmpBlk.base->buffer;
//		FLA_Obj blk = ((FLA_Obj*) ((obj.base)->buffer))[count];
//		(blk.base)->buffer = (tmpBlk.base)->buffer;
		/**///End unique branch
		
		FLA_free(symmGroupLens);
		for(i = 0; i < nSymmGroups; i++)
			FLA_free(symmGroups[i]);
		FLA_free(symmGroups);
		
		//Update
		curIndex[update_ptr]++;
		count++;
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
	dim_t stride_obj[order];
	stride_obj[0] = 1;
	for(i = 1; i < order; i++)
		stride_obj[i] = stride_obj[i-1] * FLA_Obj_dimsize(obj, i);
	
	FLA_Obj_attach_buffer_to_symm_tensor(uniqueBuffers, order, stride_obj, &obj);

	//Clear uniqueBuffers since no longer necessary
	FLA_free(uniqueBuffers);
	FLA_free(size);
	return FLA_SUCCESS;
}
