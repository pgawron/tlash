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

FLA_Error FLA_Obj_create_blocked_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj){

  FLA_Obj_create_blocked_tensor_ext( datatype, FLA_SCALAR, order, size, size, stride, blkSize, obj);

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], FLA_Obj *obj)
{

  FLA_Obj_create_tensor_ext( datatype, FLA_SCALAR, order, size, size, stride, obj );

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_blocked_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj )
{
  dim_t i,j;
  dim_t sizeObj[order];
  dim_t strideObj[order];
  dim_t nBlks = 1;
  strideObj[0] = 1;
  sizeObj[0] = size[0] / blkSize[0];

  for(i = 1; i < order; i++){
	sizeObj[i] = size[i] / blkSize[i];
	strideObj[i] = strideObj[i-1] * sizeObj[i-1];
  }
  //First set up the obj to store the FLA_Objs
  FLA_Obj_create_tensor_ext( datatype, FLA_TENSOR, order, sizeObj, sizeObj, strideObj, obj);

  FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

  for(i = 0; i < order; i++)
	nBlks *= sizeObj[i];

  //Create the blocks of obj
  for(i = 0; i < nBlks; i++){
    dim_t strideBlk[order];
    strideBlk[0] = 1;
    for(j = 1; j < order; j++)
		strideBlk[j] = strideBlk[j-1] * blkSize[j-1];
	FLA_Obj_create_tensor( datatype, order, blkSize, strideBlk, &(buf[i]));
  }

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], FLA_Obj *obj )
{
  dim_t i;
  dim_t nSecondDim = 1;
  for(i = 1; i < order; i++)
	nSecondDim *= size[i];

  //First set the 2-D object info
  FLA_Obj_create_ext( datatype, elemtype, size[0], nSecondDim, size[0], nSecondDim, stride[0], stride[0]*size[0], obj );

  //Update the tensor info
  obj->order = order;
  memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

  obj->base->order = order;
  memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof( dim_t ) );
  memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );

  obj->base->n_elem_alloc = size[0] * nSecondDim;

  for(i = 0; i < order; i++)
	obj->permutation[i] = i;
  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_tensor_without_buffer( FLA_Datatype datatype, dim_t order, dim_t size[order], FLA_Obj *obj ){
	dim_t i;
	dim_t nSecondDim = 1;
	for(i = 1; i < order; i++)
		nSecondDim *= size[i];

	//First set 2-D object info
	FLA_Obj_create_without_buffer(datatype, size[0], nSecondDim, obj);

	//Update tensor info
	obj->order = order;
    memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

    obj->base->order = order;
    memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->base->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );
    memset(&((obj->base->stride)[0]), 0, order * sizeof( dim_t ) );

	obj->base->n_elem_alloc = size[0] * nSecondDim;
	for(i = 0; i < order; i++)
		(obj->permutation)[i] = i;

	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, dim_t stride[order], FLA_Obj *obj ){
	//Omitting some things attach_buffer does because not sure how to extend yet
	obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

	FLA_Adjust_2D_info(obj);
	
	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_blocked_symm_free_buffer( FLA_Obj *obj)
{
	dim_t order = FLA_Obj_order(*obj);
	dim_t* endIndex = FLA_Obj_size(*obj);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	dim_t* stride = FLA_Obj_stride(*obj);
	FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

	dim_t update_ptr = order - 1;
	while(TRUE){
		dim_t linIndex;
		FLA_TIndex_to_LinIndex(order, curIndex, stride, &linIndex);

		FLA_Obj_free_buffer(&(buf[linIndex]));
		FLA_Obj_free_without_buffer(&(buf[linIndex]));
		
		//Update
		curIndex[update_ptr]++;
		while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
			update_ptr--;
			if(update_ptr < order)
				curIndex[update_ptr]++;
		}
		if(update_ptr >= order)
			break;
		for(dim_t i = update_ptr+1; i < order; i++)
			curIndex[i] = curIndex[update_ptr];
		update_ptr = order - 1;
	}
	FLA_Obj_free_buffer(obj);
	FLA_free(endIndex);
	FLA_free(stride);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_symm_tensor_without_buffer(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t blkSize, FLA_Obj *obj){
	dim_t i;
	dim_t nTBlks;
	dim_t nBlksPerMode;
	FLA_Obj* t_blks;
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t updateIndex;
	dim_t objLinIndex;
	dim_t sizeBlock[order];

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Obj_create_symm_tensor_without_buffer_check( datatype, order, size, blkSize, obj );

	//Figure out how many FLA_Objs we need for intermediate level
	nTBlks = 1;
	nBlksPerMode = size[0] / blkSize;
	for(i = 0; i < order; i++)
		nTBlks *= nBlksPerMode;

	for(i = 0; i < order; i++)
		sizeBlock[i] = blkSize;

	//Create the FLA_Objs
	t_blks = (FLA_Obj*)FLA_malloc(nTBlks * sizeof(FLA_Obj));

	//curIndex is our counter for each block's logical index
	//We loop over this until we hit endIndex  and create the apropriate FLA_Objs
	//updateIndex tells us which index in curIndex we need to update
	//objLinIndex tells us which linear object we are dealing with
	memset(curIndex, 0, order * sizeof(dim_t));
	for(i = 0; i < order; i++)
		endIndex[i] = nBlksPerMode;
	updateIndex = order - 1;
	objLinIndex = 0;

	while(TRUE){
		//Set up the FLA_Obj at this index
		FLA_Obj *curObj = &(t_blks[objLinIndex]);

		FLA_Obj_create_tensor_without_buffer( datatype, order, sizeBlock, curObj);
		//Set the offset array (we will use as an index identifier)
		memcpy(&((curObj->offset)[0]), &(curIndex[0]), order * sizeof(dim_t));

		//Loop update
		//Update current index
		curIndex[updateIndex]++;
		objLinIndex++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
			updateIndex--;
			if(updateIndex < order)
				curIndex[updateIndex]++;
		}
		//If we run off the edge, we know we are at the end, so break out
		if(updateIndex >= order)
			break;
		//Otherwise, update current index, and reset all others
		for(i = updateIndex + 1; i < order; i++)
			curIndex[i] = 0;
		updateIndex = order - 1;
	}

	//Buffer of tensor blocks created, set the main obj to represent this hierarchy
	
	FLA_Obj_create_tensor_without_buffer( datatype, order, endIndex, obj);
	obj->base->elemtype = FLA_TENSOR;

	dim_t* size_obj = FLA_Obj_size(*obj);
	dim_t stride_obj[order];
	stride_obj[0] = 1;
	for(i = 1; i < order; i++)
		stride_obj[i] = stride_obj[i-1]*size_obj[i-1];
	
	FLA_Obj_attach_buffer_to_tensor(t_blks, order, stride_obj, obj);

	FLA_Adjust_2D_info(obj);

	FLA_free(size_obj);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_attach_buffer_to_symm_tensor( void *buffer[], dim_t order, dim_t stride[order], FLA_Obj *obj ){
	dim_t i, updateIndex;
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t* size_obj;
	dim_t countBuffer;
	dim_t objLinIndex;
	dim_t sortedIndex[order];
	dim_t permutation[order];
	dim_t* stride_obj;
	FLA_Obj *buffer_obj;

	size_obj = FLA_Obj_size(*obj);
	buffer_obj = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
	stride_obj = FLA_Obj_stride(*obj);

	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	memcpy(&(endIndex[0]), &(size_obj[0]), order * sizeof(dim_t));
	countBuffer = 0;
	objLinIndex = 0;

	FLA_Paired_Sort index_pairs[order];
	//Loop over indices
	//If we hit a unique, set its buffer to the next in the list
	//Otherwise, find the unique index and set the buffer to that buffer (adjusting strides)
	//By looping correctly we will hit the unique before any dupes (I think)
	updateIndex = order - 1;
	while(TRUE){
		objLinIndex = 0;
		for(i = 0; i < order; i++)
			objLinIndex += curIndex[i] * stride_obj[i];

		for(i = 0; i < order; i++){
			index_pairs[i].index = i;
			index_pairs[i].val = curIndex[i];
		}

		qsort(index_pairs, order, sizeof(FLA_Paired_Sort), compare_pairwise_sort);

		for(i = 0; i < order; i++){
			permutation[i] = index_pairs[i].index;
			sortedIndex[i] = index_pairs[i].val;
		}

		//Check if this is unique or not
		dim_t uniqueIndex = TRUE;
		for(i = 0; i < order; i++){
			if(sortedIndex[i] != curIndex[i]){
				uniqueIndex = FALSE;
				break;
			}
		}
		if(uniqueIndex){
			(buffer_obj[objLinIndex].base)->buffer = buffer[countBuffer];
			//Update stride - TODO: MOVE THIS ELSEWHERE
			((buffer_obj[objLinIndex].base)->stride)[0] = 1;
			for(i = 1; i < order; i++)
				((buffer_obj[objLinIndex].base)->stride)[i] = ((buffer_obj[objLinIndex].base)->stride)[i-1] * ((buffer_obj[objLinIndex].base)->size)[i-1];
			countBuffer++;
		}else{
			dim_t ipermutation[order];
			dim_t uniqueLinIndex;
			FLA_TIndex_to_LinIndex(order, stride_obj, sortedIndex, &(uniqueLinIndex));

			//point this non-unique FLA_Obj to the correct base
			//WARNING: HACK
			FLA_free(buffer_obj[objLinIndex].base);
			(buffer_obj[objLinIndex]).base = (buffer_obj[uniqueLinIndex]).base;
			//Set the right permutation
			memcpy(&(ipermutation[0]), &(permutation[0]), order * sizeof(dim_t));
			for(i = 0; i < order; i++)
				ipermutation[permutation[i]] = i;

//			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(permutation[0]), order * sizeof(dim_t));
			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(ipermutation[0]), order * sizeof(dim_t));
		}
		FLA_Adjust_2D_info(&(buffer_obj[objLinIndex]));
		//Loop update
		//Update current index
		curIndex[updateIndex]++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
			updateIndex--;
			if(updateIndex < order)
				curIndex[updateIndex]++;
		}
		//If we run off the edge, we know we are at the end, so break out
		if(updateIndex >= order)
			break;
		//Otherwise, update current index, and reset all others
		for(i = updateIndex+1; i < order; i++)
			curIndex[i] = 0;
		updateIndex = order - 1;
	}

	//Omitting some things attach_buffer does because not sure how to extend yet
	//obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

	FLA_free(size_obj);
	FLA_free(stride_obj);
	return FLA_SUCCESS;
}
