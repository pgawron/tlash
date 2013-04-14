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

FLA_Error FLA_Obj_create_psym_tensor(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], TLA_sym sym, FLA_Obj *obj){
  FLA_Obj_create_tensor( datatype, order, size, stride, obj);

  //Update symmetries
  obj->groupToPartition = 0;
  obj->sym = sym;

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

  //Update symmetries
  obj->groupToPartition = 0;
  (obj->sym).order = order;
  (obj->sym).nSymGroups = order;
  for(i = 0; i < (obj->sym).nSymGroups; i++)
        ((obj->sym).symGroupLens)[i] = 1;
  for(i = 0; i < order; i++)
      ((obj->sym).symModes)[i] = i;

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
    obj->isStored = FALSE;

    obj->base->order = order;
    memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->base->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );
    memset(&((obj->base->stride)[0]), 0, order * sizeof( dim_t ) );

	obj->base->n_elem_alloc = size[0] * nSecondDim;
	for(i = 0; i < order; i++)
		(obj->permutation)[i] = i;

    //Update symmetries
	obj->groupToPartition = 0;
	(obj->sym).order = order;
    (obj->sym).nSymGroups = order;
    for(i = 0; i < (obj->sym).nSymGroups; i++)
        ((obj->sym).symGroupLens)[i] = 1;
    for(i = 0; i < order; i++)
        ((obj->sym).symModes)[i] = i;

	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, dim_t stride[order], FLA_Obj *obj ){
	//Omitting some things attach_buffer does because not sure how to extend yet
	obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

	FLA_Adjust_2D_info(obj);
	
	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_blocked_sym_tensor_free_buffer( FLA_Obj *obj)
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

FLA_Error FLA_Obj_blocked_psym_tensor_free_buffer( FLA_Obj *obj)
{
	dim_t order = FLA_Obj_order(*obj);
	dim_t* endIndex = FLA_Obj_size(*obj);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	dim_t* stride = FLA_Obj_stride(*obj);
	FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

	dim_t nSymGroups = (obj->sym).nSymGroups;
	dim_t symGroupLens[order];
	dim_t symModes[order];
	memcpy(&(symGroupLens[0]), &(((obj->sym).symGroupLens)[0]), nSymGroups * sizeof(dim_t));
	memcpy(&(symModes[0]), &(((obj->sym).symModes)[0]), order * sizeof(dim_t));

	dim_t update_ptr = order - 1;

	while(TRUE){
		dim_t linIndex;
		FLA_TIndex_to_LinIndex(order, curIndex, stride, &linIndex);

		dim_t isUnique = buf[linIndex].isStored;

		if(isUnique){
		    printf("at index %d ", linIndex);
		    print_array("freeing", order, curIndex);

			FLA_Obj_free_buffer(&(buf[linIndex]));
			FLA_Obj_free_without_buffer(&(buf[linIndex]));
		}
		
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
			curIndex[i] = 0;
		update_ptr = order - 1;
	}
	FLA_Obj_free_buffer(obj);
	FLA_free(endIndex);
	FLA_free(stride);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_blocked_sym_tensor_without_buffer(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t blkSize, FLA_Obj *obj){
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
		FLA_Obj_create_blocked_sym_tensor_without_buffer_check( datatype, order, size, blkSize, obj );

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
	updateIndex = 0;
	objLinIndex = 0;

	while(TRUE){
		//Set up the FLA_Obj at this index
		FLA_Obj *curObj = &(t_blks[objLinIndex]);

		FLA_Obj_create_tensor_without_buffer( datatype, order, sizeBlock, curObj);
		//Set the offset array (we will use as an index identifier)
		memset(&((curObj->offset)[0]), 0, order * sizeof(dim_t));
		//memcpy(&((curObj->offset)[0]), &(curIndex[0]), order * sizeof(dim_t));

		//Loop update
		//Update current index
		curIndex[updateIndex]++;
		objLinIndex++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
			updateIndex++;
			if(updateIndex < order)
				curIndex[updateIndex]++;
		}
		//If we run off the edge, we know we are at the end, so break out
		if(updateIndex >= order)
			break;
		//Otherwise, update current index, and reset all others
		for(i = updateIndex - 1; i < order; i--)
			curIndex[i] = 0;
		updateIndex = 0;
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

FLA_Error FLA_Obj_create_blocked_psym_tensor_without_buffer(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t blkSize, FLA_Obj *obj){
    dim_t i;
    dim_t nTBlks;
    FLA_Obj* t_blks;
    dim_t curIndex[order];
    dim_t endIndex[order];
    dim_t updateIndex;
    dim_t objLinIndex;
    dim_t sizeBlock[order];

//    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
//        FLA_Obj_create_blocked_sym_tensor_without_buffer_check( datatype, order, size, blkSize, obj );

    //Figure out how many FLA_Objs we need for intermediate level
    nTBlks = 1;

    for(i = 0; i < order; i++)
        nTBlks *= size[i] / blkSize;

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
        endIndex[i] = size[i] / blkSize;
    updateIndex = 0;
    objLinIndex = 0;

    while(TRUE){
        //Set up the FLA_Obj at this index
        FLA_Obj *curObj = &(t_blks[objLinIndex]);

        FLA_Obj_create_tensor_without_buffer( datatype, order, sizeBlock, curObj);
        //Set the offset array (we will use as an index identifier)
        memset(&((curObj->offset)[0]), 0, order * sizeof(dim_t));
        //memcpy(&((curObj->offset)[0]), &(curIndex[0]), order * sizeof(dim_t));

        //Loop update
        //Update current index
        curIndex[updateIndex]++;
        objLinIndex++;
        //If we hit the end, loop until we find the index to update
        while(updateIndex < order && curIndex[updateIndex] == endIndex[updateIndex]){
            updateIndex++;
            if(updateIndex < order)
                curIndex[updateIndex]++;
        }
        //If we run off the edge, we know we are at the end, so break out
        if(updateIndex >= order)
            break;
        //Otherwise, update current index, and reset all others
        for(i = updateIndex - 1; i < order; i--)
            curIndex[i] = 0;
        updateIndex = 0;
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

FLA_Error FLA_Obj_create_blocked_sym_tensor(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], dim_t blkSize, FLA_Obj *obj){

	//First set up the hierarchy without buffers
	FLA_Obj_create_blocked_sym_tensor_without_buffer(datatype, order, size, blkSize, obj);

	//Set up the data buffers for sym tensor
	dim_t i;
	dim_t nBlockElems = 1;
	for(i = 0; i < order; i++)
		nBlockElems *= blkSize;

	dim_t nUniques = binomial(order + (size[0] / blkSize) - 1,order);
	void** dataBuffers = (void**)FLA_malloc(nUniques * sizeof(void*));

	for(i = 0; i < nUniques; i++)
		dataBuffers[i] = (double*)FLA_malloc(nBlockElems * sizeof(double));

	//Attach empty buffers to the sym tensor
	FLA_Obj_attach_buffer_to_blocked_sym_tensor(dataBuffers, order, stride, obj);
	FLA_free(dataBuffers);

	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_blocked_psym_tensor(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], dim_t blkSize, TLA_sym sym, FLA_Obj *obj){

	//First set up the hierarchy without buffers
	FLA_Obj_create_blocked_psym_tensor_without_buffer(datatype, order, size, blkSize, obj);

	//Add symmetry to object
	obj->groupToPartition = 0;
	obj->sym = sym;
	
	//Set up the data buffers for psym tensor
	dim_t i,j;
	dim_t nBlockElems = 1;
	for(i = 0; i < order; i++)
		nBlockElems *= blkSize;

	dim_t nUniques = 1;
	dim_t modeOffset = 0;
	for(i = 0; i < (obj->sym).nSymGroups; i++){
		nUniques *= binomial((obj->sym).symGroupLens[i] + (size[modeOffset] / blkSize) - 1, (obj->sym).symGroupLens[i]);
		modeOffset += (obj->sym).symGroupLens[i];
	}

	void** dataBuffers = (void**)FLA_malloc(nUniques * sizeof(void*));
	for(i = 0; i < nUniques; i++){
		dataBuffers[i] = (double*)FLA_malloc(nBlockElems * sizeof(double));
		for(j = 0; j < nBlockElems; j++)
		    ((double*)dataBuffers[i])[j] = 0.0;
	}
	//Attach empty buffers to the sym tensor
	FLA_Obj_attach_buffer_to_blocked_psym_tensor(dataBuffers, order, stride, obj);

	FLA_free(dataBuffers);
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_attach_buffer_to_blocked_sym_tensor( void *buffer[], dim_t order, dim_t stride[order], FLA_Obj *obj ){
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

FLA_Error FLA_Obj_attach_buffer_to_blocked_psym_tensor( void *buffer[], dim_t order, dim_t stride[order], FLA_Obj *obj ){
	dim_t i, j;
	dim_t updateIndex;
	dim_t curIndex[order];
	dim_t endIndex[order];
	dim_t* size_obj;
	dim_t countBuffer;
	dim_t objLinIndex;
	dim_t sortedIndex[order];
	dim_t permutation[order];
	dim_t* stride_obj;
	FLA_Obj *buffer_obj;

	dim_t nSymGroups = (obj->sym).nSymGroups;
	dim_t symGroupLens[nSymGroups];
	dim_t symModes[order];

	memcpy(&(symGroupLens[0]), &(((obj->sym).symGroupLens)[0]), nSymGroups * sizeof(dim_t));
	memcpy(&(symModes[0]), &(((obj->sym).symModes)[0]), order * sizeof(dim_t));
	
	size_obj = FLA_Obj_size(*obj);
	buffer_obj = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
	stride_obj = FLA_Obj_stride(*obj);

	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	memcpy(&(endIndex[0]), &(size_obj[0]), order * sizeof(dim_t));
	countBuffer = 0;
	objLinIndex = 0;

	FLA_Paired_Sort index_pairs[order];
	dim_t orderedSymModes[order];
	//Loop over indices
	//If we hit a unique, set its buffer to the next in the list
	//Otherwise, find the unique index and set the buffer to that buffer (adjusting strides)
	//By looping correctly we will hit the unique before any dupes (I think)
	updateIndex = order - 1;
	
	while(TRUE){
	
	//FIX THIS FOR PSYM, ended here!!!!!!
	    print_array("attaching index", order, curIndex);

		FLA_TIndex_to_LinIndex(order, stride_obj, curIndex, &objLinIndex);

		dim_t modeOffset = 0;

		for(i = 0; i < nSymGroups; i++){
			
			for(j = 0; j < symGroupLens[i]; j++)
				orderedSymModes[modeOffset + j] = symModes[j+modeOffset];
			qsort(&(orderedSymModes[0]), symGroupLens[i], sizeof(dim_t), compare_dim_t);
		
			for(j = 0; j < symGroupLens[i]; j++){
				index_pairs[j].index = orderedSymModes[modeOffset + j];
				index_pairs[j].val = curIndex[orderedSymModes[modeOffset + j]];
			}
			qsort(index_pairs, symGroupLens[i], sizeof(FLA_Paired_Sort), compare_pairwise_sort);
			

			for(j = 0; j < symGroupLens[i]; j++){
				permutation[orderedSymModes[modeOffset + j]] = index_pairs[j].index;
				sortedIndex[orderedSymModes[modeOffset + j]] = index_pairs[j].val;
			}

			modeOffset += symGroupLens[i];
		}
		
		//Check if this is unique or not
		dim_t uniqueIndex = TRUE;
		dim_t count = 0;
		for(i = 0; i < nSymGroups; i++){
			if(symGroupLens[i] > 1){
				for(j = 0; j < symGroupLens[i] - 1; j++){
					if(curIndex[symModes[count]] > curIndex[symModes[count+1]]){
						uniqueIndex = FALSE;
						break;
					}
					count++;
				}
				if(uniqueIndex == FALSE)
					break;
			}
			count++;
		}

		if(uniqueIndex){
		    print_array("attaching to uniqueIndex", order, curIndex);
			(buffer_obj[objLinIndex].base)->buffer = buffer[countBuffer];
			//Update stride - TODO: MOVE THIS ELSEWHERE
			((buffer_obj[objLinIndex].base)->stride)[0] = 1;
			for(i = 1; i < order; i++)
				((buffer_obj[objLinIndex].base)->stride)[i] = ((buffer_obj[objLinIndex].base)->stride)[i-1] * ((buffer_obj[objLinIndex].base)->size)[i-1];
			for(i = 0; i < order; i++)
			    ((buffer_obj[objLinIndex]).permutation)[i] = i;
			(buffer_obj[objLinIndex]).isStored = TRUE;
			countBuffer++;
		}else{
			dim_t ipermutation[order];
			dim_t uniqueLinIndex;
			FLA_TIndex_to_LinIndex(order, stride_obj, sortedIndex, &(uniqueLinIndex));

			//point this non-unique FLA_Obj to the correct base
			//WARNING: HACK
			print_array("freeing base of ", order, curIndex);
			FLA_free(buffer_obj[objLinIndex].base);
			(buffer_obj[objLinIndex]).base = (buffer_obj[uniqueLinIndex]).base;
			//Set the right permutation
			memcpy(&(ipermutation[0]), &(permutation[0]), order * sizeof(dim_t));
			modeOffset = 0;
			for(i = 0; i < nSymGroups; i++){
				for(j = 0; j < symGroupLens[i]; j++){
					ipermutation[permutation[symModes[j+modeOffset]]] = orderedSymModes[j+modeOffset];
				}
				modeOffset += symGroupLens[i];
			}
			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(ipermutation[0]), order * sizeof(dim_t));
			(buffer_obj[objLinIndex]).isStored = FALSE;
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

////////////////////////////////
//More of Util functions follow
////////////////////////////////

//Note: Only splits modes within the same symmetric group...
FLA_Error TLA_split_sym_group(TLA_sym S, dim_t nSplit_modes, dim_t split_modes[nSplit_modes], TLA_sym* S1){

    if(nSplit_modes == 0){
        *S1 = S;
        return FLA_SUCCESS;
    }

    dim_t sym_group = TLA_sym_group_of_mode(S, split_modes[0]);
    if(TLA_sym_group_size(S, sym_group) == 1){
        S1->order = S.order;
	    S1->nSymGroups = S.nSymGroups;
	    memcpy(&((S1->symGroupLens)[0]), &(S.symGroupLens[0]), S.nSymGroups * sizeof(dim_t));
	    memcpy(&((S1->symModes)[0]), &(S.symModes[0]), S.order * sizeof(dim_t));
	    return FLA_SUCCESS;
	}

	dim_t i;
	dim_t symGroupOffset = TLA_sym_group_mode_offset(S, sym_group);

	//Reorder modes to indicate the split
	S1->order = S.order;
	memcpy(&((S1->symModes)[0]), &(S.symModes[0]), S.order * sizeof(dim_t));
	dim_t symModePos[nSplit_modes];
	for(i = 0; i < nSplit_modes; i++)
	    symModePos[i] = TLA_sym_pos_of_mode(*S1, split_modes[i]);

	for(i = 0; i < nSplit_modes; i++){
	    (S1->symModes)[symModePos[i]] = (S1->symModes)[symGroupOffset + i];
	    (S1->symModes)[symGroupOffset+i] = split_modes[i];
	}
	
	//Update sym_group_info
	memcpy(&((S1->symGroupLens)[0]), &(S.symGroupLens[0]), (S.nSymGroups) * sizeof(dim_t));
	(S1->symGroupLens)[sym_group] = S.symGroupLens[sym_group] - nSplit_modes;
	for(i = S.nSymGroups - 1; i >= sym_group && i < S.nSymGroups; i--)
		(S1->symGroupLens)[i+1] = (S1->symGroupLens)[i];
	(S1->symGroupLens)[sym_group] = nSplit_modes;
	(S1->nSymGroups) = S.nSymGroups + 1;

	return FLA_SUCCESS;
}

FLA_Error TLA_update_sym_based_offset(TLA_sym S, FLA_Obj* A){
    TLA_sym* Asym = &(A->sym);
    *Asym = S;

    //Update symmetries of objects
    dim_t i, j;

    for(i = 0; i < Asym->nSymGroups; i++){
        dim_t nModes_to_split = 0;
        dim_t modes_to_split[A->order];

        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(*Asym, i);
        for(j = 0; j < (Asym->symGroupLens)[i]; j++){
            if((A->offset)[(Asym->symModes)[symGroupModeOffset + j]] == 0){
                modes_to_split[nModes_to_split] = (Asym->symModes)[symGroupModeOffset + j];
                nModes_to_split++;
            }
        }
        if(nModes_to_split == 0 || nModes_to_split == (Asym->symGroupLens)[i]){

        }else{
            TLA_split_sym_group(*Asym, nModes_to_split, modes_to_split, Asym);
            //inc i by 1 since created a new group that will not split in next iteration?
        }
    }
    return FLA_SUCCESS;
}

FLA_Error TLA_update_sym_based_lin_offset(TLA_sym S, dim_t linOffset, TLA_sym* S1){
    dim_t i, j;
    *S1 = S;

    dim_t offset[S1->order];
    memset(&(offset[0]), 0, (S1->order) * sizeof(dim_t));
    dim_t stride = 1 << (S1->order - 1);
    for(i = S1->order - 1; i < S1->order; i--){
        if(linOffset >= stride){
            offset[i] = 1;
            linOffset -= stride;
        }
        stride /= 2;
    }

    //Update symmetries of objects


    for(i = 0; i < S1->nSymGroups; i++){
        dim_t nModes_to_split = 0;
        dim_t modes_to_split[S1->order];

        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(*S1, i);
        for(j = 0; j < (S1->symGroupLens)[i]; j++){
            if((offset)[(S1->symModes)[symGroupModeOffset + j]] == 0){
                modes_to_split[nModes_to_split] = (S1->symModes)[symGroupModeOffset + j];
                nModes_to_split++;
            }
        }
        if(nModes_to_split == 0 || nModes_to_split == (S1->symGroupLens)[i]){

        }else{
            TLA_split_sym_group(*S1, nModes_to_split, modes_to_split, S1);
            //inc i by 1 since created a new group that will not split in next iteration?
        }
    }
    return FLA_SUCCESS;
}
