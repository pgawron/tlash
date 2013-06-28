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

FLA_Error FLA_Obj_create_tensor( FLA_Datatype datatype, dim_t order, const dim_t size[], const dim_t stride[], FLA_Obj *obj)
{

  FLA_Obj_create_tensor_ext( datatype, FLA_SCALAR, order, size, size, stride, obj );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_psym_tensor(FLA_Datatype datatype, dim_t order, const dim_t size[], const dim_t stride[], TLA_sym sym, FLA_Obj *obj){
  FLA_Obj_create_tensor( datatype, order, size, stride, obj);

  //Update symmetries
  obj->sym = sym;

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_blocked_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, const dim_t flat_size[], const dim_t size_inner[], const dim_t stride[], const dim_t blk_size[], FLA_Obj *obj )
{
  dim_t i;
  dim_t size_obj[FLA_MAX_ORDER];
  dim_t stride_obj[FLA_MAX_ORDER];
  dim_t stride_blk[FLA_MAX_ORDER];
  FLA_Obj* buf;

  //First set up the obj to store the FLA_Objs
  //Set the blocked-object size & stride
  FLA_array_elemwise_quotient(order, flat_size, blk_size, size_obj);
  FLA_Set_tensor_stride(order, size_obj, stride_obj);

  FLA_Obj_create_tensor_ext( datatype, FLA_TENSOR, order, size_obj, size_obj, stride_obj, obj);

  //Create the blocks of obj
  buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
  FLA_Set_tensor_stride(order, blk_size, stride_blk);

  for(i = 0; i < FLA_Obj_num_elem_alloc(*obj); i++){
	FLA_Obj_create_tensor( datatype, order, blk_size, stride_blk, &(buf[i]));
  }

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, const dim_t size[], const dim_t size_inner[], const dim_t stride[], FLA_Obj *obj )
{
  dim_t i;
  dim_t nSecondDim;

  //First set the 2-D object info
  nSecondDim = FLA_array_product(order - 1, &(size[1]));
  FLA_Obj_create_ext( datatype, elemtype, size[0], nSecondDim, size[0], nSecondDim, stride[0], stride[0]*size[0], obj );

  //Update the tensor info
  //View object
  obj->order = order;
  memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

  //Base object
  obj->base->order = order;
  memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->size_inner)[0]), &(size_inner[0]), order * sizeof( dim_t ) );
  memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof( dim_t ) );
  memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );
  obj->base->n_elem_alloc = size[0] * nSecondDim;

  //View metadata (permutation & isStored)
  obj->isStored = TRUE;
  for(i = 0; i < order; i++)
	obj->permutation[i] = i;

  //Update symmetries (Assumed dense)
  TLA_Sym_init_nonsymmetric(order, &(obj->sym));

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_tensor_without_buffer( FLA_Datatype datatype, dim_t order, const dim_t size[], FLA_Obj *obj ){
	dim_t i;
	dim_t nSecondDim;

	//First set 2-D object info
	nSecondDim = FLA_array_product(order - 1, &(size[1]));
	FLA_Obj_create_without_buffer(datatype, size[0], nSecondDim, obj);

	//Update tensor info
	//View object
	obj->order = order;
    memcpy(&((obj->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->offset)[0]), 0, order * sizeof( dim_t ) );

    //Base object
    obj->base->order = order;
    memcpy(&((obj->base->size)[0]), &(size[0]), order * sizeof( dim_t ) );
    memcpy(&((obj->base->size_inner)[0]), &(size[0]), order * sizeof( dim_t ) );
    memset(&((obj->base->index)[0]), 0, order * sizeof( dim_t ) );
    memset(&((obj->base->stride)[0]), 0, order * sizeof( dim_t ) );
	obj->base->n_elem_alloc = size[0] * nSecondDim;

	//View metadata (permutation & isStored)
	obj->isStored = FALSE;
	for(i = 0; i < order; i++)
		(obj->permutation)[i] = i;

    //Update symmetries
	TLA_Sym_init_nonsymmetric(order, &(obj->sym));

	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, const dim_t stride[], FLA_Obj *obj ){
	//Omitting some things attach_buffer does because not sure how to extend yet

	//Attach and adjust meta-data
	obj->base->buffer = buffer;
	memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));
	obj->isStored = TRUE;

	//Potentially unneeded
	FLA_Adjust_2D_info(obj);
	
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_blocked_tensor_free_buffer( FLA_Obj *obj)
{
	if(FLA_Obj_elemtype(*obj) == FLA_TENSOR || FLA_Obj_elemtype(*obj) == FLA_MATRIX){
		dim_t i;
		FLA_Obj* buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);
		for(i = 0; i < FLA_Obj_num_elem_alloc(*obj); i++){
			FLA_Obj_free_buffer(&(buf[i]));
			FLA_Obj_free_without_buffer(&(buf[i]));
		}
		FLA_Obj_free_buffer(obj);
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_blocked_psym_tensor_free_buffer( FLA_Obj *obj)
{
	dim_t order = FLA_Obj_order(*obj);
	const dim_t* stride = obj->base->stride;
	FLA_Obj* buf;

	dim_t update_ptr;
	const dim_t* endIndex = obj->size;
	dim_t curIndex[FLA_MAX_ORDER];
	dim_t linIndex;

	//Init obj data
	buf = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

	//Init loop data
	update_ptr = order - 1;
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));

	//Loop over all blocked indices, if block is stored, free it
	while(TRUE){
		dim_t i;
		dim_t isUnique;
		//For this index, find which block it is
		linIndex = FLA_TIndex_to_LinIndex(order, curIndex, stride);

		isUnique = buf[linIndex].isStored;
		if(isUnique){
			FLA_Obj_free_buffer(&(buf[linIndex]));
			FLA_Obj_free_without_buffer(&(buf[linIndex]));
		}
		
		//Loop Update
		curIndex[update_ptr]++;
		while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
			update_ptr--;
			if(update_ptr < order)
				curIndex[update_ptr]++;
		}
		if(update_ptr >= order)
			break;
		for(i = update_ptr+1; i < order; i++)
			curIndex[i] = 0;
		update_ptr = order - 1;
	}

	//Free alloc'd data
	FLA_Obj_free_buffer(obj);
	return FLA_SUCCESS;
}


//NOTE: This function doesn't set the offset array of the FLA_View
FLA_Error FLA_Obj_create_blocked_tensor_without_buffer(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blk_size[], FLA_Obj *obj){
    dim_t i;
    dim_t blked_size[FLA_MAX_ORDER];
    dim_t nTBlks;
    FLA_Obj* t_blks;
    dim_t stride_obj[FLA_MAX_ORDER];

    //Determine blocked size of tensor (size of tensor whose elements are the blocks)
    FLA_array_elemwise_quotient(order, flat_size, blk_size, blked_size);
    nTBlks = FLA_array_product(order, blked_size);

    //Create the tensor blocks
    t_blks = (FLA_Obj*)FLA_malloc(nTBlks * sizeof(FLA_Obj));
    for(i = 0; i < nTBlks; i++){
        FLA_Obj *curObj = &(t_blks[i]);
        FLA_Obj_create_tensor_without_buffer( datatype, order, blk_size, curObj);
    }

    //Buffer of tensor blocks created, set the main obj to represent this hierarchy
    //Create object
    FLA_Obj_create_tensor_without_buffer( datatype, order, blked_size, obj );
    obj->base->elemtype = FLA_TENSOR;

    //Attach blocks to object
    FLA_Set_tensor_stride(order, blked_size, stride_obj);
    FLA_Obj_attach_buffer_to_tensor( t_blks, order, stride_obj, obj );

    return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_blocked_psym_tensor_without_buffer(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blk_size[], TLA_sym sym, FLA_Obj *obj){
    dim_t i;
    dim_t blked_size[FLA_MAX_ORDER];
    dim_t nTBlks;
    FLA_Obj* t_blks;

    dim_t curIndex[FLA_MAX_ORDER];
    dim_t updateIndex;
    dim_t objLinIndex;

    dim_t stride_obj[FLA_MAX_ORDER];

//    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
//        FLA_Obj_create_blocked_sym_tensor_without_buffer_check( datatype, order, size, blkSize, obj );

    //Determine blocked size of tensor
    FLA_array_elemwise_quotient(order, flat_size, blk_size, blked_size);
    nTBlks = FLA_array_product(order, blked_size);

    //Create the tensor blocks
    t_blks = (FLA_Obj*)FLA_malloc(nTBlks * sizeof(FLA_Obj));

    memset(curIndex, 0, order * sizeof(dim_t));
    updateIndex = 0;
    objLinIndex = 0;

    //Loop over block-indexes in tensor.  Create blocks and set offsets
    while(TRUE){
        //Set up the FLA_Obj at this index
        FLA_Obj *curObj = &(t_blks[objLinIndex]);
        FLA_Obj_create_tensor_without_buffer( datatype, order, blk_size, curObj);

        //Set the offset array
        memset(&((curObj->offset)[0]), 0, order * sizeof(dim_t));

        //Loop update
        curIndex[updateIndex]++;
        objLinIndex++;
        //If we hit the end, loop until we find the index to update
        while(updateIndex < order && curIndex[updateIndex] == blked_size[updateIndex]){
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
    //TODO: See if datatype can be swapped with FLA_TENSOR
    FLA_Obj_create_tensor_without_buffer( datatype, order, blked_size, obj);
    obj->sym = sym;
    obj->base->elemtype = FLA_TENSOR;


    FLA_Set_tensor_stride(order, blked_size, stride_obj);
    FLA_Obj_attach_buffer_to_tensor(t_blks, order, stride_obj, obj);

    return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_blocked_tensor(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blocked_stride[], const dim_t blk_size[], FLA_Obj *obj){
    dim_t i;
    dim_t blked_size[FLA_MAX_ORDER];
    dim_t nBlocks;
    dim_t nBlockElems;

    void** dataBuffers;

    //Set up the sym struct
    TLA_sym sym;
    TLA_Sym_init_nonsymmetric(order, &sym);

    //First set up the hierarchy without buffers
    FLA_Obj_create_blocked_psym_tensor_without_buffer(datatype, order, flat_size, blk_size, sym, obj);

    //Create each block for tensor
    //NOTE: probably already know how many blocks there are...
    FLA_array_elemwise_quotient(order, flat_size, blk_size, blked_size);
    nBlocks = FLA_array_product(order, blked_size);
    nBlockElems = FLA_array_product(order, blk_size);

    //Create dataBuffers for each block
    dataBuffers = (void**) FLA_malloc( nBlocks * sizeof(void*) );
    for (i = 0; i < nBlocks; i++)
    {
        dataBuffers[i] = (double*) FLA_malloc( nBlockElems * sizeof(double) );
        memset( &(((double* )dataBuffers[i])[0]), 0, nBlockElems * sizeof(double) );
    }

    //Attach empty buffers to the tensor blocks
    FLA_Obj_attach_buffer_to_blocked_tensor( dataBuffers, order, blocked_stride,
            obj );

    //Free local arrays
    FLA_free(dataBuffers);
    return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create_blocked_psym_tensor(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blocked_stride[], const dim_t blk_size[], TLA_sym sym, FLA_Obj *obj){
	dim_t i;
	dim_t nBlockElems;
	dim_t blked_size[FLA_MAX_ORDER];
	dim_t nUniques = 1;
	dim_t modeOffset = 0;

	void** dataBuffers;

	//First set up the hierarchy without buffers
	FLA_Obj_create_blocked_psym_tensor_without_buffer(datatype, order, flat_size, blk_size, sym, obj);
	
	//Set up the data buffers for psym tensor
	nBlockElems = FLA_array_product(order, blk_size);

	//Determine number of unique blocks in tensor
	FLA_array_elemwise_quotient(order, flat_size, blk_size, blked_size);

	for(i = 0; i < (obj->sym).nSymGroups; i++){
	    dim_t mode = (obj->sym).symModes[modeOffset];
	    dim_t symGroupLen = (obj->sym).symGroupLens[i];
	    dim_t blkedDim = blked_size[mode];

		nUniques *= binomial(symGroupLen + blkedDim - 1, symGroupLen);
		modeOffset += symGroupLen;
	}

	//Create data arrays for each block
	dataBuffers = (void**)FLA_malloc(nUniques * sizeof(void*));
	for(i = 0; i < nUniques; i++){
		dataBuffers[i] = (double*)FLA_malloc(nBlockElems * sizeof(double));
		memset(&(((double*)dataBuffers[i])[0]), 0, nBlockElems * sizeof(double));
	}

	//Attach empty buffers to the sym tensor
	FLA_Obj_attach_buffer_to_blocked_psym_tensor(dataBuffers, order, blocked_stride, obj);

	//Free local data
	FLA_free(dataBuffers);
	return FLA_SUCCESS;
}


FLA_Error FLA_Obj_attach_buffer_to_blocked_tensor( void *buffer[], dim_t order, const dim_t stride[], FLA_Obj *obj ){
    dim_t i, j;
    const dim_t* size_obj = obj->size;
    //const dim_t* stride_obj = obj->base->stride;
    FLA_Obj *buffer_obj;

    dim_t nBlocks;

    //Get needed info from object
    buffer_obj = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

    nBlocks = FLA_array_product(order, size_obj);

    //Attach buffer to block and adjust stride/permutation
    for(i = 0; i < nBlocks; i++){
        (buffer_obj[i].base)->buffer = buffer[i];
        FLA_Set_tensor_stride(order, (buffer_obj[i].base)->size, (buffer_obj[i].base)->stride);

        for(j = 0; j < order; j++)
            ((buffer_obj[i]).permutation)[j] = j;
        (buffer_obj[i]).isStored = TRUE;
    }

    //Omitting some things attach_buffer does because not sure how to extend yet
    memcpy(&((obj->base->stride)[0]), &(stride[0]), order * sizeof(dim_t));

    return FLA_SUCCESS;
}

FLA_Error FLA_Obj_attach_buffer_to_blocked_psym_tensor( void *buffer[], dim_t order, const dim_t stride[], FLA_Obj *obj ){
	dim_t i;
	//FLA_Obj-related data
    const dim_t* size_obj = obj->size;
    const dim_t* stride_obj = obj->base->stride;
	FLA_Obj *buffer_obj;

	//Loop-related data
	dim_t updateIndex;
	dim_t curIndex[FLA_MAX_ORDER];

	//Attach buffer-related data
	dim_t countBuffer;
	dim_t objLinIndex;
	dim_t sortedIndex[FLA_MAX_ORDER];
	dim_t permutation[FLA_MAX_ORDER];
	dim_t ipermutation[FLA_MAX_ORDER];
	dim_t isUniqueIndex;

	//Set needed info of obj
    buffer_obj = (FLA_Obj*)FLA_Obj_base_buffer(*obj);

    //Initialize loop information
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	countBuffer = 0;
	objLinIndex = 0;

	//Loop over indices
	//If we hit a unique, set its buffer to the next block in the list
	//Otherwise, find the unique index and point buffer to unique buffer (adjusting strides)
	updateIndex = order - 1;
	while(TRUE){
		objLinIndex = FLA_TIndex_to_LinIndex(order, stride_obj, curIndex);

		//For this index, get the needed unique-index information
		isUniqueIndex = FLA_get_unique_info(obj->sym, curIndex, sortedIndex, permutation, ipermutation);
		
		if(isUniqueIndex){
			//Index is unique:
			//Attach block, update metadata (stride, permutation, & isStored
			(buffer_obj[objLinIndex].base)->buffer = buffer[countBuffer];
			FLA_Set_tensor_stride(order, (buffer_obj[objLinIndex].base)->size, (buffer_obj[objLinIndex].base)->stride);

			for(i = 0; i < order; i++)
			    ((buffer_obj[objLinIndex]).permutation)[i] = i;
			(buffer_obj[objLinIndex]).isStored = TRUE;
			//Ensure next unique index uses next block
			countBuffer++;
		}else{
			//Index is non-unique:
			//Find unique linear index, point buffer to unique linear index
			//Adjust meta-data (permutation & isStored)
			dim_t uniqueLinIndex;
			uniqueLinIndex = FLA_TIndex_to_LinIndex(order, stride_obj, sortedIndex);

			//point this non-unique FLA_Obj to the correct base
			//WARNING: HACK
			FLA_free(buffer_obj[objLinIndex].base);
			(buffer_obj[objLinIndex]).base = (buffer_obj[uniqueLinIndex]).base;

			memcpy(&(((buffer_obj[objLinIndex]).permutation)[0]), &(ipermutation[0]), order * sizeof(dim_t));
			(buffer_obj[objLinIndex]).isStored = FALSE;
		}

		//Loop update
		curIndex[updateIndex]++;
		//If we hit the end, loop until we find the index to update
		while(updateIndex < order && curIndex[updateIndex] == size_obj[updateIndex]){
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

	return FLA_SUCCESS;
}

////////////////////////////////
//More of Util functions follow
////////////////////////////////

//Note: Only splits modes within the same symmetric group...
FLA_Error TLA_split_sym_group(TLA_sym S, dim_t nSplit_modes, const dim_t split_modes[], TLA_sym* S1){
	dim_t i;
	dim_t sym_group;
	dim_t symGroupOffset;

	dim_t split_modes_copied;
	dim_t index_to_fill;

	dim_t nextGroupOffset;
	dim_t nRemainingModes;

	//0 modes to split, done
    if(nSplit_modes == 0){
        *S1 = S;
        return FLA_SUCCESS;
    }

    //If group to split is size 1, done
    sym_group = TLA_sym_group_of_mode(S, split_modes[0]);
    if(TLA_sym_group_size(S, sym_group) == 1){
    	*S1 = S;
	    return FLA_SUCCESS;
	}

    //Otherwise, do more work
    S1->order = S.order;

	symGroupOffset = TLA_sym_group_mode_offset(S, sym_group);

	//Reorder modes to indicate the split
	//Copy all modes before the split to output
	memcpy(&((S1->symModes)[0]), &(S.symModes[0]), symGroupOffset * sizeof(dim_t));
	
	//Copy the modes to be split next
	memcpy(&((S1->symModes)[symGroupOffset]), &(split_modes[0]), nSplit_modes * sizeof(dim_t));

	//Copy the non-split modes in the sym group
	split_modes_copied = 0;
	index_to_fill = symGroupOffset + nSplit_modes;
	for(i = 0; i < S.symGroupLens[sym_group]; i++){
		if(split_modes_copied < nSplit_modes && (S.symModes[symGroupOffset + i] == split_modes[split_modes_copied])){
			split_modes_copied++;
		}else{
			(S1->symModes)[index_to_fill] = S.symModes[symGroupOffset + i];
			index_to_fill++;
		}
	}

	//Copy remaining modes in other sym groups
	nextGroupOffset = symGroupOffset + S.symGroupLens[sym_group];
	nRemainingModes = S.order - nextGroupOffset;
	memcpy(&((S1->symModes)[nextGroupOffset]), &(S.symModes[nextGroupOffset]), (nRemainingModes) * sizeof(dim_t));

	
	//Update sym_group_info
	//Update number of groups
	(S1->nSymGroups) = S.nSymGroups + 1;

	//Copy lengths of groups before group in question
	memcpy(&((S1->symGroupLens)[0]), &(S.symGroupLens[0]), (sym_group) * sizeof(dim_t));
	//Copy lengths of groups after group in question
	//leaving space for the update to the split group
	memcpy(&((S1->symGroupLens)[sym_group + 2]), &(S.symGroupLens[sym_group + 1]), (S.nSymGroups - (sym_group + 1)) * sizeof(dim_t));
	//Fill in entries for split group
	S1->symGroupLens[sym_group] = nSplit_modes;
	S1->symGroupLens[sym_group+1] = S.symGroupLens[sym_group] - nSplit_modes;


	return FLA_SUCCESS;
}

//Might be incorrect
//Updates the symmetry of an object based on its offset index & passed symmetry.
//Result should be an object with symmetry being subset of S
FLA_Error TLA_update_sym_based_offset(TLA_sym S, FLA_Obj* A){
    dim_t i, j;
	//Adjust A symmetry to be same as S
    TLA_sym* Asym = &(A->sym);
    *Asym = S;
    
    //Based on the offset index, break symmetry of modes in groups
    //For example (0 1 2) at [0,1,1] -> (0) (1 2)

    for(i = 0; i < Asym->nSymGroups; i++){
        dim_t nModes_to_split = 0;
        //NOTE: Switch to S.order?
        //NOTE: move modes_to_split out of the loop?
        dim_t modes_to_split[FLA_MAX_ORDER];

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
