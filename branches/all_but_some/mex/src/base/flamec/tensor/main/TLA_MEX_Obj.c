#include "TLA_MEX.h"
#include "FLAME.h"

void TLA_mxa_to_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;
    dim_t order;
    const mxArray* size_mxa;
    const mxArray* data_mxa;

    double* temp;
    double* data;

    dim_t size[FLA_MAX_ORDER];
    dim_t stride[FLA_MAX_ORDER];

    if(!mxIsStruct(mxa)){
        mexErrMsgTxt("Input must be a tensor (must be a structure).");
    }

    // Check that the struct has the right fields
    int size_field_num = mxGetFieldNumber(mxa, "size");
    int data_field_num = mxGetFieldNumber(mxa, "data");

    if ((size_field_num == -1) || (data_field_num == -1))
    {
        mexErrMsgTxt("Input must be a tensor (must have size and data fields).");
    }

    size_mxa = mxGetFieldByNumber(mxa, 0, size_field_num);
    data_mxa = mxGetFieldByNumber(mxa, 0, data_field_num);

    if(!mxIsDouble(size_mxa)){
        mexErrMsgTxt("Input must be a tensor (size must be a double).");
    }
    if(!mxIsDouble(data_mxa)){
        mexErrMsgTxt("Input must be a tensor (data must be a double).");
    }

    //Data extracted from mxArray structures
    order = mxGetNumberOfElements(size_mxa);
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays

    temp = mxGetPr(size_mxa);
    for(i = 0; i < order; i++)
        size[i] = (dim_t)(temp[i]);

    FLA_Set_tensor_stride(order, size, stride);

    data = mxGetPr(data_mxa);
    FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size, A);
    FLA_Obj_attach_buffer_to_tensor(data, order, stride, A);
}

/****
 * CAUTION: MATLAB deals with doubles by default.  Fear for some rounding issue when converting sizes and what not
 *
 */
void TLA_mxa_to_blocked_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;
    int flat_size_field_num;
    int block_size_field_num;
    int data_blocks_field_num;

    const mxArray* flat_size_mxa;
    const mxArray* block_size_mxa;
    const mxArray* data_blocks_mxa;

    int numBlocks;

    //Data extracted from mxArray structures
    dim_t order;
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays
    double* temp;

    dim_t flat_size[FLA_MAX_ORDER];
    dim_t blkSize[FLA_MAX_ORDER];
    dim_t blked_size[FLA_MAX_ORDER];
    double** dataBlocks;
    dim_t stride[FLA_MAX_ORDER];

    if (!mxIsStruct(mxa))
    {
        mexErrMsgTxt("Input must be a blocked tensor (must be a structure).");
    }

    // Check that the struct has the right fields
    flat_size_field_num = mxGetFieldNumber(mxa, "flat_size");
    block_size_field_num = mxGetFieldNumber(mxa, "block_size");
    data_blocks_field_num = mxGetFieldNumber(mxa, "data_blocks");

    if ((flat_size_field_num == -1) || (block_size_field_num == -1) || (data_blocks_field_num == -1))
    {
        mexErrMsgTxt("Input must be a blocked tensor (must have flat_size, block_size, and data_blocks fields).");
    }

    flat_size_mxa = mxGetFieldByNumber(mxa, 0, flat_size_field_num);
    block_size_mxa = mxGetFieldByNumber(mxa, 0, block_size_field_num);
    data_blocks_mxa = mxGetFieldByNumber(mxa, 0, data_blocks_field_num);

    if(!mxIsDouble(flat_size_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (flat_size must be a double).");
    }
    if(!mxIsDouble(block_size_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (block_size must be a double).");
    }
    if(!mxIsCell(data_blocks_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (data_blocks must be a cell).");
    }

    numBlocks = mxGetNumberOfElements(data_blocks_mxa);

    for(i = 0; i < numBlocks; i++){
        if(!mxIsStruct(mxGetCell(data_blocks_mxa, i))){
            mexErrMsgTxt("Data blocks must each be tensors");
        }
    }


    order = (dim_t)mxGetNumberOfElements(flat_size_mxa);

    temp = mxGetPr(flat_size_mxa);
    for(i = 0; i < order; i++)
        flat_size[i] = (dim_t)(temp[i]);
    temp = mxGetPr(block_size_mxa);
    for(i = 0; i < order; i++)
        blkSize[i] = (dim_t)(temp[i]);

    dataBlocks = (double**)FLA_malloc(numBlocks * sizeof(double*));
    for(i = 0; i < numBlocks; i++){
        FLA_Obj curObj;
        TLA_mxa_to_tensor(mxGetCell(data_blocks_mxa, i), &curObj);

        dataBlocks[i] = (double*) FLA_malloc(curObj.base->n_elem_alloc * sizeof(double));
        memcpy(&(dataBlocks[i][0]), FLA_Obj_base_buffer(curObj), curObj.base->n_elem_alloc * sizeof(double));
    }

    FLA_array_elemwise_quotient(order, flat_size, blkSize, blked_size);
    FLA_Set_tensor_stride(order, blked_size, stride);

    FLA_Obj_create_blocked_tensor_without_buffer(FLA_DOUBLE, order, flat_size, blkSize, A);
    FLA_Obj_attach_buffer_to_blocked_tensor((void**)dataBlocks, order, stride, A);
}

void TLA_mxa_to_blocked_psym_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;
    int flat_size_field_num;
    int block_size_field_num;
    int data_blocks_field_num;
    int sym_field_num;

    const mxArray* flat_size_mxa;
    const mxArray* block_size_mxa;
    const mxArray* data_blocks_mxa;
    const mxArray* sym_mxa;

    int numBlocks;

    //Data extracted from mxArray structures
    dim_t order;
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays
    double* temp;

    dim_t flat_size[FLA_MAX_ORDER];
    dim_t blkSize[FLA_MAX_ORDER];
    dim_t blked_size[FLA_MAX_ORDER];
    double** dataBlocks;
    dim_t stride[FLA_MAX_ORDER];

    TLA_sym sym;

    if (!mxIsStruct(mxa))
    {
        mexErrMsgTxt("Input must be a blocked psym tensor (must be a structure).");
    }

    // Check that the struct has the right fields
    flat_size_field_num = mxGetFieldNumber(mxa, "flat_size");
    block_size_field_num = mxGetFieldNumber(mxa, "block_size");
    data_blocks_field_num = mxGetFieldNumber(mxa, "data_blocks");
    sym_field_num = mxGetFieldNumber(mxa, "sym");

    if ((flat_size_field_num == -1) || (block_size_field_num == -1) || (data_blocks_field_num == -1) || (sym_field_num == -1))
    {
        mexErrMsgTxt("Input must be a blocked psym tensor (must have flat_size, block_size, data_blocks, and sym fields).");
    }


    flat_size_mxa = mxGetFieldByNumber(mxa, 0, flat_size_field_num);
    block_size_mxa = mxGetFieldByNumber(mxa, 0, block_size_field_num);
    data_blocks_mxa = mxGetFieldByNumber(mxa, 0, data_blocks_field_num);
    sym_mxa = mxGetFieldByNumber(mxa, 0, sym_field_num);

    if(!mxIsDouble(flat_size_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (flat_size must be a double).");
    }
    if(!mxIsDouble(block_size_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (block_size must be a double).");
    }
    if(!mxIsCell(data_blocks_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (data_blocks must be a cell).");
    }

    if(!mxIsCell(sym_mxa)){
        mexErrMsgTxt("Input must be a blocked tensor (sym must be a cell).");
    }

    numBlocks = mxGetNumberOfElements(data_blocks_mxa);

    for(i = 0; i < numBlocks; i++){
        if(!mxIsStruct(mxGetCell(data_blocks_mxa, i))){
            mexErrMsgTxt("Data blocks must each be tensors");
        }
    }

    order = (dim_t)mxGetNumberOfElements(flat_size_mxa);
    temp = mxGetPr(flat_size_mxa);
    for(i = 0; i < order; i++)
        flat_size[i] = (dim_t)(temp[i]);
    temp = mxGetPr(block_size_mxa);
    for(i = 0; i < order; i++)
        blkSize[i] = (dim_t)(temp[i]);

    dataBlocks = (double**)FLA_malloc(numBlocks * sizeof(double*));
    for(i = 0; i < numBlocks; i++){
        FLA_Obj curObj;
        TLA_mxa_to_tensor(mxGetCell(data_blocks_mxa, i), &curObj);

        dataBlocks[i] = (double*) FLA_malloc(curObj.base->n_elem_alloc * sizeof(double));
        memcpy(&(dataBlocks[i][0]), FLA_Obj_base_buffer(curObj), curObj.base->n_elem_alloc * sizeof(double));
    }


    TLA_mxa_to_sym(sym_mxa, &sym);

    FLA_array_elemwise_quotient(order, flat_size, blkSize, blked_size);
    FLA_Set_tensor_stride(order, blked_size, stride);

    FLA_Obj_create_blocked_psym_tensor_without_buffer(FLA_DOUBLE, order, flat_size, blkSize, sym, A);
    FLA_Obj_attach_buffer_to_blocked_psym_tensor((void**)dataBlocks, order, stride, A);
}

void TLA_mxa_to_sym(const mxArray * mxa, TLA_sym* sym){
    dim_t i, j;
    dim_t count = 0;
    int nSymGroups;

    if(!mxIsCell(mxa)){
        mexErrMsgTxt("sym must be a cell");
    }

    nSymGroups = mxGetNumberOfElements(mxa);
    (sym->nSymGroups) = nSymGroups;
    sym->order = 0;
    for(i = 0; i < nSymGroups; i++){
    	double* symModes_mxa;
    	dim_t symGroupLen;
        mxArray* symGroup = mxGetCell(mxa, i);
        if(!mxIsDouble(symGroup)){
            mexErrMsgTxt("symGroups must be double arrays");
        }

        symModes_mxa = mxGetPr(symGroup);
        symGroupLen = mxGetNumberOfElements(symGroup);
        for(j = 0; j < symGroupLen; j++){
            //Subtract 1 because MATLAB 1-indexes
            (sym->symModes)[count++] = (dim_t)(symModes_mxa[j]) - 1;
        }
        (sym->symGroupLens)[i] = symGroupLen;
        (sym->order) += symGroupLen;
    }

}

void TLA_tensor_to_mxa(FLA_Obj A, mxArray ** mxa){
    int i;
    mxArray* mxa_size = mxCreateDoubleMatrix(1, A.order, mxREAL);
    double* size = mxGetPr(mxa_size);
    mxArray* mxa_data;
    double* data;
    double* dataBuf;

    for(i = 0; i < A.order; i++)
        size[i] = (double)A.size[i];

    mxa_data = mxCreateDoubleMatrix(1, A.base->n_elem_alloc, mxREAL);
    data = mxGetPr(mxa_data);
    dataBuf = (double*)FLA_Obj_base_buffer(A);

    memcpy(&(data[0]), &(dataBuf[0]), A.base->n_elem_alloc * sizeof(double));

    mxArray* plhs[1];
    mxArray* prhs[2];
    prhs[0] = mxa_data;
    prhs[1] = mxa_size;

    mexCallMATLAB(1, plhs, 2, prhs, "tensor");
    (*mxa) = plhs[0];
}

void TLA_blocked_tensor_to_mxa(FLA_Obj A, mxArray ** mxa ){
    int i;

    mxArray* mxa_flat_size = mxCreateDoubleMatrix(1, A.order, mxREAL);
    mxArray* mxa_blk_size = mxCreateDoubleMatrix(1, A.order, mxREAL);
    double* flat_size = mxGetPr(mxa_flat_size);
    double* blk_size = mxGetPr(mxa_blk_size);

    mxArray* mxa_blocks;

    mxArray* plhs[1];
    mxArray* prhs[3];

    for(i = 0; i < A.order; i++){
        dim_t blkDim = (((FLA_Obj*)FLA_Obj_base_buffer(A))[0]).size[i];
        dim_t blkedDim = A.size[i];
        dim_t flatDim = blkDim * blkedDim;

        blk_size[i] = (double) blkDim;
        flat_size[i] = (double) flatDim;
    }
    mxa_blocks = mxCreateCellMatrix(1, A.base->n_elem_alloc);

    for(i = 0; i < A.base->n_elem_alloc; i++){
        FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(A))[i];
        const char* fields[] = {"size", "data"};
        mxArray *mxa_DataBlock = mxCreateStructMatrix(1, 1, 2, fields);

        TLA_tensor_to_mxa(curObj, &mxa_DataBlock);

        mxSetCell(mxa_blocks, i, mxDuplicateArray(mxa_DataBlock));
    }

    //Set up for creating blocked tensor in MATLAB

    prhs[0] = mxa_blocks;
    prhs[1] = mxa_flat_size;
    prhs[2] = mxa_blk_size;

    mexCallMATLAB(1, plhs, 3, prhs, "blockedtensor");
    (*mxa) = plhs[0];
}

void TLA_blocked_psym_tensor_to_mxa(FLA_Obj A, mxArray ** mxa ){
    int i, j;

    dim_t order = A.order;
    mxArray* mxa_order = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray* mxa_flat_size = mxCreateDoubleMatrix(1, A.sym.nSymGroups, mxREAL);
    mxArray* mxa_blk_size = mxCreateDoubleMatrix(1, A.sym.nSymGroups, mxREAL);
    mxArray* mxa_sym = mxCreateCellMatrix(1, A.sym.nSymGroups);

    double* order_ptr = mxGetPr(mxa_order);
    double* flat_size = mxGetPr(mxa_flat_size);
    double* blk_size = mxGetPr(mxa_blk_size);

    dim_t nUniqueBlocks;
    TLA_sym sym;

    mxArray* mxa_blocks;

    dim_t* endIndex;
    dim_t curIndex[FLA_MAX_ORDER];
    dim_t update_ptr = 0;
    dim_t uniqueCount = 0;

    /**
     *  All this for psym
     */
    dim_t permutation[FLA_MAX_ORDER];
    dim_t sortedIndex[FLA_MAX_ORDER];
    dim_t* stride_obj = FLA_Obj_stride(A);
    dim_t objLinIndex;
    dim_t nSymGroups = A.sym.nSymGroups;
    dim_t symGroupLens[FLA_MAX_ORDER];
    dim_t symModes[FLA_MAX_ORDER];

    FLA_Paired_Sort index_pairs[FLA_MAX_ORDER];
    dim_t orderedSymModes[FLA_MAX_ORDER];
    /**
     * End necessary for psym
     */
    mxArray* plhs[1];
    mxArray* prhs[5];


    order_ptr[0] = (double)order;
    for(i = 0; i < A.sym.nSymGroups; i++){
        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(A.sym, i);
        dim_t symMode = A.sym.symModes[symGroupModeOffset];
        dim_t blkDim;
        dim_t blkedDim;
        dim_t flatDim;

        FLA_Obj block = ((FLA_Obj*)FLA_Obj_base_buffer(A))[0];

        blkDim = block.size[symMode];
        blkedDim = A.size[symMode];
        flatDim = blkDim * blkedDim;

        blk_size[i] = (double)blkDim;
        flat_size[i] = (double)flatDim;
    }

    nUniqueBlocks = 1;
    sym = A.sym;
    for(i = 0; i < sym.nSymGroups; i++){
        dim_t symGroupLen = sym.symGroupLens[i];
        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(sym, i);
        dim_t symGroupDim = A.size[sym.symModes[symGroupModeOffset]];

        nUniqueBlocks *= binomial(symGroupLen + symGroupDim - 1, symGroupLen);
    }

    mxa_blocks = mxCreateCellMatrix(1, nUniqueBlocks);

    //Works for now
    endIndex = A.size;
    memset(&(curIndex[0]), 0, order * sizeof(dim_t));



    memcpy(&(symGroupLens[0]), &(A.sym.symGroupLens[0]), nSymGroups * sizeof(dim_t));
    memcpy(&(symModes[0]), &(A.sym.symModes[0]), order* sizeof(dim_t));


    FLA_Obj* buffer = (FLA_Obj*)FLA_Obj_base_buffer(A);
    while(TRUE){
    	dim_t modeOffset = 0;
		dim_t uniqueIndex = TRUE;
		dim_t count = 0;

		objLinIndex = FLA_TIndex_to_LinIndex(order, stride_obj, curIndex);

		for(i = 0; i < nSymGroups; i++){
			for(j = 0; j < symGroupLens[i]; j++){
				orderedSymModes[j+modeOffset] = symModes[j+modeOffset];
			}
			qsort(&(orderedSymModes[modeOffset]), symGroupLens[i], sizeof(dim_t), compare_dim_t);

			for(j = 0; j < symGroupLens[i]; j++){
				index_pairs[j].index = orderedSymModes[j+modeOffset];
				index_pairs[j].val = curIndex[orderedSymModes[j+modeOffset]];
			}
			qsort(index_pairs, symGroupLens[i], sizeof(FLA_Paired_Sort), compare_pairwise_sort);


			for(j = 0; j < symGroupLens[i]; j++){
				permutation[orderedSymModes[j+modeOffset]] = index_pairs[j].index;
				sortedIndex[orderedSymModes[j+modeOffset]] = index_pairs[j].val;
			}

			modeOffset += symGroupLens[i];
		}

		//Check if this is unique or not

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

        if(buffer[objLinIndex].isStored){
            const char* fields[] = {"size", "data"};
            mxArray *mxa_DataBlock = mxCreateStructMatrix(1, 1, 2, fields);

            print_array("unique index found", order, curIndex);
            FLA_Obj_print_matlab("convertBlock", buffer[objLinIndex]);
            TLA_tensor_to_mxa(buffer[objLinIndex], &mxa_DataBlock);
            mxSetCell(mxa_blocks, uniqueCount, mxDuplicateArray(mxa_DataBlock));
            uniqueCount++;
        }

        //Loop Update
        curIndex[update_ptr]++;
        while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
            update_ptr++;
            if(update_ptr < order)
                curIndex[update_ptr]++;
        }
        if(update_ptr >= order)
            break;

        memset(&(curIndex[0]), 0, (update_ptr) * sizeof(dim_t));
        update_ptr = 0;
    }

    TLA_sym_to_mxa(A.sym, &mxa_sym);

    //Set up for creating blocked tensor in MATLAB
    prhs[0] = mxa_blocks;
    prhs[1] = mxa_order;
    prhs[2] = mxa_flat_size;
    prhs[3] = mxa_blk_size;
    prhs[4] = mxa_sym;

    mexCallMATLAB(1, plhs, 5, prhs, "blockedpsymtensor");
    (*mxa) = plhs[0];
}

void TLA_sym_to_mxa(TLA_sym sym, mxArray** mxa){
    dim_t i, j;
    dim_t nSymGroups = sym.nSymGroups;
    dim_t count = 0;
    for(i = 0; i < nSymGroups; i++){
        dim_t symGroupLen = sym.symGroupLens[i];
        mxArray* symGroup = mxCreateDoubleMatrix(1,symGroupLen, mxREAL);
        double* symModes_mxa = mxGetPr(symGroup);

        for(j = 0; j < symGroupLen; j++){
            //Add 1 because MATLAB 1-indexes
            symModes_mxa[j] = (double)(sym.symModes[count++] + 1);
        }

        mxSetCell(*mxa, i, mxDuplicateArray(symGroup));
    }
}
