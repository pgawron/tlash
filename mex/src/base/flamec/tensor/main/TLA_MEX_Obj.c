#include "TLA_MEX.h"
#include "FLAME.h"

void TLA_mxa_to_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;

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

    const mxArray* size_mxa;
    const mxArray* data_mxa;

    size_mxa = mxGetFieldByNumber(mxa, 0, size_field_num);
    data_mxa = mxGetFieldByNumber(mxa, 0, data_field_num);

    if(!mxIsDouble(size_mxa)){
        mexErrMsgTxt("Input must be a tensor (size must be a double).");
    }
    if(!mxIsDouble(data_mxa)){
        mexErrMsgTxt("Input must be a tensor (data must be a double).");
    }

    //Data extracted from mxArray structures
    dim_t order;
    order = mxGetNumberOfElements(size_mxa);
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays
    double* temp;

    dim_t size[order];
    temp = mxGetPr(size_mxa);
    for(i = 0; i < order; i++)
        size[i] = (dim_t)(temp[i]);

    dim_t stride[order];
    stride[0] = 1;
    for(i = 1; i < order; i++)
        stride[i] = stride[i-1] * size[i-1];

    double* data = mxGetPr(data_mxa);
    FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size, A);
    FLA_Obj_attach_buffer_to_tensor(data, order, stride, A);
}

/****
 * CAUTION: MATLAB deals with doubles by default.  Fear for some rounding issue when converting sizes and what not
 *
 */
void TLA_mxa_to_blocked_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;

    if (!mxIsStruct(mxa))
    {
        mexErrMsgTxt("Input must be a blocked tensor (must be a structure).");
    }

    // Check that the struct has the right fields
    int flat_size_field_num = mxGetFieldNumber(mxa, "flat_size");
    int block_size_field_num = mxGetFieldNumber(mxa, "block_size");
    int data_blocks_field_num = mxGetFieldNumber(mxa, "data_blocks");

    if ((flat_size_field_num == -1) || (block_size_field_num == -1) || (data_blocks_field_num == -1))
    {
        mexErrMsgTxt("Input must be a blocked tensor (must have flat_size, block_size, and data_blocks fields).");
    }

    const mxArray* flat_size_mxa;
    const mxArray* block_size_mxa;
    const mxArray* data_blocks_mxa;

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

    int numBlocks = mxGetNumberOfElements(data_blocks_mxa);

    for(i = 0; i < numBlocks; i++){
        if(!mxIsStruct(mxGetCell(data_blocks_mxa, i))){
            mexErrMsgTxt("Data blocks must each be tensors");
        }
    }

    //Data extracted from mxArray structures
    dim_t order;
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays
    double* temp;

    order = (dim_t)mxGetNumberOfElements(flat_size_mxa);
    dim_t flat_size[order];
    dim_t blkSize[order];
    temp = mxGetPr(flat_size_mxa);
    for(i = 0; i < order; i++)
        flat_size[i] = (dim_t)(temp[i]);
    temp = mxGetPr(block_size_mxa);
    for(i = 0; i < order; i++)
        blkSize[i] = (dim_t)(temp[i]);

    double* dataBlocks[numBlocks];
    for(i = 0; i < numBlocks; i++){
        FLA_Obj curObj;
        TLA_mxa_to_tensor(mxGetCell(data_blocks_mxa, i), &curObj);

        dataBlocks[i] = (double*) FLA_malloc(curObj.base->n_elem_alloc * sizeof(double));
        memcpy(&(dataBlocks[i][0]), FLA_Obj_base_buffer(curObj), curObj.base->n_elem_alloc * sizeof(double));
    }

    dim_t stride[order];
    stride[0] = 1;
    for(i = 1; i < order; i++)
        stride[i] = (flat_size[i-1]/blkSize[i-1]) * stride[i-1];
    FLA_Obj_create_blocked_tensor_without_buffer(FLA_DOUBLE, order, flat_size, blkSize, A);
    FLA_Obj_attach_buffer_to_blocked_tensor((void**)dataBlocks, order, stride, A);
}

void TLA_mxa_to_blocked_psym_tensor(const mxArray * mxa, FLA_Obj* A){
    int i;

    if (!mxIsStruct(mxa))
    {
        mexErrMsgTxt("Input must be a blocked psym tensor (must be a structure).");
    }

    // Check that the struct has the right fields
    int flat_size_field_num = mxGetFieldNumber(mxa, "flat_size");
    int block_size_field_num = mxGetFieldNumber(mxa, "block_size");
    int data_blocks_field_num = mxGetFieldNumber(mxa, "data_blocks");
    int sym_field_num = mxGetFieldNumber(mxa, "sym");

    if ((flat_size_field_num == -1) || (block_size_field_num == -1) || (data_blocks_field_num == -1) || (sym_field_num == -1))
    {
        mexErrMsgTxt("Input must be a blocked psym tensor (must have flat_size, block_size, data_blocks, and sym fields).");
    }

    const mxArray* flat_size_mxa;
    const mxArray* block_size_mxa;
    const mxArray* data_blocks_mxa;
    const mxArray* sym_mxa;

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

    int numBlocks = mxGetNumberOfElements(data_blocks_mxa);

    for(i = 0; i < numBlocks; i++){
        if(!mxIsStruct(mxGetCell(data_blocks_mxa, i))){
            mexErrMsgTxt("Data blocks must each be tensors");
        }
    }

    //Data extracted from mxArray structures
    dim_t order;
    //Temp used as intermediate for converting MATLAB double arrays to dim_t arrays
    double* temp;

    order = (dim_t)mxGetNumberOfElements(flat_size_mxa);
    dim_t flat_size[order];
    dim_t blkSize[order];
    temp = mxGetPr(flat_size_mxa);
    for(i = 0; i < order; i++)
        flat_size[i] = (dim_t)(temp[i]);
    temp = mxGetPr(block_size_mxa);
    for(i = 0; i < order; i++)
        blkSize[i] = (dim_t)(temp[i]);

    double* dataBlocks[numBlocks];
    for(i = 0; i < numBlocks; i++){
        FLA_Obj curObj;
        TLA_mxa_to_tensor(mxGetCell(data_blocks_mxa, i), &curObj);

        dataBlocks[i] = (double*) FLA_malloc(curObj.base->n_elem_alloc * sizeof(double));
        memcpy(&(dataBlocks[i][0]), FLA_Obj_base_buffer(curObj), curObj.base->n_elem_alloc * sizeof(double));
    }

    TLA_sym sym;
    TLA_mxa_to_sym(sym_mxa, &sym);
    printf("sym\n\n");
    printf("order: %d\n", sym.order);
    printf("symgrouplens: ");
    for(i = 0; i < sym.nSymGroups; i++){
        printf(" %d", (sym.symGroupLens)[i]);
    }
    printf("\n");

    dim_t stride[order];
    stride[0] = 1;
    for(i = 1; i < order; i++)
        stride[i] = (flat_size[i-1]/blkSize[i-1]) * stride[i-1];
    FLA_Obj_create_blocked_psym_tensor_without_buffer(FLA_DOUBLE, order, flat_size, blkSize, sym, A);
    printf("Asym\n\n");
    printf("order: %d\n", A->sym.order);
    printf("symgrouplens: ");
    for(i = 0; i < (A->sym).nSymGroups; i++){
        printf(" %d", ((A->sym).symGroupLens)[i]);
    }
    printf("\n");

    FLA_Obj_attach_buffer_to_blocked_psym_tensor((void**)dataBlocks, order, stride, A);
}

void TLA_mxa_to_sym(const mxArray * mxa, TLA_sym* sym){
    dim_t i, j;

    dim_t count = 0;

    if(!mxIsCell(mxa)){
        mexErrMsgTxt("sym must be a cell");
    }

    int nSymGroups = mxGetNumberOfElements(mxa);
    (sym->nSymGroups) = nSymGroups;
    sym->order = 0;
    for(i = 0; i < nSymGroups; i++){
        mxArray* symGroup = mxGetCell(mxa, i);
        if(!mxIsDouble(symGroup)){
            mexErrMsgTxt("symGroups must be double arrays");
        }

        double* symModes_mxa = mxGetPr(symGroup);
        dim_t symGroupLen = mxGetNumberOfElements(symGroup);
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

    for(i = 0; i < A.order; i++)
        size[i] = (double)A.size[i];

    mxArray* mxa_data = mxCreateDoubleMatrix(1, A.base->n_elem_alloc, mxREAL);
    double* data = mxGetPr(mxa_data);
    double* dataBuf = (double*)FLA_Obj_base_buffer(A);

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

    for(i = 0; i < A.order; i++){
        dim_t blkDim = (((FLA_Obj*)FLA_Obj_base_buffer(A))[0]).size[i];
        dim_t blkedDim = A.size[i];
        dim_t flatDim = blkDim * blkedDim;

        blk_size[i] = (double) blkDim;
        flat_size[i] = (double) flatDim;
    }
    mxArray* mxa_blocks = mxCreateCellMatrix(1, A.base->n_elem_alloc);

    for(i = 0; i < A.base->n_elem_alloc; i++){
        FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(A))[i];
        const char* fields[] = {"size", "data"};

        mxArray *mxa_DataBlock = mxCreateStructMatrix(1, 1, 2, fields);

        TLA_tensor_to_mxa(curObj, &mxa_DataBlock);

        mxSetCell(mxa_blocks, i, mxDuplicateArray(mxa_DataBlock));
    }

    //Set up for creating blocked tensor in MATLAB
    mxArray* plhs[1];
    mxArray* prhs[3];
    prhs[0] = mxa_blocks;
    prhs[1] = mxa_flat_size;
    prhs[2] = mxa_blk_size;

    mexCallMATLAB(1, plhs, 3, prhs, "blockedtensor");
    (*mxa) = plhs[0];
}

void TLA_blocked_psym_tensor_to_mxa(FLA_Obj A, mxArray ** mxa ){
    int i;

    printf("getting fields\n");
    dim_t order = A.order;
    mxArray* mxa_order = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxArray* mxa_flat_size = mxCreateDoubleMatrix(1, A.sym.nSymGroups, mxREAL);
    mxArray* mxa_blk_size = mxCreateDoubleMatrix(1, A.sym.nSymGroups, mxREAL);
    mxArray* mxa_sym = mxCreateCellMatrix(1, A.sym.nSymGroups);

    double* order_ptr = mxGetPr(mxa_order);
    double* flat_size = mxGetPr(mxa_flat_size);
    double* blk_size = mxGetPr(mxa_blk_size);

    printf("getting sizes\n");
    order_ptr[0] = (double)order;
    for(i = 0; i < A.sym.nSymGroups; i++){
        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(A.sym, i);
        dim_t symMode = A.sym.symModes[symGroupModeOffset];
        FLA_Obj block = ((FLA_Obj*)FLA_Obj_base_buffer(A))[0];

        dim_t blkDim = block.size[symMode];
        dim_t blkedDim = A.size[symMode];
        dim_t flatDim = blkDim * blkedDim;

        blk_size[i] = (double)blkDim;
        flat_size[i] = (double)flatDim;
    }

    printf("calcing uniques\n");
    dim_t nUniqueBlocks = 1;
    TLA_sym sym = A.sym;
    printf("nSymGroups: %d", sym.nSymGroups);
    for(i = 0; i < sym.nSymGroups; i++){
        dim_t symGroupLen = sym.symGroupLens[i];
        dim_t symGroupModeOffset = TLA_sym_group_mode_offset(sym, i);
        dim_t symGroupDim = A.size[sym.symModes[symGroupModeOffset]];

        nUniqueBlocks *= binomial(symGroupLen + symGroupDim - 1, symGroupLen);
    }

    printf("nUniqueBlocks: %d\n", nUniqueBlocks);
    mxArray* mxa_blocks = mxCreateCellMatrix(1, nUniqueBlocks);

    //Works for now
    dim_t* endIndex = A.size;
    dim_t curIndex[order];
    memset(&(curIndex[0]), 0, order * sizeof(dim_t));
    dim_t update_ptr = 0;
    dim_t uniqueCount = 0;
    dim_t linCount = 0;

    FLA_Obj* buffer = (FLA_Obj*)FLA_Obj_base_buffer(A);
    printf("beginning loop\n");
    while(TRUE){
        print_array("curIndex", order, curIndex);
        printf("linCount: %d\n", linCount);
        if(buffer[linCount].isStored){
            printf("found unique\n");
            const char* fields[] = {"size", "data"};
            mxArray *mxa_DataBlock = mxCreateStructMatrix(1, 1, 2, fields);

            printf("ja?\n");
            FLA_Obj_print_matlab("to convert", buffer[linCount]);
            printf("converting tensor to mxa\n");

            TLA_tensor_to_mxa(buffer[linCount], &mxa_DataBlock);
            mxSetCell(mxa_blocks, uniqueCount, mxDuplicateArray(mxa_DataBlock));
            uniqueCount++;
        }

        printf("looping\n");
        //Loop Update
        linCount++;
        curIndex[update_ptr]++;
        while(update_ptr < order && curIndex[update_ptr] == endIndex[update_ptr]){
            update_ptr++;
            if(update_ptr < order)
                curIndex[update_ptr]++;
        }
        if(update_ptr >= order)
            break;
        printf("update_ptr %d\n", update_ptr);

        memset(&(curIndex[0]), 0, (update_ptr) * sizeof(dim_t));
        update_ptr = 0;
    }

    TLA_sym_to_mxa(A.sym, &mxa_sym);

    //Set up for creating blocked tensor in MATLAB
    mxArray* plhs[1];
    mxArray* prhs[5];
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
