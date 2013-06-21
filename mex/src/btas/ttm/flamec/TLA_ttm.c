#include "TLA_MEX.h"
#include "FLAME.h"

void mexFunction(int nargout, mxArray * pargout[],
                 int nargin, const mxArray * pargin[]){

	int i;
	dim_t order;
	FLA_Obj A, B, C;
	dim_t mode;
    dim_t C_size[FLA_MAX_ORDER];
    dim_t C_blked_size[FLA_MAX_ORDER];
    dim_t C_block_size[FLA_MAX_ORDER];
    dim_t C_stride[FLA_MAX_ORDER];

    FLA_Obj alpha, beta;

	if (nargin > 3)
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    FLA_Init();

    //Parse input
    mexPrintf("Parsing A\n");
    TLA_mxa_to_blocked_tensor(pargin[0], &A);
    FLA_Obj_print_matlab("A",A);

    mexPrintf("Parsing B\n");
    TLA_mxa_to_blocked_tensor(pargin[1], &B);
    FLA_Obj_print_matlab("B",B);

    mexPrintf("Parsing mode\n");
    TLA_mxa_to_dim_t(pargin[2], &mode);
    mode--;
    printf("mode: %d\n", mode);

    order = A.order;

    memcpy(&(C_blked_size[0]), &(A.size[0]), A.order * sizeof(dim_t));
    memcpy(&(C_block_size[0]), &((((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size)[0]), A.order * sizeof(dim_t));
    C_blked_size[mode] = B.size[0];
    C_block_size[mode] = (((FLA_Obj*)FLA_Obj_base_buffer(B))[0]).size[0];

    FLA_array_elemwise_product(order, C_blked_size, C_block_size, C_size);
    FLA_Set_tensor_stride(order, C_blked_size, C_stride);

    mexPrintf("Creating C\n");
    FLA_Obj_create_blocked_tensor(FLA_DOUBLE, A.order, C_size, C_stride, C_block_size, &C);
    FLA_Set_zero_tensor(C);
    FLA_Obj_print_matlab("C",C);

    //Compute
    alpha = FLA_ONE;
    beta = FLA_ONE;

    FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
    //Pass output
    mexPrintf("Passing C back\n");
    TLA_blocked_tensor_to_mxa(C, &(pargout[0]));

    //DELETE EVERYTHING
    mexPrintf("finalizing\n");
    FLA_Finalize();
}

