#include "TLA_MEX.h"
#include "FLAME.h"

void mexFunction(int nargout, mxArray * pargout[],
                 int nargin, const mxArray * pargin[]){
    int i;
    FLA_Obj A, B, C;
    dim_t order;
    dim_t C_size[FLA_MAX_ORDER];
    dim_t C_blked_size[FLA_MAX_ORDER];
    dim_t C_block_size[FLA_MAX_ORDER];
    dim_t C_stride[FLA_MAX_ORDER];

    FLA_Obj alpha, beta;

    if (nargin > 2)
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    FLA_Init();


    //Parse input
    mexPrintf("Parsing A\n");
    TLA_mxa_to_blocked_psym_tensor(pargin[0], &A);
    FLA_Obj_print_matlab("A",A);

    mexPrintf("Parsing B\n");
    TLA_mxa_to_blocked_tensor(pargin[1], &B);
    FLA_Obj_print_matlab("B",B);

    order = A.order;

    for(i = 0; i < order; i++){
        C_blked_size[i] = B.size[0];
        C_block_size[i] = (((FLA_Obj*)FLA_Obj_base_buffer(B))[0]).size[0];
        C_size[i] = C_blked_size[i] * C_block_size[i];
    }


    FLA_Set_tensor_stride(order, C_blked_size, C_stride);

    mexPrintf("Creating sym C\n");
    FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, C_size, C_stride, C_block_size, A.sym, &C);
    FLA_Set_zero_tensor(C);
    FLA_Obj_print_matlab("preC",C);

    //Compute
    alpha = FLA_ONE;
    beta = FLA_ONE;

    FLA_Sttsm_with_psym_temps(alpha, A, beta, B, C);
    FLA_Obj_print_matlab("C", C);
    //Pass output
    mexPrintf("Passing C back\n");
    TLA_blocked_psym_tensor_to_mxa(C, &(pargout[0]));

    //DELETE EVERYTHING
    mexPrintf("finalizing\n");
    FLA_Finalize();
}

