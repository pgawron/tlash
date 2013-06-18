#include "TLA_MEX.h"
#include "FLAME.h"

void mexFunction(int nargout, mxArray * pargout[],
                 int nargin, const mxArray * pargin[]){
    if (nargin > 2)
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    FLA_Init();
    int i;

    //Parse input
    mexPrintf("Parsing A\n");
    FLA_Obj A;
    TLA_mxa_to_blocked_psym_tensor(pargin[0], &A);
    FLA_Obj_print_matlab("A",A);

    mexPrintf("Parsing B\n");
    FLA_Obj B;
    TLA_mxa_to_blocked_tensor(pargin[1], &B);
    FLA_Obj_print_matlab("B",B);

//    mexPrintf("Parsing mode\n");
//    dim_t mode;
//    TLA_mxa_to_dim_t(pargin[2], &mode);
//    mode--;
//    printf("mode: %d\n", mode);

    FLA_Obj C;
    dim_t order = A.order;
    dim_t C_size[order];
    dim_t C_blked_size[order];
    dim_t C_block_size[order];

    for(i = 0; i < order; i++){
        C_blked_size[i] = B.size[0];
        C_block_size[i] = (((FLA_Obj*)FLA_Obj_base_buffer(B))[0]).size[0];
        C_size[i] = C_blked_size[i] * C_block_size[i];
    }

    dim_t C_stride[order];
    C_stride[0] = 1;
    for(i = 1; i < order; i++)
        C_stride[i] = C_stride[i-1] * C_blked_size[i-1];

    mexPrintf("Creating sym C\n");
    FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, C_size, C_stride, C_block_size, A.sym, &C);
    FLA_Set_zero_tensor(C);
    FLA_Obj_print_matlab("preC",C);

    //Compute
    FLA_Obj alpha = FLA_ONE;
    FLA_Obj beta = FLA_ONE;

    FLA_Sttsm_with_psym_temps(alpha, A, beta, B, C);
    FLA_Obj_print_matlab("C", C);
    //Pass output
    mexPrintf("Passing C back\n");
    TLA_blocked_psym_tensor_to_mxa(C, &(pargout[0]));

    //DELETE EVERYTHING
    mexPrintf("finalizing\n");
    FLA_Finalize();
}

