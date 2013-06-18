#include "TLA_MEX.h"
#include "FLAME.h"
void mexFunction(int nargout, mxArray * pargout[],
                 int nargin, const mxArray * pargin[]){
    if (nargin > 1)
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    FLA_Init();
    int i;

    //Parse input
    FLA_Obj A;
    mexPrintf("Creating blocked psym tensor\n");
    TLA_mxa_to_blocked_psym_tensor(pargin[0], &A);

    FLA_Obj_print_matlab("A", A);
    //Pass output
    mexPrintf("Passing A back\n");
    TLA_blocked_psym_tensor_to_mxa(A, &(pargout[0]));

    //DELETE EVERYTHING
    mexPrintf("finalizing\n");
    FLA_Finalize();
}
