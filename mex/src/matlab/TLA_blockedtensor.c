#include "TLA_MEX.h"
#include "FLAME.h"
void mexFunction(int nargout, mxArray * pargout[],
                 int nargin, const mxArray * pargin[]){
	int i;
	FLA_Obj A;
    if (nargin > 1)
    {
        mexErrMsgTxt("Too many input arguments.");
    }

    FLA_Init();


    //Parse input

    mexPrintf("Creating blocked tensor\n");
    TLA_mxa_to_blocked_tensor(pargin[0], &A);

    FLA_Obj_print_matlab("A", A);
    //Pass output
    mexPrintf("Passing A back\n");
    TLA_blocked_tensor_to_mxa(A, &(pargout[0]));

    //DELETE EVERYTHING
    mexPrintf("finalizing\n");
    FLA_Finalize();
}
