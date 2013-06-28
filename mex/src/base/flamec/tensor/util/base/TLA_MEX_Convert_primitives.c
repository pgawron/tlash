#include "TLA_MEX.h"

void TLA_mxa_to_dim_t(const mxArray* mxa, dim_t* val){
	double* data;
    mexPrintf("mxa class %s\n", mxGetClassName(mxa));
    if(!mxIsDouble(mxa)){
        mexErrMsgTxt("Mode multiplication must be double");
    }

    data = mxGetPr(mxa);
    (*val) = (dim_t)(*data);
}
