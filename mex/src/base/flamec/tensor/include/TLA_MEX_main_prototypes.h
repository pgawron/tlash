//TLASH -> MEX
void TLA_tensor_to_mxa(FLA_Obj A, mxArray ** mxa);
void TLA_blocked_tensor_to_mxa(FLA_Obj A, mxArray ** mxa );

//MEX -> TLASH
void TLA_mxa_to_tensor(const mxArray * mxa, FLA_Obj* A);
void TLA_mxa_to_blocked_tensor(const mxArray * mxa, FLA_Obj* A);
