#include "FLAME.h"
#include "stdio.h"

void Usage() {
	printf("Test sttsm_but_one operation.\n\n");
	printf("  test_sttsm_but_one <m> <nA> <nC> <bA> <bC> <ignore_mode>\n\n");
	printf("  m: order of symmetric tensors\n");
	printf("  nA: mode-length of tensor A\n");
	printf("  nC: mode-length of symmetric tensor C\n");
	printf("  bA: mode-length of block of A\n");
	printf("  bC: mode-length of block of C\n");
	printf(
			"  ignore_mode: the mode which will be left out of the multiply\n\n");
}

void initSymmTensor(dim_t order, dim_t nA, dim_t bA, FLA_Obj* obj) {
	dim_t i;
	dim_t flat_size[FLA_MAX_ORDER];
	dim_t blocked_stride[FLA_MAX_ORDER];
	dim_t blocked_size[FLA_MAX_ORDER];
	dim_t block_size[FLA_MAX_ORDER];
	TLA_sym sym;

	for (i = 0; i < order; i++) {
		flat_size[i] = nA;
		block_size[i] = bA;
	}

	FLA_array_elemwise_quotient(order, flat_size, block_size, blocked_size);
	FLA_Set_tensor_stride(order, blocked_size, blocked_stride);

	sym.order = order;
	sym.nSymGroups = 1;
	sym.symGroupLens[0] = sym.order;
	for (i = 0; i < sym.order; i++)
		(sym.symModes)[i] = i;
	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size,
			blocked_stride, block_size, sym, obj);
	FLA_Random_psym_tensor(*obj);
}

void initMatrix(dim_t nC, dim_t nA, dim_t bC, dim_t bA, FLA_Obj* obj) {
	dim_t order = 2;
	dim_t size[] = {nC, nA};
	dim_t sizeObj[] = {nC / bC, nA / bA};
	dim_t strideObj[] = { 1, sizeObj[0] };
	dim_t sizeBlk[] = { bC, bA };

	FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size, strideObj, sizeBlk, obj);
	FLA_Random_tensor(*obj);
}

void initOutputTensor(dim_t order, dim_t nA, dim_t nC, dim_t bA, dim_t bC,
		dim_t ignore_mode, FLA_Obj* obj)
{
	dim_t i;
	dim_t flat_size[FLA_MAX_ORDER];
	dim_t blocked_size[FLA_MAX_ORDER];
	dim_t blocked_stride[FLA_MAX_ORDER];
	dim_t block_size[FLA_MAX_ORDER];
	TLA_sym sym;

	for (i = 0; i < order; i++) {
		flat_size[i] = nC;
		block_size[i] = bC;
	}
	flat_size[ignore_mode] = nA;
	block_size[ignore_mode] = bA;

	FLA_array_elemwise_quotient(order, flat_size, block_size, blocked_size);
	FLA_Set_tensor_stride(order, blocked_size, blocked_stride);

	sym.order = order;
	sym.nSymGroups = 2;
	sym.symGroupLens[0] = 1;
	sym.symGroupLens[1] = sym.order - 1;
	dim_t symIndex = 1;
	for (i = 0; i < order; i++) {
		if (i == ignore_mode)
			(sym.symModes)[0] = i;
		else
			(sym.symModes[symIndex++]) = i;
	}
	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size,
			blocked_stride, block_size, sym, obj);
	FLA_Random_psym_tensor(*obj);
}

void test_sttsm(int m, int nA, int nC, int bA, int bC, dim_t ignore_mode,
		double* elapsedTime) {
	dim_t i;

	FLA_Obj alpha = FLA_ONE;
	FLA_Obj beta = FLA_ONE;

	FLA_Obj A, B, C;

	double startTime;
	double endTime;

	initSymmTensor(m, nA, bA, &A);
	initMatrix(nC, nA, bC, bA, &B);
	initOutputTensor(m, nA, nC, bA, bC, ignore_mode, &C);

	FLA_Obj_print_matlab("A", A);
	FLA_Obj_print_matlab("B", B);
	FLA_Obj_print_matlab("preC", C);

	startTime = FLA_Clock();

	FLA_Sttsm_but_one(alpha, A, ignore_mode, beta, B, C);

	endTime = FLA_Clock();

	*elapsedTime = endTime - startTime;

	FLA_Obj_print_matlab("C", C);
	printf("diff = C - (preC + ttm(A,{B");
	for (i = 1; i < m - 1; i++)
		printf(",B");
	printf("},[");
	for (i = 0; i < m; i++) {
		if (i != ignore_mode)
			printf(" %d", i + 1);
	}
	printf("]));\n");
	printf("max(diff(:))\n");

	FLA_Obj_blocked_tensor_free_buffer(&B);
	FLA_Obj_free_without_buffer(&B);

	FLA_Obj_blocked_psym_tensor_free_buffer(&A);
	FLA_Obj_free_without_buffer(&A);
	FLA_Obj_blocked_psym_tensor_free_buffer(&C);
	FLA_Obj_free_without_buffer(&C);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* nC, dim_t* bA, dim_t* bC, dim_t* ignore_mode){
    int argNum = 0;

    if(argc != 7){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *nC = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);
    *bC = atoi(argv[++argNum]);
    *ignore_mode = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t nC, dim_t bA, dim_t bC, dim_t ignore_mode){

	if(order <= 0 || nA <= 0 || nC <= 0 || bA <= 0 || bC <= 0){
		printf("m, nA, nC, bA, and bC must be greater than 0\n");
		return FLA_FAILURE;
	}

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		return FLA_FAILURE;
	}

	if(nC % bC != 0){
		printf("bC must evenly divide nC\n");
		return FLA_FAILURE;
	}

	if(ignore_mode >= order){
		printf("ignore_mode must be less than order\n");
		return FLA_FAILURE;
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]) {
	dim_t order;
	dim_t nA;
	dim_t nC;
	dim_t bA;
	dim_t bC;
	dim_t ignore_mode;

	double elapsedTime;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &nC, &bA, &bC, &ignore_mode) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, nC, bA, bC, ignore_mode) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_sttsm(order, nA, nC, bA, bC, ignore_mode, &elapsedTime);

	FLA_Finalize();

	//Print out results
	printf("ARGS COMPACT %d %d %d %d %d\n", order, nA, nC, bA, bC);
	printf("TIME %.6f\n", elapsedTime);
	return 0;
}
