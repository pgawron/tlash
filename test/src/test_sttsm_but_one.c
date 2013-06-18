#include "FLAME.h"
#include "stdio.h"

void Usage() {
	printf("Test sttsm_but_one operation.\n\n");
	printf("  test_sttsm_but_one <m> <nA> <nC> <bA> <bC> <ignore_mode>\n\n");
	printf("  m: order of hyper-symmetric tensors\n");
	printf("  nA: mode-length of tensor A\n");
	printf("  nC: mode-length of hyper-symmetric tensor C\n");
	printf("  bA: mode-length of block of A\n");
	printf("  bC: mode-length of block of C\n");
	printf(
			"  ignore_mode: the mode which will be left out of the multiply\n\n");
}

void initSymmTensor(dim_t order, dim_t nA, dim_t bA, FLA_Obj* obj) {
	dim_t i;
	dim_t flat_size[order];
	dim_t blocked_stride[order];
	dim_t blocked_size[order];
	dim_t block_size[order];
	for (i = 0; i < order; i++) {
		flat_size[i] = nA;
		block_size[i] = bA;
		blocked_size[i] = nA / bA;
	}

	blocked_stride[0] = 1;
	for (i = 1; i < order; i++) {
		blocked_stride[i] = blocked_stride[i - 1] * blocked_size[i - 1];
	}

	TLA_sym sym;
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
	dim_t i, j;
	dim_t order = 2;
	dim_t sizeObj[] = { nC / bC, nA / bA };
	dim_t strideObj[] = { 1, sizeObj[0] };
	dim_t sizeBlk[] = { bC, bA };
	dim_t strideBlk[] = { 1, bC };

	FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, sizeObj, obj);
	obj->base->elemtype = FLA_TENSOR;

	FLA_Obj* buf = (FLA_Obj*) FLA_malloc(
			sizeObj[0] * sizeObj[1] * sizeof(FLA_Obj));
	FLA_Obj_attach_buffer_to_tensor(buf, order, strideObj, obj);

	FLA_Adjust_2D_info(obj);

	for (i = 0; i < sizeObj[1]; i++)
		for (j = 0; j < sizeObj[0]; j++) {
			FLA_Obj* curObj = FLA_Obj_base_buffer(*obj);
			FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeBlk, strideBlk,
					&(curObj[j + (sizeObj[0] * i)]));
			FLA_Random_matrix(curObj[j + (sizeObj[0] * i)]);
			FLA_Adjust_2D_info(&(curObj[j + (sizeObj[0] * i)]));
		}
}

void initOutputTensor(dim_t order, dim_t nA, dim_t nC, dim_t bA, dim_t bC,
		dim_t ignore_mode, FLA_Obj* obj) {

	dim_t i;
	dim_t flat_size[order];
	dim_t blocked_size[order];
	dim_t blocked_stride[order];
	dim_t block_size[order];

	for (i = 0; i < order; i++) {
		flat_size[i] = nC;
		block_size[i] = bC;
	}
	flat_size[ignore_mode] = nA;
	block_size[ignore_mode] = bA;

	for (i = 0; i < order; i++) {
		blocked_size[i] = flat_size[i] / block_size[i];
	}

	blocked_stride[0] = 1;
	for (i = 1; i < order; i++) {
		blocked_stride[i] = blocked_stride[i - 1] * blocked_size[i - 1];
	}

	TLA_sym sym;
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

	initSymmTensor(m, nA, bA, &A);
	initMatrix(nC, nA, bC, bA, &B);
	initOutputTensor(m, nA, nC, bA, bC, ignore_mode, &C);

	FLA_Obj_print_matlab("A", A);
	FLA_Obj_print_matlab("B", B);
	FLA_Obj_print_matlab("preC", C);

	double startTime = FLA_Clock();

	FLA_Sttsm_but_one(alpha, A, ignore_mode, beta, B, C);

	double endTime = FLA_Clock();

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

	FLA_Obj_blocked_free_buffer(&B);
	FLA_Obj_free_without_buffer(&B);

	FLA_Obj_blocked_psym_tensor_free_buffer(&A);
	FLA_Obj_free_without_buffer(&A);
	FLA_Obj_blocked_psym_tensor_free_buffer(&C);
	FLA_Obj_free_without_buffer(&C);
}
int main(int argc, char* argv[]) {
	FLA_Init();

	if (argc < 7) {
		Usage();
		FLA_Finalize();
		return 0;
	}

	int argNum = 0;
	const dim_t m = atoi(argv[++argNum]);
	const dim_t nA = atoi(argv[++argNum]);
	const dim_t nC = atoi(argv[++argNum]);
	const dim_t bA = atoi(argv[++argNum]);
	const dim_t bC = atoi(argv[++argNum]);
	const dim_t ignore_mode = atoi(argv[++argNum]);

	if (nA % bA != 0 || nC % bC != 0) {
		printf("bA must evenly divide nA and bC must evenly divide nC\n");
		FLA_Finalize();
		return 0;
	}
	if (m <= 0 || bA <= 0 || bC < 0) {
		printf("m and b must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	double elapsedTime;
	test_sttsm(m, nA, nC, bA, bC, ignore_mode, &elapsedTime);

	FLA_Finalize();

	//Print out results
//	double gflops = Sttsm_GFlops(m, nA, nC, bA, bC, elapsedTime);
	printf("ARGS COMPACT %d %d %d %d %d\n", m, nA, nC, bA, bC);
	printf("TIME %.6f\n", elapsedTime);
//	printf("GFLOPS %.6f\n", gflops);
	return 0;
}
