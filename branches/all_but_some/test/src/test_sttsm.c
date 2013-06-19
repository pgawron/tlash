#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sttsm operation.\n\n");
    printf("  test_sttsm <m> <nA> <nC> <bA> <bC> <psym_temps>\n\n");
    printf("  m: order of hyper-symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  nC: mode-length of hyper-symmetric tensor C\n");
    printf("  bA: mode-length of block of A\n");
    printf("  bC: mode-length of block of C\n");
	printf("  psym_temps: Use psym temps or not\n\n");
}

void initSymmTensor(dim_t order, dim_t size[order], dim_t b, FLA_Obj* obj){
    dim_t i;
    dim_t blocked_stride[order];
	dim_t block_size[order];
    blocked_stride[0] = 1;
	block_size[0] = b;
    for(i = 1; i < order; i++){
        blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / block_size[i-1]);
		block_size[i] = b;
	}
		
    TLA_sym sym;
    sym.order = order;
    sym.nSymGroups = 1;
    sym.symGroupLens[0] = sym.order;
    for(i = 0; i < sym.order; i++)
        (sym.symModes)[i] = i;
  FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size, blocked_stride, block_size, sym, obj);
  FLA_Random_psym_tensor(*obj);
}

void initMatrix(dim_t size[2], dim_t bC, dim_t bA, FLA_Obj* obj){
  dim_t i,j;
  dim_t order = 2;
  dim_t sizeObj[] = {size[0] / bC, size[1] / bA};
  dim_t strideObj[] = {1, sizeObj[0]};
  dim_t sizeBlk[] = {bC, bA};
  dim_t strideBlk[] = {1, bC};

  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, sizeObj, obj);
  obj->base->elemtype = FLA_TENSOR;

  FLA_Obj* buf = (FLA_Obj*)FLA_malloc(sizeObj[0] * sizeObj[1] * sizeof(FLA_Obj));
  FLA_Obj_attach_buffer_to_tensor(buf, order, strideObj, obj);

  FLA_Adjust_2D_info(obj);

  for(i = 0; i < sizeObj[1]; i++)
	for(j = 0; j < sizeObj[0]; j++){
		FLA_Obj* curObj = FLA_Obj_base_buffer(*obj);
		FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeBlk, strideBlk, &(curObj[j + (sizeObj[0]*i)]));
		FLA_Random_matrix(curObj[j+ (sizeObj[0]*i)]);
		FLA_Adjust_2D_info(&(curObj[j+(sizeObj[0]*i)]));
	}
}

void setSymmTensorToZero(FLA_Obj obj){
    dim_t i,j;
	dim_t order = FLA_Obj_order(obj);
	dim_t numel = 1;
	for(i = 0; i < order; i++)
		numel *= FLA_Obj_dimsize(obj, i);

	FLA_Obj* obj_buf = (FLA_Obj*)FLA_Obj_base_buffer(obj);
	for(i = 0; i < numel; i++){
		dim_t inner_numel = 1;
		for(j = 0; j < order; j++)
			inner_numel *= FLA_Obj_dimsize(obj_buf[i], j);
		memset(&(((double*)FLA_Obj_base_buffer(obj_buf[i]))[0]), 0, inner_numel * sizeof(double));
	}
}

void test_sttsm(int m, int nA, int nC, int bA, int bC, int psym_temps, double* elapsedTime){
	dim_t i;
	dim_t aSize[m];
	for(i = 0; i < m; i++)
		aSize[i] = nA;
	dim_t bSize[] = {nC, nA};
	dim_t cSize[m];
	for(i = 0; i < m; i++)
		cSize[i] = nC;

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj A, B, C;

  initSymmTensor(m, aSize, bA, &A);

  initMatrix(bSize, bC, bA, &B);

  initSymmTensor(m, cSize, bC, &C);
//  setSymmTensorToZero(C);

	FLA_Obj_print_matlab("A", A);

	FLA_Obj_print_matlab("B", B);

	FLA_Obj_print_matlab("preC", C);
	
  double startTime = FLA_Clock();
	if(psym_temps == 1)
		FLA_Sttsm_with_psym_temps(alpha, A, beta, B, C);
	else {
		FLA_Sttsm_without_psym_temps(alpha, A, beta, B, C);
	}

  double endTime = FLA_Clock();

  *elapsedTime = endTime - startTime;


	FLA_Obj_print_matlab("C", C);

	printf("diff = C - (preC + ttm(A,{B");
	for(i = 1; i < m; i++)
		printf(",B");
	printf("}));\n");
	printf("max(diff(:))\n");
	
  FLA_Obj_blocked_free_buffer(&B);
  FLA_Obj_free_without_buffer(&B);

  FLA_Obj_blocked_psym_tensor_free_buffer(&A);
  FLA_Obj_free_without_buffer(&A);
  FLA_Obj_blocked_psym_tensor_free_buffer(&C);
  FLA_Obj_free_without_buffer(&C);
}

double Sttsm_GFlops(dim_t m, dim_t nA, dim_t nC, dim_t bA, dim_t bC, double elapsedTime){
	dim_t d, i;
	dim_t ops = 0;

	dim_t barP = nC / bC;

	for(d=0;d < m; d++){
		dim_t qTerm = 1;
		for(i = 0; i < d+1; i++)
			qTerm = qTerm * barP;
		dim_t nTerm = 1;
		for(i = 0; i < m-d; i++)
			nTerm *= nA;

		dim_t factTerm = 1;
		for(i = 0; i < d+1; i++){
			factTerm *= barP + i;
			factTerm /= d+1+i;
		}

		ops += qTerm * nTerm * factTerm;
	}
	printf("flops: %d\n", 2*ops); 
	return 2*ops/(1.e9*elapsedTime);
}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 7){
		Usage();
		FLA_Finalize();
		return 0;
	}
		

	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	const int nA = atoi(argv[++argNum]);
	const int nC = atoi(argv[++argNum]);
	const int bA = atoi(argv[++argNum]);
	const int bC = atoi(argv[++argNum]);
	const int psym_temps = atoi(argv[++argNum]);
	
	if(nA % bA != 0 || nC % bC != 0){
		printf("bA must evenly divide nA and bC must evenly divide nC\n");
		FLA_Finalize();
		return 0;
	}
	if(m <= 0 || bA <= 0 || bC < 0){
		printf("m and b must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	double elapsedTime;
	test_sttsm(m, nA, nC, bA, bC, psym_temps, &elapsedTime);

	FLA_Finalize();

	//Print out results
//	double gflops = Sttsm_GFlops(m, nA, nC, bA, bC, elapsedTime);
	printf("ARGS COMPACT %d %d %d %d %d\n", m, nA, nC, bA, bC);
    printf("TIME %.6f\n", elapsedTime);
//	printf("GFLOPS %.6f\n", gflops);
	return 0;
}
