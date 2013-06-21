#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sttsm operation.\n\n");
    printf("  test_sttsm <m> <nA> <nC> <bA> <bC> <psym_temps>\n\n");
    printf("  m: order of symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  nC: mode-length of symmetric tensor C\n");
    printf("  bA: mode-length of block of A\n");
    printf("  bC: mode-length of block of C\n");
	printf("  psym_temps: Use psym temps or not\n\n");
}

void initSymmTensor(dim_t order, dim_t size[], dim_t b, FLA_Obj* obj){
    dim_t i;
    dim_t blocked_stride[FLA_MAX_ORDER];
	dim_t block_size[FLA_MAX_ORDER];
	TLA_sym sym;

    blocked_stride[0] = 1;
	block_size[0] = b;
    for(i = 1; i < order; i++){
        blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / block_size[i-1]);
		block_size[i] = b;
	}
		

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
  FLA_Obj* buf;

  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, sizeObj, obj);
  obj->base->elemtype = FLA_TENSOR;

  buf = (FLA_Obj*)FLA_malloc(sizeObj[0] * sizeObj[1] * sizeof(FLA_Obj));
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
	FLA_Obj* obj_buf;

	for(i = 0; i < order; i++)
		numel *= FLA_Obj_dimsize(obj, i);

	obj_buf = (FLA_Obj*)FLA_Obj_base_buffer(obj);
	for(i = 0; i < numel; i++){
		dim_t inner_numel = 1;
		for(j = 0; j < order; j++)
			inner_numel *= FLA_Obj_dimsize(obj_buf[i], j);
		memset(&(((double*)FLA_Obj_base_buffer(obj_buf[i]))[0]), 0, inner_numel * sizeof(double));
	}
}

void test_sttsm(int m, int nA, int nC, int bA, int bC, int psym_temps, double* elapsedTime){
	dim_t i;
	dim_t aSize[FLA_MAX_ORDER];
	dim_t bSize[] = {nC, nA};
	dim_t cSize[FLA_MAX_ORDER];
	FLA_Obj alpha = FLA_ONE;
	FLA_Obj beta = FLA_ONE;
	FLA_Obj A, B, C;
	double startTime, endTime;

	for(i = 0; i < m; i++)
		aSize[i] = nA;
	for(i = 0; i < m; i++)
		cSize[i] = nC;

  initSymmTensor(m, aSize, bA, &A);

  initMatrix(bSize, bC, bA, &B);

  initSymmTensor(m, cSize, bC, &C);
//  setSymmTensorToZero(C);

	FLA_Obj_print_matlab("A", A);

	FLA_Obj_print_matlab("B", B);

	FLA_Obj_print_matlab("preC", C);
	
  startTime = FLA_Clock();
	if(psym_temps == 1)
		FLA_Sttsm_with_psym_temps(alpha, A, beta, B, C);
	else {
		FLA_Sttsm_without_psym_temps(alpha, A, beta, B, C);
	}

  endTime = FLA_Clock();

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
		dim_t nTerm = 1;
		dim_t factTerm = 1;
		for(i = 0; i < d+1; i++)
			qTerm = qTerm * barP;

		for(i = 0; i < m-d; i++)
			nTerm *= nA;

		for(i = 0; i < d+1; i++){
			factTerm *= barP + i;
			factTerm /= d+1+i;
		}

		ops += qTerm * nTerm * factTerm;
	}
	printf("flops: %d\n", 2*ops); 
	return 2*ops/(1.e9*elapsedTime);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* nC, dim_t* bA, dim_t* bC, dim_t* psym_temps){
    int argNum = 0;

    if(argc != 7){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *nC = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);
    *bC = atoi(argv[++argNum]);
    *psym_temps = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t nC, dim_t bA, dim_t bC, dim_t psym_temps){

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

	if(psym_temps >= 2){
		printf("psym_temps must be an element of {0,1}\n");
		return FLA_FAILURE;
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t nC;
	dim_t bA;
	dim_t bC;
	dim_t psym_temps;
	double elapsedTime;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &nC, &bA, &bC, &psym_temps) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, nC, bA, bC, psym_temps) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_sttsm(order, nA, nC, bA, bC, psym_temps, &elapsedTime);

	FLA_Finalize();

	//Print out results
//	double gflops = Sttsm_GFlops(m, nA, nC, bA, bC, elapsedTime);
	printf("ARGS COMPACT %d %d %d %d %d\n", order, nA, nC, bA, bC);
    printf("TIME %.6f\n", elapsedTime);
//	printf("GFLOPS %.6f\n", gflops);
	return 0;
}
