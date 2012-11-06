#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sttsm operation.\n\n");
    printf("  test_sttsm <m> <nA> <nC><b>\n\n");
    printf("  m: order of hyper-symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  nC: mode-length of hyper-symmetric tensor C\n");
    printf("  b: mode-length of block\n\n");
}

void initSymmTensor(dim_t order, dim_t size[order], dim_t b, FLA_Obj* obj){
  FLA_Obj_create_symm_tensor_without_buffer(FLA_DOUBLE, order, size, b, obj);
  FLA_Obj_create_Random_symm_tensor_data(b, *obj);
}

void initMatrix(dim_t size[2], dim_t b, FLA_Obj* obj){
  dim_t i,j;
  dim_t order = 2;
  dim_t sizeObj[] = {size[0] / b, size[1] / b};
  dim_t strideObj[] = {1, sizeObj[0]};
  dim_t sizeBlk[] = {b, b};
  dim_t strideBlk[] = {1, b};

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

void test_sttsm(int m, int nA, int nC, int b){
	//Setup parameters
/*
	dim_t stOrder = 3;
	dim_t stSize[] = {4, 4, 4};
	dim_t smOrder = 2;
	dim_t smSize[] = {8, 4};
	dim_t smStride[] = {1, 8};
	dim_t scOrder = 3;
	dim_t scSize[] = {8, 8, 8};
	dim_t sblkSize = 2;
*/
	dim_t i;
	dim_t aSize[m];
	for(i = 0; i < m; i++)
		aSize[i] = nA;
	dim_t bSize[] = {nC, nA};
	dim_t cSize[m];
	for(i = 0; i < m; i++)
		cSize[i] = nC;
	dim_t blkSize = b;
	//End setup parameters

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj A, B, C;

  initSymmTensor(m, aSize, blkSize, &A);
  initMatrix(bSize, blkSize, &B);

  initSymmTensor(m, cSize, blkSize, &C);
  setSymmTensorToZero(C);

    //check identity multiply
/*
	FLA_Obj* t_buf = (FLA_Obj*)FLA_Obj_base_buffer(t);
	((double*)t_buf[0].base->buffer)[0] = 1.00;
	((double*)t_buf[1].base->buffer)[0] = 0;
	((double*)t_buf[2].base->buffer)[0] = 0;
	((double*)t_buf[3].base->buffer)[0] = 1.00;
*/
  
//	printf("t tensor\n");
//	FLA_Obj_print_flat_tensor(t);
	
//	printf("m matrix\n");
//	FLA_Obj_print_flat_tensor(m);

  FLA_Sttsm(alpha, A, beta, B, C);

//	printf("c tensor\n");
//	FLA_Obj_print_flat_tensor(c);
}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 5){
		Usage();
		FLA_Finalize();
		return 0;
	}
		

	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	const int nA = atoi(argv[++argNum]);
	const int nC = atoi(argv[++argNum]);
	const int b = atoi(argv[++argNum]);

	if(nA % b != 0 || nC % b != 0){
		printf("b must evenly divide nA and nC\n");
		FLA_Finalize();
		return 0;
	}
	if(m <= 0 || b <= 0){
		printf("m and b must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	test_sttsm(m, nA, nC, b);

	FLA_Finalize();
	return 0;
}
