#include "FLAME.h"
#include "stdio.h"

dim_t stOrder = 3;
dim_t stSize[] = {4, 4, 4};
dim_t smOrder = 2;
dim_t smSize[] = {8, 4};
dim_t smStride[] = {1, 8};
dim_t scOrder = 3;
dim_t scSize[] = {8, 8, 8};
dim_t sblkSize = 2;

void initSymmTensor(dim_t order, dim_t size[order], dim_t b, FLA_Obj* obj){
  FLA_Obj_create_symm_tensor_without_buffer(FLA_DOUBLE, order, size, b, obj);
  FLA_Obj_create_Random_symm_tensor_data(b, *obj);
}

void initMatrix(FLA_Obj* obj){
  dim_t i,j;
  dim_t sizeObj[] = {smSize[0] / sblkSize, smSize[1] / sblkSize};
  dim_t strideObj[] = {1, sizeObj[0]};
  dim_t sizeBlk[] = {sblkSize, sblkSize};
  dim_t strideBlk[] = {1, sblkSize};

  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, smOrder, sizeObj, obj);
  obj->base->elemtype = FLA_TENSOR;

  FLA_Obj* buf = (FLA_Obj*)FLA_malloc(sizeObj[0] * sizeObj[1] * sizeof(FLA_Obj));
  FLA_Obj_attach_buffer_to_tensor(buf, smOrder, strideObj, obj);

  FLA_Adjust_2D_info(obj);

  for(i = 0; i < sizeObj[1]; i++)
	for(j = 0; j < sizeObj[0]; j++){
		FLA_Obj* curObj = FLA_Obj_base_buffer(*obj);
		FLA_Obj_create_tensor(FLA_DOUBLE, smOrder, sizeBlk, strideBlk, &(curObj[j + (sizeObj[0]*i)]));
		FLA_Random_matrix(curObj[j+ (sizeObj[0]*i)]);
		FLA_Adjust_2D_info(&(curObj[j+(sizeObj[0]*i)]));
	}
}

void printSymmObj(FLA_Obj obj){
	dim_t i;
	dim_t order = FLA_Obj_order(obj);
	dim_t* idx = &((obj.base)->index[0]);
	dim_t n_elem;

	if(FLA_Obj_elemtype(obj) == FLA_TENSOR || FLA_Obj_elemtype(obj) == FLA_MATRIX){
		printf("block idx:");
		for(i = 0; i < order; i++)
			printf(" %d", idx[i]);
		printf("\n");

		n_elem = 1;
		for(i = 0; i < order; i++)
			n_elem *= FLA_Obj_dimsize(obj, i);
		FLA_Obj* buf = FLA_Obj_base_buffer(obj);
		for(i = 0; i < n_elem; i++)
			printSymmObj(buf[i]);
	}
	else{
		printf("data:");
		n_elem = 1;
		for(i = 0; i < order; i++)
			n_elem *= FLA_Obj_dimsize(obj, i);
		double* buf = (double*)FLA_Obj_base_buffer(obj);
		for(i = 0; i < n_elem; i++)
			printf(" %.3f", buf[i]);
		printf("\n");
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

void test_sttsm(){
  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj t, m, c;

  initSymmTensor(stOrder, stSize, sblkSize, &t);
  initMatrix(&m);

  initSymmTensor(scOrder, scSize, sblkSize, &c);
  setSymmTensorToZero(c);

	printf("t tensor\n");
	FLA_Obj_print_tensor(t);
	
	printf("m matrix\n");
	FLA_Obj_print_tensor(m);
	printf("c tensor before\n");
	FLA_Obj_print_tensor(c);
	

  FLA_Sttsm(alpha, t, beta, m, c);

	printf("c tensor\n");
	FLA_Obj_print_tensor(c);
}
