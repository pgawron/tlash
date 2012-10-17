#include "FLAME.h"
#include "stdio.h"

dim_t tOrder = 3;
dim_t tSize[] = {4, 4, 4};
dim_t mOrder = 2;
dim_t mSize[] = {8, 4};
dim_t cOrder = 3;
dim_t cSize[] = {8, 8, 8};
dim_t blksize = 2;

void initSymmTensor(dim_t order, dim_t size[order], dim_t b, FLA_Obj* obj){
  FLA_Obj_create_symm_tensor_without_buffer(FLA_DOUBLE, order, size, b, obj);
  FLA_Obj_create_Random_symm_tensor_data(b, *obj);
}

void initMatrix(FLA_Obj* obj){
  FLA_Obj_create(FLA_DOUBLE, mSize[0], mSize[1], 1, mSize[0], obj);
  FLA_Random_matrix(*obj);
}

void printObj(FLA_Obj obj){
	dim_t i;
	dim_t order = FLA_Obj_order(obj);
	dim_t* idx = &((obj.base)->index[0]);
	dim_t n_elem;

	if(FLA_Obj_elemtype(obj) == FLA_TENSOR){
		printf("block idx:");
		for(i = 0; i < order; i++)
			printf(" %d", idx[i]);
		printf("\n");

		n_elem = 1;
		for(i = 0; i < order; i++)
			n_elem *= FLA_Obj_dimsize(obj, i);
		FLA_Obj* buf = FLA_Obj_base_buffer(obj);
		for(i = 0; i < n_elem; i++)
			printObj(buf[i]);
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
		for(j = 0; j < inner_numel; j++)
			inner_numel *= FLA_Obj_dimsize(obj_buf[i], j);
		memset(&(((double*)FLA_Obj_base_buffer(obj_buf[i]))[0]), 0, inner_numel * sizeof(double));
	}
}

void test_sttsm(){
  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj t, m, c;

  initSymmTensor(tOrder, tSize, blkSize);
  initMatrix(m);
  initSymmTensor(cOrder, cSize, blkSize);
  setSymmTensorToZero(c);

  FLA_Ttm_single(alpha, t, mode_mult, beta, m, c);
  double* res = (double*)((c.base)->buffer);
}
