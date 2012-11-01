#include "FLAME.h"
#include "stdio.h"


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

void test_sttsm(){
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
	dim_t tOrder = 3;
	dim_t tSize[] = {4, 4, 4};
	dim_t mSize[] = {8, 4};
	dim_t cOrder = 3;
	dim_t cSize[] = {8, 8, 8};
	dim_t blkSize = 2;
	//End setup parameters

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj t, m, c;

  initSymmTensor(tOrder, tSize, blkSize, &t);
  initMatrix(mSize, blkSize, &m);

  initSymmTensor(cOrder, cSize, blkSize, &c);
  setSymmTensorToZero(c);

    //check identity multiply
/*
	FLA_Obj* t_buf = (FLA_Obj*)FLA_Obj_base_buffer(t);
	((double*)t_buf[0].base->buffer)[0] = 1.00;
	((double*)t_buf[1].base->buffer)[0] = 0;
	((double*)t_buf[2].base->buffer)[0] = 0;
	((double*)t_buf[3].base->buffer)[0] = 1.00;
*/
  
	printf("t tensor\n");
	FLA_Obj_print_tensor(t);
	
	printf("m matrix\n");
	FLA_Obj_print_tensor(m);

  FLA_Sttsm(alpha, t, beta, m, c);

	printf("c tensor\n");
	FLA_Obj_print_tensor(c);
}
