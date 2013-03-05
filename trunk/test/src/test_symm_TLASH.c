#include "FLAME.h"
#include "stdio.h"

void initSymmObj(FLA_Obj* obj){
  dim_t i;
  dim_t order = 3;
  dim_t size[] = {8,8,8};
  dim_t b = 2;
//  dim_t stride[] = {1, 2, 4};
  dim_t nData = 1;

  for(i = 0; i < order; i++)
    nData *= size[i];

  int* data = FLA_malloc(nData * sizeof(int));
  for(i = 0; i < nData; i++)
  	data[i] = i;
  FLA_Obj_create_symm_tensor_without_buffer(FLA_INT, order, size, b, obj);
  FLA_Obj_create_Random_symm_tensor_data(b, *obj);
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

void test_symm_tlash(){
  FLA_Obj A;
  initSymmObj(&A);
	dim_t i;
	for (i = 0; i < 64; i++){
		if ((((FLA_Obj*)A.base->buffer)[i]).base->buffer == 0) {
			printf("error on %d\n", i);
		}
	}
  printObj(A);
}
