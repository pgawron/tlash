#include "FLAME.h"
#include "stdio.h"

void initSymmObj(FLA_Obj* obj){
  dim_t i;
  dim_t order = 3;
  dim_t size[] = {8,8,8};
  dim_t blkSize[] = {2,2,2};
  dim_t blocked_stride[FLA_MAX_ORDER];
  TLA_sym sym;

  blocked_stride[0] = 1;
  for(i = 1; i < order; i++)
      blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / blkSize[i-1]);


  sym.order = FLA_Obj_order(*obj);
  sym.nSymGroups = 1;
  sym.symGroupLens[0] = sym.order;
  for(i = 0; i < sym.order; i++)
      (sym.symModes)[i] = i;
  FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size, blocked_stride, blkSize, sym, obj);
  FLA_Random_psym_tensor(*obj);
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
  dim_t i;

  initSymmObj(&A);

	for (i = 0; i < 64; i++){
		if ((((FLA_Obj*)A.base->buffer)[i]).base->buffer == 0) {
			printf("error on %d\n", i);
		}
	}
  printObj(A);
}
