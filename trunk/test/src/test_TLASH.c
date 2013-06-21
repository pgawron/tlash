#include "FLAME.h"
#include "stdio.h"

void initObj(FLA_Obj* obj){
  dim_t i;
  dim_t order = 3;
  dim_t size[] = {3,4,5};
  dim_t stride[] = {1, 3, 12};
  dim_t nData = 1;
  int* data;
  for(i = 0; i < order; i++)
    nData *= size[i];

  data = (int*)FLA_malloc(nData * sizeof(int));
  for(i = 0; i < nData; i++)
  	data[i] = i;
  FLA_Obj_create_tensor_without_buffer(FLA_INT, order, size, obj);
  FLA_Obj_attach_buffer_to_tensor(data, order, stride, obj);
}

void printView(char* msg, FLA_Obj obj){
  dim_t i;
  dim_t order   = FLA_Obj_order(obj);
  dim_t* size   = FLA_Obj_size(obj);
  dim_t* stride = FLA_Obj_stride(obj);
  printf("View: %s\n", msg);
  printf("order: %d\n", order);
  printf("size: ");
  for(i = 0; i < order; i++)
    printf("%d ", size[i]);
  printf("\n");
  printf("stride: ");
  for(i = 0; i < order; i++)
    printf("%d ", stride[i]);
  printf("\n");
  printf("\n");
}

void repart_test(){
  FLA_Obj A;
  FLA_Obj A0, A1;
  FLA_Obj A00, A10, A20;

  dim_t* A00_buf;
  dim_t* A10_buf;
  dim_t* A20_buf;

  initObj(&A);

  printf("Hello world!\n");
  printView("A", A);



  FLA_Part_1xmode2(A, &A0, 
					  &A1,
				   1, 2, FLA_TOP);

  printView("A0", A0);
  printView("A1", A1);



  FLA_Repart_1xmode2_to_1xmode3(A0, &A00,
								    &A10,
								A1, &A20,
								1, 2, FLA_TOP);

  A00_buf = ((dim_t*)FLA_Obj_buffer_at_view(A00));
  A10_buf = ((dim_t*)FLA_Obj_buffer_at_view(A10));
  A20_buf = ((dim_t*)FLA_Obj_buffer_at_view(A20));

  printf("%d %d %d\n", A00_buf[0], A10_buf[0], A20_buf[0]);

  printView("A00", A00);
  printView("A10", A10);
  printView("A20", A20);

  FLA_Cont_with_1xmode3_to_1xmode2(&A0,A00,
								       A10,
							       &A1,A20,
								   1, FLA_TOP);

  printView("A0", A0);
  printView("A1", A1);

  FLA_Merge_1xmode2(A0,
                    A1, &A, 1);

  printView("A", A);
}

void permute_test(){
  FLA_Obj A;
  FLA_Obj A0, A1;
  FLA_Obj A00, A01,  A10, A11;

  FLA_Obj A000, A010,    A001, A011,
          A100, A110,    A101, A111;

  FLA_Obj blk_B[8];

  initObj(&A);


  FLA_Part_1xmode2(A, &A0, &A1, 1, 2, FLA_TOP);

  FLA_Part_1xmode2(A0, &A00, &A01, 2, 2, FLA_TOP);
  FLA_Part_1xmode2(A1, &A10, &A11, 2, 2, FLA_TOP);

  FLA_Part_1xmode2(A00, &A000, &A001, 3, 2, FLA_TOP);
  FLA_Part_1xmode2(A01, &A010, &A011, 3, 2, FLA_TOP);
  FLA_Part_1xmode2(A10, &A100, &A101, 3, 2, FLA_TOP);
  FLA_Part_1xmode2(A11, &A110, &A111, 3, 2, FLA_TOP);

  printView("A000", A000);
  FLA_Obj_show("data:", A000, "%d", "");
  printView("A001", A001);
  FLA_Obj_show("data:", A001, "%d", "");
  printView("A010", A010);
  FLA_Obj_show("data:", A010, "%d", "");
  printView("A011", A011);
  FLA_Obj_show("data:", A011, "%d", "");
  printView("A100", A100);
  FLA_Obj_show("data:", A100, "%d", "");
  printView("A101", A101);
  FLA_Obj_show("data:", A101, "%d", "");
  printView("A110", A110);
  FLA_Obj_show("data:", A110, "%d", "");
  printView("A111", A111);
  FLA_Obj_show("data:", A111, "%d", "");


  //FLA_Permute(permutation, A.order, blk_A_size, blk_A, &B_order, &B_size, &blk_B);
  
  printView("B000", blk_B[0]);
  FLA_Obj_show("data:", blk_B[0], "%d", "");
  printView("B000", blk_B[1]);
  FLA_Obj_show("data:", blk_B[1], "%d", "");
  printView("B000", blk_B[2]);
  FLA_Obj_show("data:", blk_B[2], "%d", "");
  printView("B000", blk_B[3]);
  FLA_Obj_show("data:", blk_B[3], "%d", "");
  printView("B000", blk_B[4]);
  FLA_Obj_show("data:", blk_B[4], "%d", "");
  printView("B000", blk_B[5]);
  FLA_Obj_show("data:", blk_B[5], "%d", "");
  printView("B000", blk_B[6]);
  FLA_Obj_show("data:", blk_B[6], "%d", "");
  printView("B000", blk_B[7]);
  FLA_Obj_show("data:", blk_B[7], "%d", "");
}

void test_tlash(){
  repart_test();
  permute_test();
}
