#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test TLASH Permute operation.\n\n");
    printf("  tensor_permute <m> <nA> <bA> <perm0> <perm1> ...\n\n");
    printf("  m: order of tensor\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  permK: permutation index K\n");
}

void initTensor(dim_t order, dim_t nA, FLA_Obj* A){
	dim_t i;	
	dim_t size[order];
	dim_t stride[order];
	stride[0] = 1;

	for(i = 0; i < order; i++)
		size[i] = nA;

	for(i = 1; i < order; i++)
		stride[i] = stride[i-1] * size[i-1];

	FLA_Obj_create_tensor(FLA_DOUBLE, order, size, stride, A);
}

void test_permute_tensor(dim_t order, dim_t nA, dim_t permutation[order]){
  dim_t i;
  FLA_Obj A, P;

  initTensor(order, nA, &A);
  FLA_Random_tensor(A);

  initTensor(order, nA, &P);

  printf("A = tensor([");
  FLA_Obj_print_tensor(A);
  printf("],[");
  for(i = 0; i < FLA_Obj_order(A); i++)
	printf("%d ", FLA_Obj_dimsize(A,i));
  printf("]);\n\n");

  FLA_Permute(A, permutation, P);

  printf("postA = tensor([");
  FLA_Obj_print_tensor(A);
  printf("],[");
  for(i = 0; i < FLA_Obj_order(A); i++)
	printf("%d ", FLA_Obj_dimsize(A,i));
  printf("]);\n\n");

  printf("P = tensor([");
  FLA_Obj_print_tensor(P);
  printf("],[");
  for(i = 0; i < FLA_Obj_order(P); i++)
	printf("%d ", FLA_Obj_dimsize(P,i));
  printf("]);\n\n");

}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 4){
		Usage();
		FLA_Finalize();
		return 0;
	}

	dim_t i;		
	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	const int nA = atoi(argv[++argNum]);
	dim_t permutation[m];

    for(i = 0; i < m; i++)
		permutation[i] = atoi(argv[++argNum]);

	if(m <= 0 || nA <= 0 ){
		printf("m, nA must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	test_permute_tensor(m, nA, permutation);

	FLA_Finalize();

	return 0;
}