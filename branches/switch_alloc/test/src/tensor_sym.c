#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sym tensor operations.\n\n");
    printf("  tensor_sym <m> <nA> <bA>\n\n");
    printf("  m: order of hyper-symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
}

void test_create_sym_tensor(int order, int nA, int bA){
	dim_t i;
	dim_t size_A[order];
	dim_t blk_size[order];
	for(i = 0; i < order; i++){
		size_A[i] = nA;
		blk_size[i] = bA;
	}

	dim_t stride_A[order];
	stride_A[0] = 1;

	for(i = 1; i < order; i++)
		stride_A[i] = size_A[i-1]/bA * stride_A[i-1];

	FLA_Obj A;

	TLA_sym sym;
	TLA_Sym_init_nonsymmetric(order, &sym);
	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size_A, stride_A, blk_size, sym, &A);
	FLA_Random_psym_tensor(A);
  
	printf("A tensor\n");
	printf("a = tensor([");
	FLA_Obj_print_tensor(A);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(A); i++)
		printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(A)))[0],i) * FLA_Obj_dimsize(A,i));
	printf("]);\n\n");

  FLA_Obj_blocked_psym_tensor_free_buffer(&A);
  FLA_Obj_free_without_buffer(&A);
}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 4){
		Usage();
		FLA_Finalize();
		return 0;
	}
		
	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	const int nA = atoi(argv[++argNum]);
	const int bA = atoi(argv[++argNum]);

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		FLA_Finalize();
		return 0;
	}
	if(m <= 0 || bA <= 0){
		printf("m and bA must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	test_create_sym_tensor(m, nA, bA);

	FLA_Finalize();

	return 0;
}
