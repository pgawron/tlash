#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test TLASH Print operations.\n\n");
    printf("  tlash_part <m> <nA> <bA>\n\n");
    printf("  m: order of tensor\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
}

void test_tlash_print_routines(dim_t order, dim_t nA, dim_t bA){
  dim_t i;
  FLA_Obj Ascalar;
  FLA_Obj Ablocked;

  dim_t scalar_sizeA[order];
  dim_t block_sizeA[order];
  dim_t stride_A[order];

  for(i = 0; i < order; i++){
	scalar_sizeA[i] = nA;
	block_sizeA[i] = bA;
  }

  stride_A[0] = 1;
  for(i = 1; i < order; i++)
	stride_A[i] = scalar_sizeA[i-1]*stride_A[i-1];

  printf("scalar tensor test\n");
  FLA_Obj_create_tensor(FLA_DOUBLE, order, scalar_sizeA, stride_A, &Ascalar);
  FLA_Random_tensor(Ascalar);
  FLA_Obj_print_tensor2(Ascalar);
	
  printf("\n\n");
  printf("blocked tensor test\n");
  FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, scalar_sizeA, stride_A, block_sizeA, &Ablocked);
  FLA_Random_tensor(Ablocked);
	
  FLA_Obj_print_tensor2(Ablocked);

  printf("\n\ncheck\n");
  FLA_Obj* buffer = FLA_Obj_base_buffer(Ablocked);
  for(i = 0; i < FLA_Obj_num_elem_alloc(Ablocked); i++){
  	FLA_Obj_print_tensor2(buffer[i]);
  	printf("\n");
  	}

  printf("\n\n");
  printf("blocked symtensor test\n");

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

	if(m <= 0 || nA <= 0 || bA <= 0){
		printf("m, nA, and bA must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	test_tlash_print_routines(m, nA, bA);

	FLA_Finalize();

	return 0;
}
