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

  dim_t scalar_sizeA[FLA_MAX_ORDER];
  dim_t block_sizeA[FLA_MAX_ORDER];
  dim_t blocked_sizeA[FLA_MAX_ORDER];
  dim_t stride_A[FLA_MAX_ORDER];

  FLA_Obj* buf;

  for(i = 0; i < order; i++){
	scalar_sizeA[i] = nA;
	block_sizeA[i] = bA;
  }

  FLA_array_elemwise_quotient(order, scalar_sizeA, block_sizeA, blocked_sizeA);
  FLA_Set_tensor_stride(order, blocked_sizeA, stride_A);

  printf("scalar tensor test\n");
  FLA_Obj_create_tensor(FLA_DOUBLE, order, scalar_sizeA, stride_A, &Ascalar);
  FLA_Random_tensor(Ascalar);
  FLA_Obj_print_matlab("Ascalar", Ascalar);
	
  printf("\n\n");
  printf("blocked tensor test\n");
  FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, scalar_sizeA, stride_A, block_sizeA, &Ablocked);
  FLA_Random_tensor(Ablocked);
  FLA_Obj_print_matlab("Ablocked", Ablocked);
  
  printf("\n\n");

  buf = FLA_Obj_base_buffer(Ablocked);

  (*buf).permutation[0] = order - 1;
  for(i = 1; i < order; i++)
  	(*buf).permutation[i] = i-1;

  print_array("permutation", order, buf->permutation);
  printf("permuted blocked tensor test\n");  

  FLA_Obj_print_matlab("Ablocked", Ablocked);

  FLA_Obj_free_buffer(&Ascalar);
  FLA_Obj_free_without_buffer(&Ascalar);

  FLA_Obj_blocked_tensor_free_buffer(&Ablocked);
  FLA_Obj_free_without_buffer(&Ablocked);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* bA){
    int argNum = 0;

    if(argc != 4){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t bA){

	if(order <= 0 || nA <= 0 || bA <= 0){
		printf("m, nA, and bA must be greater than 0\n");
		return FLA_FAILURE;
	}

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		return FLA_FAILURE;
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t bA;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &bA) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, bA) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_tlash_print_routines(order, nA, bA);

	FLA_Finalize();

	return 0;
}
