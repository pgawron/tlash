#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sym tensor operations.\n\n");
    printf("  tensor_sym <m> <nA> <bA>\n\n");
    printf("  m: order of symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
}

void test_create_sym_tensor(int order, int nA, int bA){
	dim_t i;
	dim_t size_A[FLA_MAX_ORDER];
	dim_t blk_size[FLA_MAX_ORDER];
	dim_t blked_size[FLA_MAX_ORDER];
	dim_t stride_A[FLA_MAX_ORDER];
	FLA_Obj A;
	TLA_sym sym;

	for (i = 0; i < order; i++) {
		size_A[i] = nA;
		blk_size[i] = bA;
		blked_size[i] = nA / bA;
	}

	FLA_Set_tensor_stride(order, blked_size, stride_A);

	sym.order = order;
	sym.nSymGroups = 1;
	sym.symGroupLens[0] = order;
	for (i = 0; i < order; i++)
		sym.symModes[i] = i;

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size_A, stride_A,
			blk_size, sym, &A);
	FLA_Random_psym_tensor(A);

	FLA_Obj_print_matlab("A", A);

	FLA_Obj_blocked_psym_tensor_free_buffer(&A);
	FLA_Obj_free_without_buffer(&A);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* bA){
    int argNum = 0;

    if(argc < 3){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t bA){

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		return FLA_FAILURE;
	}
	if(order <= 0 || nA <= 0 || bA <= 0){
		printf("m, nA and bA must be greater than 0\n");
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
		
	test_create_sym_tensor(order, nA, bA);

	FLA_Finalize();

	return 0;
}
