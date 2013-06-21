#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test TLASH Permute operation.\n\n");
    printf("  tensor_permute <m> <nA> <perm0> <perm1> ...\n\n");
    printf("  m: order of tensor\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  permK: permutation index K\n");
}

void initTensor(dim_t order, dim_t nA, FLA_Obj* A){
	dim_t i;	
	dim_t size[FLA_MAX_ORDER];
	dim_t stride[FLA_MAX_ORDER];
	stride[0] = 1;

	for(i = 0; i < order; i++)
		size[i] = nA;

	for(i = 1; i < order; i++)
		stride[i] = stride[i-1] * size[i-1];

	FLA_Obj_create_tensor(FLA_DOUBLE, order, size, stride, A);
}

void test_permute_tensor(dim_t order, dim_t nA, dim_t permutation[]){
  FLA_Obj A, P;

  initTensor(order, nA, &A);
  FLA_Random_tensor(A);

  initTensor(order, nA, &P);

  FLA_Obj_print_matlab("A", A);

  FLA_Permute(A, permutation, &P);

  FLA_Obj_print_matlab("P", P);

}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* permutation){
    dim_t i;
    int argNum = 0;

    if(argc < 3){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);

    if(argc != 3 + (*order)){
    	return FLA_FAILURE;
    }

    for(i = 0; i < *order; i++)
		permutation[i] = atoi(argv[++argNum]);


    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t permutation[]){
	dim_t i, j;
	if(order <= 0 || nA <= 0 ){
			printf("m, nA must be greater than 0\n");
			return FLA_FAILURE;
	}

	for(i = 0; i < order; i++){
		if(permutation[i] >= order){
			printf("permutation[i] must be less than order\n");
			return FLA_FAILURE;
		}
		for(j = 0; j < order; j++){
			if((i != j) && (permutation[i] == permutation[j])){
				printf("permutation must not contain duplicates\n");
				return FLA_FAILURE;
			}
		}
	}
	return FLA_SUCCESS;
}


int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t permutation[FLA_MAX_ORDER];

	FLA_Init();

	//Parse input
	if(parse_input(argc, argv, &order, &nA, permutation) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, permutation) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_permute_tensor(order, nA, permutation);

	FLA_Finalize();

	return 0;
}
