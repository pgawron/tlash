#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test ttm operation.\n\n");
    printf("  tensor_ttm <m> <nA1> ... <nAm> <nCk> <mode>\n\n");
    printf("  m: order of tensor\n");
    printf("  nAk: dimension of kth mode of A\n");
    printf("  nCk: dimension of kth mode of C\n");
    printf("  mode: mode to multiply in\n");
}

void test_permute_ttm_scalar(dim_t order, dim_t sizeA[], dim_t nC, dim_t mode){
	dim_t i;
	dim_t sizeB[2] = {nC, sizeA[mode]};
	dim_t sizeC[FLA_MAX_ORDER];
	dim_t permutation[FLA_MAX_ORDER];
	dim_t stored_sizeA[FLA_MAX_ORDER];

	dim_t strideA[FLA_MAX_ORDER];
	dim_t strideB[2];
	dim_t strideC[FLA_MAX_ORDER];

	FLA_Obj A,B,C;

	printf("\nPERMUTE TEST\n\n");
	memcpy(&(sizeC[0]), &(sizeA[0]), order * sizeof(dim_t));

	sizeC[mode] = nC;


	for(i = 1; i < order; i++)
		permutation[i] = i - 1;
	permutation[0] = order - 1;


	for(i = 0; i < order-1; i++)
		stored_sizeA[i] = sizeA[i+1];
	stored_sizeA[order - 1] = sizeA[0];

	FLA_Set_tensor_stride(order, stored_sizeA, strideA);
	FLA_Set_tensor_stride(2, sizeB, strideB);
	FLA_Set_tensor_stride(order, sizeC, strideC);


	FLA_Obj_create_tensor(FLA_DOUBLE, order, stored_sizeA, strideA, &A);
	FLA_Obj_create_tensor(FLA_DOUBLE, 2, sizeB, strideB, &B);
	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeC, strideC, &C);

	FLA_Random_tensor(A);
	FLA_Random_tensor(B);
	FLA_Random_tensor(C);

	memcpy(&(A.permutation[0]), &(permutation[0]), order * sizeof(dim_t));

	print_array("permutation", order, A.permutation);

	FLA_Obj_print_matlab("A", A);
	FLA_Obj_print_matlab("B", B);
	FLA_Obj_print_matlab("preC", C);

	FLA_Ttm_single_mode(FLA_ONE, A, mode, FLA_ONE, B, C);

	FLA_Obj_print_matlab("postC", C);

	FLA_Obj_free_buffer(&A);
	FLA_Obj_free_buffer(&B);
	FLA_Obj_free_buffer(&C);

	FLA_Obj_free_without_buffer(&A);
	FLA_Obj_free_without_buffer(&B);
	FLA_Obj_free_without_buffer(&C);
}

void test_ttm_scalar(dim_t order, dim_t sizeA[], dim_t nC, dim_t mode){
	dim_t sizeB[2] = {nC, sizeA[mode]};
	dim_t sizeC[FLA_MAX_ORDER];

	dim_t strideA[FLA_MAX_ORDER];
	dim_t strideB[2];
	dim_t strideC[FLA_MAX_ORDER];

	FLA_Obj A,B,C;

	memcpy(&(sizeC[0]), &(sizeA[0]), order * sizeof(dim_t));

	sizeC[mode] = nC;

	FLA_Set_tensor_stride(order, sizeA, strideA);
	FLA_Set_tensor_stride(2, sizeB, strideB);
	FLA_Set_tensor_stride(order, sizeC, strideC);

	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeA, strideA, &A);
	FLA_Obj_create_tensor(FLA_DOUBLE, 2, sizeB, strideB, &B);
	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeC, strideC, &C);

	FLA_Random_tensor(A);
	FLA_Random_tensor(B);
	FLA_Random_tensor(C);

	FLA_Obj_print_matlab("A", A);
	FLA_Obj_print_matlab("B", B);
	FLA_Obj_print_matlab("preC", C);

	FLA_Ttm_single_mode(FLA_ONE, A, mode, FLA_ONE, B, C);
	FLA_Obj_print_matlab("C", C);

	FLA_Obj_free_buffer(&A);
	FLA_Obj_free_buffer(&B);
	FLA_Obj_free_buffer(&C);

	FLA_Obj_free_without_buffer(&A);
	FLA_Obj_free_without_buffer(&B);
	FLA_Obj_free_without_buffer(&C);
}


FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* nC, dim_t* mode){
    dim_t i;
    int argNum = 0;

    if(argc < 2){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);

    if(argc < 2 + (*order)){
    	return FLA_FAILURE;
    }

    for(i = 0; i < *order; i++){
    	nA[i] = atoi(argv[++argNum]);
    }

    if(argc != 2 + (*order) + 2){
    	return FLA_FAILURE;
    }

    *nC = atoi(argv[++argNum]);
    *mode = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA[], dim_t nC, dim_t mode){
	dim_t i;

	if(order <= 0 || nC <= 0){
		printf("m and nC must be greater than 0\n");
		return FLA_FAILURE;
	}

	if(mode >= order){
		printf("mode must be less than order\n");
		return FLA_FAILURE;
	}

	for(i = 0; i < order; i++){
		if(nA[i] <= 0){
			printf("nA[i] must be greater than 0\n");
			return FLA_FAILURE;
		}
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA[FLA_MAX_ORDER];
	dim_t nC;
	dim_t mode;

	FLA_Init();

	if(parse_input(argc, argv, &order, nA, &nC, &mode) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, nC, mode) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	//test_ttm_scalar(order, nA, nC, mode);
	test_permute_ttm_scalar(order, nA, nC, mode);

	FLA_Finalize();

	return 0;
}
