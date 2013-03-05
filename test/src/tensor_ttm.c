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

void test_permute_ttm_scalar(dim_t order, dim_t sizeA[order], dim_t nC, dim_t mode){
	printf("\nPERMUTE TEST\n\n");
	dim_t i;
	dim_t sizeB[2] = {nC, sizeA[mode]};
	dim_t sizeC[order];
	memcpy(&(sizeC[0]), &(sizeA[0]), order * sizeof(dim_t));

	sizeC[mode] = nC;

	dim_t permutation[order];
	for(i = 1; i < order; i++)
		permutation[i] = i - 1;
	permutation[0] = order - 1;

	dim_t stored_sizeA[order];
	for(i = 0; i < order-1; i++)
		stored_sizeA[i] = sizeA[i+1];
	stored_sizeA[order - 1] = sizeA[0];

	dim_t strideA[order];
	dim_t strideB[2] = {1, nC};
	dim_t strideC[order];

	strideA[0] = 1;
	strideC[0] = 1;
	for(i = 1; i < order; i++){
		strideA[i] = stored_sizeA[i-1] * strideA[i-1];
		strideC[i] = sizeC[i-1] * strideC[i-1];
	}

	FLA_Obj A,B,C;
	FLA_Obj_create_tensor(FLA_DOUBLE, order, stored_sizeA, strideA, &A);
	FLA_Obj_create_tensor(FLA_DOUBLE, 2, sizeB, strideB, &B);
	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeC, strideC, &C);

	FLA_Random_tensor(A);
	FLA_Random_tensor(B);
	FLA_Random_tensor(C);

	memcpy(&(A.permutation[0]), &(permutation[0]), order * sizeof(dim_t));

	printf("permutation[");
	for(i = 0; i < order; i++)
		printf("%d ", A.permutation[i]);
	printf("]\n");

	printf("\nA = permute(tensor([");
	FLA_Obj_print_tensor(A);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(A); i++)
		printf("%d ", FLA_Obj_dimsize(A,i));
	printf("]),[");
	for(i = 0; i < order; i++)
		printf("%d ", A.permutation[i]+1);
	printf("]);\n");
	

	printf("\nB = reshape([");
	double* buf_B = FLA_Obj_base_buffer(B);
	for(i = 0; i < nC * sizeA[mode]; i++)
		printf("%.3f ", buf_B[i]); 
	printf("],[%d %d]);", nC, sizeA[mode]);

	printf("\nC = tensor([");
	FLA_Obj_print_tensor(C);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(C); i++)
		printf("%d ", FLA_Obj_dimsize(C,i));
	printf("]);");

	FLA_Ttm_single_mode(FLA_ONE, A, mode, FLA_ONE, B, C);

	printf("\npostC = tensor([");
	FLA_Obj_print_tensor(C);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(C); i++)
		printf("%d ", FLA_Obj_dimsize(C,i));
	printf("]);");
}
void test_ttm_scalar(dim_t order, dim_t sizeA[order], dim_t nC, dim_t mode){
	dim_t sizeB[2] = {nC, sizeA[mode]};
	dim_t sizeC[order];
	memcpy(&(sizeC[0]), &(sizeA[0]), order * sizeof(dim_t));

	sizeC[mode] = nC;

	dim_t strideA[order];
	dim_t strideB[2] = {1, nC};
	dim_t strideC[order];

	strideA[0] = 1;
	strideC[0] = 1;
	dim_t i;
	for(i = 1; i < order; i++){
		strideA[i] = sizeA[i-1] * strideA[i-1];
		strideC[i] = sizeC[i-1] * strideC[i-1];
	}

	FLA_Obj A,B,C;
	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeA, strideA, &A);
	FLA_Obj_create_tensor(FLA_DOUBLE, 2, sizeB, strideB, &B);
	FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeC, strideC, &C);

	FLA_Random_tensor(A);
	FLA_Random_tensor(B);
	FLA_Random_tensor(C);

	printf("\nA = tensor([");
	FLA_Obj_print_tensor(A);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(A); i++)
		printf("%d ", FLA_Obj_dimsize(A,i));
	printf("]);");

	printf("\nB = reshape([");
	double* buf_B = FLA_Obj_base_buffer(B);
	for(i = 0; i < nC * sizeA[mode]; i++)
		printf("%.3f ", buf_B[i]); 
	printf("],[%d %d]);", nC, sizeA[mode]);

	printf("\nC = tensor([");
	FLA_Obj_print_tensor(C);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(C); i++)
		printf("%d ", FLA_Obj_dimsize(C,i));
	printf("]);");


	FLA_Ttm_single_mode(FLA_ONE, A, mode, FLA_ONE, B, C);
	printf("\npostC = tensor([");
	FLA_Obj_print_tensor(C);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(C); i++)
		printf("%d ", FLA_Obj_dimsize(C,i));
	printf("]);");
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
	dim_t nA[m];
	for(i = 0; i < m; i++)
		nA[i] = atoi(argv[++argNum]);
	const int nC = atoi(argv[++argNum]);
	const int mode = atoi(argv[++argNum]);

/*
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
*/

	test_ttm_scalar(m, nA, nC, mode);
	test_permute_ttm_scalar(m, nA, nC, mode);

	FLA_Finalize();

	return 0;
}
