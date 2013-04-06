#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test psym tensor operations.\n\n");
    printf("  tensor_psym <m> <nA> <bA> <nGroups> <group1Len> ... <groupNLen> <symMode1> ... <symModeM>\n\n");
    printf("  m: order of symmetric tensor\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
	printf("  nGroups: number of symmetric groups\n");
	printf("  groupKLen: length of symmetric group K\n");
	printf("  symModeK: Kth mode in symmetric group ordering\n");
}

void print_tensor(const char* varName, FLA_Obj A){
	dim_t i;
	printf("%s tensor\n", varName);
	printf("%s = tensor([", varName);
	FLA_Obj_print_tensor(A);
	printf("],[");
	for(i = 0; i < FLA_Obj_order(A); i++)
		printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(A)))[0],i) * FLA_Obj_dimsize(A,i));
	printf("]);\n\n");
}

void create_psym_tensor(dim_t order, dim_t nA, dim_t bA, TLA_sym sym, FLA_Obj* A){
	dim_t i;
	dim_t flat_size[order];
	for(i = 0; i < order; i++)
		flat_size[i] = nA;

	dim_t blocked_stride[order];
	blocked_stride[0] = 1;

	for(i = 1; i < order; i++)
		blocked_stride[i] = flat_size[i-1]/bA * blocked_stride[i-1];

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, bA, sym, A);
	FLA_Random_psym_tensor(*A);
}

int main(int argc, char* argv[]){
	dim_t i;

	FLA_Init();

	if(argc < 6){
		Usage();
		FLA_Finalize();
		return 0;
	}
		
	TLA_sym sym;
	int argNum = 0;
	const dim_t m = atoi(argv[++argNum]);
	sym.order = m;
	const dim_t nA = atoi(argv[++argNum]);
	const dim_t bA = atoi(argv[++argNum]);
	sym.nSymGroups = atoi(argv[++argNum]);

	if(argc < 5 + sym.nSymGroups + m){
        Usage();
        FLA_Finalize();
        return 0;
    }

	for(i = 0; i < sym.nSymGroups; i++)
		sym.symGroupLens[i] = atoi(argv[++argNum]);
	for(i = 0; i < m; i++)
		sym.symModes[i] = atoi(argv[++argNum]);


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

	FLA_Obj T;

	create_psym_tensor(m, nA, bA, sym, &T);
	print_tensor("T", T);
	
	FLA_Obj_blocked_psym_tensor_free_buffer(&T);
	FLA_Obj_free_without_buffer(&T);

	FLA_Finalize();

	return 0;
}
