#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test psttm operation.\n\n");
    printf("  tensor_psttm <m> <nSym> <symGroupLen1> ... <symGroupK> <symMode1> ... <symModeK> <nA1> ...<nAk> <mode> <nC> <bA> \n\n");
    printf("  m: order of tensor\n");
    printf("  nSym: number symmetric groups\n");
    printf("  KgroupLen: length of kth symmetric group\n");
    printf("  symModeK: kth symmetric mode\n");
    printf("  nAk: dimension of kth symmetric group of A\n");
    printf("  mode: mode to multiply in\n");
    printf("  nC: dimension of kth mode of C\n");
    printf("  bA: block-size of A\n");
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

void create_random_psym_tensor(dim_t order, dim_t size[order], dim_t b, TLA_sym sym, FLA_Obj* obj){
    dim_t i;
    dim_t blocked_stride[order];
    blocked_stride[0] = 1;
    for(i = 1; i < order; i++)
        blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / b);
    FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size, blocked_stride, b, sym, obj);
    FLA_Random_psym_tensor(*obj);
}

void create_random_matrix(dim_t size[2], dim_t b, FLA_Obj* obj){
    dim_t i,j;
    dim_t order = 2;
    dim_t sizeObj[] = {size[0] / b, size[1] / b};
    dim_t strideObj[] = {1, sizeObj[0]};
    dim_t sizeBlk[] = {b, b};
    dim_t strideBlk[] = {1, b};

    FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, sizeObj, obj);
    obj->base->elemtype = FLA_TENSOR;

    FLA_Obj* buf = (FLA_Obj*)FLA_malloc(sizeObj[0] * sizeObj[1] * sizeof(FLA_Obj));
    FLA_Obj_attach_buffer_to_tensor(buf, order, strideObj, obj);

    FLA_Adjust_2D_info(obj);

    for(i = 0; i < sizeObj[1]; i++)
      for(j = 0; j < sizeObj[0]; j++){
          FLA_Obj* curObj = FLA_Obj_base_buffer(*obj);
          FLA_Obj_create_tensor(FLA_DOUBLE, order, sizeBlk, strideBlk, &(curObj[j + (sizeObj[0]*i)]));
          FLA_Random_matrix(curObj[j+ (sizeObj[0]*i)]);
          FLA_Adjust_2D_info(&(curObj[j+(sizeObj[0]*i)]));
      }
}

void test_psttm(dim_t order, TLA_sym sym, dim_t size_A[order], dim_t bA, dim_t mode, dim_t nC){
    FLA_Obj A, B, C;

    create_random_psym_tensor(order, size_A, bA, sym, &A);

    TLA_sym symC;
    symC = sym;

    dim_t split_mode_arr[1];
    split_mode_arr[0] = mode;
    TLA_split_sym_group(sym, 1, split_mode_arr, &symC);
    dim_t size_C[order];
    memcpy(&(size_C[0]), &(size_A[0]), order * sizeof(dim_t));
    size_C[mode] = nC;

    create_random_psym_tensor(order, size_C, bA, symC, &C);

    dim_t size_B[] = {size_C[mode], size_A[mode]};
    create_random_matrix(size_B, bA, &B);

    print_tensor("A", A);
    print_tensor("B", B);
    print_tensor("C", C);

    FLA_Psttm(FLA_ONE, A, mode, FLA_ONE, B, C);
    print_tensor("C", C);

}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 3){
		Usage();
		FLA_Finalize();
		return 0;
	}
	
	dim_t i,j;
	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	TLA_sym sym;

	sym.nSymGroups = atoi(argv[++argNum]);
	sym.order = m;

	if(argc != 3 + sym.nSymGroups + m + sym.nSymGroups + 3){
	    Usage();
        FLA_Finalize();
        return 0;
	}

    dim_t nA[m];

    for(i = 0; i < sym.nSymGroups; i++)
		sym.symGroupLens[i] = atoi(argv[++argNum]);

	for(i = 0; i < m; i++)
	    sym.symModes[i] = atoi(argv[++argNum]);

	dim_t count = 0;
	for(i = 0; i < sym.nSymGroups; i++){
	    ++argNum;
	    for(j = 0; j < sym.symGroupLens[i]; j++)
	        nA[sym.symModes[count++]] = atoi(argv[argNum]);
	}


	const int mode = atoi(argv[++argNum]);
	const int nC = atoi(argv[++argNum]);
	const int bA = atoi(argv[++argNum]);



	test_psttm(m, sym, nA, bA, mode, nC);

	FLA_Finalize();

	return 0;
}
