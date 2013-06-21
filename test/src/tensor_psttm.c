#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test psttm operation.\n\n");
    printf("  tensor_psttm <m> <nSym> <symGroupLen1> ... <symGroupK> <symMode1> ... <symModeK> <nA1> ...<nAk> <bA1> ... <bAK> <mode> <nC> <bC> \n\n");
    printf("  m: order of tensor\n");
    printf("  nSym: number symmetric groups\n");
    printf("  symGroupLenK: length of kth symmetric group\n");
    printf("  symModeK: kth symmetric mode\n");
    printf("  nAk: dimension of kth symmetric group of A\n");
    printf("  bAk: block dimension of kth symmetric group of A\n");
    printf("  mode: mode to multiply in\n");
    printf("  nC: dimension of C's <mode>\n");
    printf("  bC: block-size of C's <mode> dimension\n");
}

void create_random_psym_tensor(dim_t order, TLA_sym sym, dim_t n[], dim_t b[], FLA_Obj* A){
    dim_t blocked_size[FLA_MAX_ORDER];
    dim_t blocked_stride[FLA_MAX_ORDER];

    FLA_array_elemwise_quotient(order, n, b, blocked_size);
    FLA_Set_tensor_stride(order, blocked_size, blocked_stride);

    FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, n, blocked_stride, b, sym, A);
    FLA_Random_psym_tensor(*A);
}

void create_random_matrix(dim_t size[2], dim_t blk_size[2], FLA_Obj* obj){
    dim_t order = 2;
    dim_t blked_size[2];
    dim_t blked_stride[2];

    FLA_array_elemwise_quotient(order, size, blk_size, blked_size);
    FLA_Set_tensor_stride(order, blked_size, blked_stride);

    FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size, blked_stride, blk_size, obj);
    FLA_Random_tensor(*obj);
}

void test_psttm(dim_t order, TLA_sym sym, dim_t nA[], dim_t bA[], dim_t mode, dim_t nCMode, dim_t bCMode){
    FLA_Obj A, B, C;

    dim_t i,j;
    dim_t flat_size_A[FLA_MAX_ORDER];
    dim_t blk_size_A[FLA_MAX_ORDER];
    dim_t flat_size_C[FLA_MAX_ORDER];
    dim_t blk_size_C[FLA_MAX_ORDER];
    TLA_sym symC;

    dim_t split_mode_arr[1];

    dim_t bB[2];
    dim_t nB[2];
    for(i = 0; i < sym.nSymGroups; i++){
        dim_t groupOffset = TLA_sym_group_mode_offset(sym, i);
        for(j = 0; j < sym.symGroupLens[i]; j++){
            flat_size_A[sym.symModes[groupOffset + j]] = nA[i];
            blk_size_A[sym.symModes[groupOffset + j]] = bA[i];
        }
    }
    create_random_psym_tensor(order, sym, flat_size_A, blk_size_A, &A);

    symC = sym;

    split_mode_arr[0] = mode;
    TLA_split_sym_group(sym, 1, split_mode_arr, &symC);


    memcpy(&(flat_size_C[0]), &(flat_size_A[0]), order * sizeof(dim_t));
    flat_size_C[mode] = nCMode;
    memcpy(&(blk_size_C[0]), &(blk_size_A[0]), order * sizeof(dim_t));
    blk_size_C[mode] = bCMode;

    create_random_psym_tensor(order, symC, flat_size_C, blk_size_C, &C);

    bB[0] = FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(C))[0], mode);
    bB[1] = FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(A))[0], mode);
    nB[0] = FLA_Obj_dimsize(C, mode) * bB[0];
    nB[1] = FLA_Obj_dimsize(A, mode) * bB[1];

    create_random_matrix(nB, bB, &B);

    FLA_Obj_print_matlab("A", A);
    FLA_Obj_print_matlab("B", B);
    FLA_Obj_print_matlab("preC", C);

    FLA_Psttm(FLA_ONE, A, mode, FLA_ONE, B, C);
    FLA_Obj_print_matlab("C", C);

    FLA_Obj_blocked_psym_tensor_free_buffer(&A);
    FLA_Obj_blocked_tensor_free_buffer(&B);
    FLA_Obj_blocked_psym_tensor_free_buffer(&C);

    FLA_Obj_free_without_buffer(&A);
    FLA_Obj_free_without_buffer(&B);
    FLA_Obj_free_without_buffer(&C);
}

FLA_Error check_errors(dim_t order, TLA_sym sym, dim_t nA[], dim_t bA[], dim_t mode, dim_t nC, dim_t bC){
    dim_t i, j;
    dim_t count = 0;
    dim_t symModesFound[FLA_MAX_ORDER];
    if(order <= 0){
        printf("m must be greater than 0\n");
        return FLA_FAILURE;
    }

    if(sym.nSymGroups <= 0){
        printf("nSymGroups must be greater than 0\n");
        return FLA_FAILURE;
    }

    for(i = 0; i < sym.nSymGroups; i++)
        count += sym.symGroupLens[i];
    if(count != order){
        printf("symGroupLens must sum to m\n");
        return FLA_FAILURE;
    }

    memset(&(symModesFound[0]), 0, order * sizeof(dim_t));
    for(i = 0; i < order; i++){
        for(j = 0; j < sym.order; j++)
            if(sym.symModes[j] == i)
                symModesFound[i] = 1;
    }
    for(i = 0; i < order; i++)
        if(symModesFound[i] == 0){
            printf("symModes missing mode %d\n", i);
            return FLA_FAILURE;
        }

    for(i = 0; i < sym.nSymGroups; i++){
        if(nA[i] <= 0){
            printf("n[%d] must be greater than 0\n", i);
            return FLA_FAILURE;
        }
        if(bA[i] <= 0){
            printf("b[%d] must be greater than 0\n", i);
            return FLA_FAILURE;
        }
        if(nA[i] % bA[i] != 0){
            printf("b[%d] must evenly divide n[%d]\n", i, i);
            return FLA_FAILURE;
        }
    }

    if(nC <= 0){
        printf("nC must be greater than 0\n");
        return FLA_FAILURE;
    }
    if(bC <= 0){
        printf("bC must be greater than 0\n");
        return FLA_FAILURE;
    }
    if(nC % bC != 0){
        printf("bC must evenly divide nC\n");
        return FLA_FAILURE;
    }

    return FLA_SUCCESS;
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, TLA_sym* sym, dim_t* nA, dim_t* bA, dim_t* mode, dim_t* nC, dim_t* bC){
    dim_t i;
    int argNum = 0;

    if(argc < 3){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    sym->order = *order;
    sym->nSymGroups = atoi(argv[++argNum]);

    if(argc != 3 + sym->nSymGroups + *order + sym->nSymGroups + sym->nSymGroups + 3){
        return FLA_FAILURE;
    }

    for(i = 0; i < sym->nSymGroups; i++)
        (sym->symGroupLens)[i] = atoi(argv[++argNum]);
    for(i = 0; i < *order; i++)
        (sym->symModes)[i] = atoi(argv[++argNum]);


    for(i = 0; i < sym->nSymGroups; i++)
        nA[i] = atoi(argv[++argNum]);

    for(i = 0; i < sym->nSymGroups; i++)
        bA[i] = atoi(argv[++argNum]);

    *mode = atoi(argv[++argNum]);
    *nC = atoi(argv[++argNum]);
    *bC = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	TLA_sym sym;
	dim_t nA[FLA_MAX_ORDER];
	dim_t bA[FLA_MAX_ORDER];
	dim_t mode;
	dim_t nC;
	dim_t bC;

	FLA_Init();


    //Parse inputs
    if(parse_input(argc, argv, &order, &sym, nA, bA, &mode, &nC, &bC) == FLA_FAILURE){
        Usage();
        FLA_Finalize();
        return 0;
    }

    //Error check
    if(check_errors(order, sym, nA, bA, mode, nC, bC) == FLA_FAILURE){
        FLA_Finalize();
        return 0;
    }

	test_psttm(order, sym, nA, bA, mode, nC, bC);

	FLA_Finalize();

	return 0;
}
