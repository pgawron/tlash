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

void create_random_psym_tensor(dim_t order, TLA_sym sym, dim_t n[sym.nSymGroups], dim_t b[sym.nSymGroups], FLA_Obj* A){
    dim_t i;
    dim_t blocked_size[order];

    for(i = 0; i < order; i++)
        blocked_size[i] = n[i] / b[i];

    dim_t blocked_stride[order];
    blocked_stride[0] = 1;
    for(i = 1; i < order; i++)
        blocked_stride[i] = blocked_size[i-1] * blocked_stride[i-1];

    FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, n, blocked_stride, b, sym, A);
    FLA_Random_psym_tensor(*A);
}

void create_random_matrix(dim_t size[2], dim_t blk_size[2], FLA_Obj* obj){
    dim_t i,j;
    dim_t order = 2;
    dim_t blked_size[] = {size[0] / blk_size[0], size[1] / blk_size[1]};
    dim_t blked_stride[] = {1, blked_size[0]};
    dim_t blk_stride[] = {1, blk_size[0]};

    FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, blked_size, obj);
    obj->base->elemtype = FLA_TENSOR;

    FLA_Obj* buf = (FLA_Obj*)FLA_malloc(blked_size[0] * blked_size[1] * sizeof(FLA_Obj));
    FLA_Obj_attach_buffer_to_tensor(buf, order, blked_stride, obj);

    FLA_Adjust_2D_info(obj);

    for(i = 0; i < blked_size[1]; i++)
      for(j = 0; j < blked_size[0]; j++){
          FLA_Obj* curObj = FLA_Obj_base_buffer(*obj);
          FLA_Obj_create_tensor(FLA_DOUBLE, order, blk_size, blk_stride, &(curObj[j + (blked_size[0]*i)]));
          FLA_Random_matrix(curObj[j+ (blked_size[0]*i)]);
          FLA_Adjust_2D_info(&(curObj[j+(blked_size[0]*i)]));
      }
}

void test_psttm(dim_t order, TLA_sym sym, dim_t nA[sym.nSymGroups], dim_t bA[sym.nSymGroups], dim_t mode, dim_t nCMode, dim_t bCMode){
    FLA_Obj A, B, C;

    dim_t i,j;
    dim_t flat_size_A[order];
    dim_t blk_size_A[order];

    for(i = 0; i < sym.nSymGroups; i++){
        dim_t groupOffset = TLA_sym_group_mode_offset(sym, i);
        for(j = 0; j < sym.symGroupLens[i]; j++){
            flat_size_A[sym.symModes[groupOffset + j]] = nA[i];
            blk_size_A[sym.symModes[groupOffset + j]] = bA[i];
        }
    }
    create_random_psym_tensor(order, sym, flat_size_A, blk_size_A, &A);

    TLA_sym symC;
    symC = sym;

    dim_t split_mode_arr[1];
    split_mode_arr[0] = mode;
    TLA_split_sym_group(sym, 1, split_mode_arr, &symC);

    dim_t flat_size_C[order];
    dim_t blk_size_C[order];

    memcpy(&(flat_size_C[0]), &(flat_size_A[0]), order * sizeof(dim_t));
    flat_size_C[mode] = nCMode;
    memcpy(&(blk_size_C[0]), &(blk_size_A[0]), order * sizeof(dim_t));
    blk_size_C[mode] = bCMode;

    create_random_psym_tensor(order, symC, flat_size_C, blk_size_C, &C);

    dim_t bB[] = {FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(C))[0], mode), FLA_Obj_dimsize(((FLA_Obj*)FLA_Obj_base_buffer(A))[0], mode)};
    dim_t nB[] = {FLA_Obj_dimsize(C, mode) * bB[0], FLA_Obj_dimsize(A, mode) * bB[1]};
    create_random_matrix(nB, bB, &B);

    FLA_Obj_print_matlab("A", A);
    FLA_Obj_print_matlab("B", B);
    FLA_Obj_print_matlab("preC", C);

    FLA_Psttm(FLA_ONE, A, mode, FLA_ONE, B, C);
    FLA_Obj_print_matlab("C", C);

}

FLA_Error check_errors(dim_t order, TLA_sym sym, dim_t nA[], dim_t bA[], dim_t mode, dim_t nC, dim_t bC){
    dim_t i, j;
    if(order <= 0){
        printf("m must be greater than 0\n");
        return FLA_FAILURE;
    }

    if(sym.nSymGroups <= 0){
        printf("nSymGroups must be greater than 0\n");
        return FLA_FAILURE;
    }
    dim_t count = 0;
    for(i = 0; i < sym.nSymGroups; i++)
        count += sym.symGroupLens[i];
    if(count != order){
        printf("symGroupLens must sum to m\n");
        return FLA_FAILURE;
    }
    dim_t symModesFound[order];
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
            FLA_Finalize();
        }
        if(bA[i] <= 0){
            printf("b[%d] must be greater than 0\n", i);
            FLA_Finalize();
        }
        if(nA[i] % bA[i] != 0){
            printf("b[%d] must evenly divide n[%d]\n", i, i);
            FLA_Finalize();
            return 0;
        }
    }

    if(nC <= 0){
        printf("nC must be greater than 0\n");
        FLA_Finalize();
    }
    if(bC <= 0){
        printf("bC must be greater than 0\n");
        FLA_Finalize();
    }
    if(nC % bC != 0){
        printf("bC must evenly divide nC\n");
        FLA_Finalize();
        return 0;
    }

    return FLA_SUCCESS;
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, TLA_sym* sym, dim_t** nA, dim_t** bA, dim_t* mode, dim_t* nC, dim_t* bC){
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


    *nA = FLA_malloc(sym->nSymGroups * sizeof(dim_t));
    for(i = 0; i < sym->nSymGroups; i++)
        (*nA)[i] = atoi(argv[++argNum]);

    *bA = FLA_malloc(sym->nSymGroups * sizeof(dim_t));
    for(i = 0; i < sym->nSymGroups; i++)
            (*bA)[i] = atoi(argv[++argNum]);

    *mode = atoi(argv[++argNum]);
    *nC = atoi(argv[++argNum]);
    *bC = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	FLA_Init();

    dim_t order;
    TLA_sym sym;
    dim_t* nA;
    dim_t* bA;
    dim_t mode;
    dim_t nC;
    dim_t bC;

    //Parse inputs
    if(parse_input(argc, argv, &order, &sym, &nA, &bA, &mode, &nC, &bC) == FLA_FAILURE){
        Usage();
        FLA_Finalize();
        return 0;
    }

    //Error check
    if(check_errors(order, sym, nA, bA, mode, nC, bC) == FLA_FAILURE){
        FLA_free(nA);
        FLA_free(bA);
        FLA_Finalize();
        return 0;
    }

	test_psttm(order, sym, nA, bA, mode, nC, bC);

    FLA_free(nA);
    FLA_free(bA);
	FLA_Finalize();

	return 0;
}
