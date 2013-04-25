#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test psym tensor operations.\n\n");
    printf("  tensor_psym <m> <nGroups> <group1Len> ... <groupNLen> <symMode1> ... <symModeM> <n1> ... <nK> <b1> ... <bK> <perm1> ... <permK>\n\n");
    printf("  m: order of symmetric tensor\n");
	printf("  nGroups: number of symmetric groups\n");
	printf("  groupKLen: length of symmetric group K\n");
	printf("  symModeK: Kth mode in symmetric group ordering\n");
    printf("  nK: mode-length of symmetric group K\n");
    printf("  bK: mode-length of block of A\n");
    printf("  permK: permutation index K\n");
}

void create_psym_tensor(dim_t order, TLA_sym sym, dim_t n[sym.nSymGroups], dim_t b[sym.nSymGroups], FLA_Obj* A){
	dim_t i, j;
	dim_t flat_size[order];
	dim_t blocked_size[order];
	dim_t blk_size[order];
	for(i = 0; i < sym.nSymGroups; i++){
	    dim_t groupOffset = TLA_sym_group_mode_offset(sym, i);
	    for(j = 0; j < sym.symGroupLens[i]; j++){
	        flat_size[sym.symModes[groupOffset + j]] = n[i];
	        blk_size[sym.symModes[groupOffset + j]] = b[i];
	        blocked_size[sym.symModes[groupOffset + j]] = n[i] / b[i];
	    }
	}

	dim_t blocked_stride[order];
	blocked_stride[0] = 1;
	for(i = 1; i < order; i++)
		blocked_stride[i] = blocked_size[i-1] * blocked_stride[i-1];

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, blk_size, sym, A);
	FLA_Random_psym_tensor(*A);
}

FLA_Error check_errors(dim_t m, TLA_sym sym, dim_t n[], dim_t b[], dim_t permutation[]){
    dim_t i, j;
    if(m <= 0){
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
    if(count != m){
        printf("symGroupLens must sum to m\n");
        return FLA_FAILURE;
    }
    dim_t symModesFound[m];
    memset(&(symModesFound[0]), 0, m * sizeof(dim_t));
    for(i = 0; i < m; i++){
        for(j = 0; j < sym.order; j++)
            if(sym.symModes[j] == i)
                symModesFound[i] = 1;
    }
    for(i = 0; i < m; i++)
        if(symModesFound[i] == 0){
            printf("symModes missing mode %d\n", i);
            return FLA_FAILURE;
        }

    for(i = 0; i < sym.nSymGroups; i++){
        if(n[i] <= 0){
            printf("n[%d] must be greater than 0\n", i);
            FLA_Finalize();
        }
        if(b[i] <= 0){
            printf("b[%d] must be greater than 0\n", i);
            FLA_Finalize();
        }
        if(n[i] % b[i] != 0){
            printf("b[%d] must evenly divide n[%d]\n", i, i);
            FLA_Finalize();
            return 0;
        }
    }

    dim_t permutationModesFound[m];
    memset(&(permutationModesFound[0]), 0, m * sizeof(dim_t));
    for(i = 0; i < m; i++){
        for(j = 0; j < sym.order; j++)
            if(permutation[j] == i)
                permutationModesFound[i] = 1;
    }
    for(i = 0; i < m; i++)
        if(permutationModesFound[i] == 0){
            printf("permutationModes missing mode %d\n", i);
            return FLA_FAILURE;
        }

    return FLA_SUCCESS;
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, TLA_sym* sym, dim_t** n, dim_t** b, dim_t** permutation){
    dim_t i;
    int argNum = 0;
    if(argc < 3){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    sym->order = *order;
    sym->nSymGroups = atoi(argv[++argNum]);

    if(argc != 3 + sym->nSymGroups + *order + sym->nSymGroups + sym->nSymGroups + *order){
        return FLA_FAILURE;
    }

    for(i = 0; i < sym->nSymGroups; i++)
        (sym->symGroupLens)[i] = atoi(argv[++argNum]);
    for(i = 0; i < *order; i++)
        (sym->symModes)[i] = atoi(argv[++argNum]);


    *n = FLA_malloc(sym->nSymGroups * sizeof(dim_t));
    for(i = 0; i < sym->nSymGroups; i++)
        (*n)[i] = atoi(argv[++argNum]);

    *b = FLA_malloc(sym->nSymGroups * sizeof(dim_t));
    for(i = 0; i < sym->nSymGroups; i++)
        (*b)[i] = atoi(argv[++argNum]);

    *permutation = FLA_malloc((*order) * sizeof(dim_t));
    for(i = 0; i < *order; i++)
        (*permutation)[i] = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

void test_permute_tensor(dim_t permutation[], FLA_Obj A){
    printf("Testing permutation routine\n");
    print_array("permutation", A.order, permutation);
    printf("-----------------------------\n");

    dim_t i,j;
    dim_t num_elem_alloc_A = FLA_Obj_num_elem_alloc(A);
    for(i = 0; i < num_elem_alloc_A; i++){
        FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(A))[i];
        printf("A[%d] ", i);
        print_array("stored permutation", curObj.order, curObj.permutation);
        print_array("stored size", curObj.order, curObj.size);
        FLA_Obj_print_matlab("data", curObj);
        printf("raw data: ");
        dim_t n_obj_alloc = FLA_Obj_num_elem_alloc(curObj);
        double* obj_base_buffer = (double*)FLA_Obj_base_buffer(curObj);
        for(j = 0; j < n_obj_alloc; j++)
            printf("%.3f ", obj_base_buffer[j]);
        printf("\n");

        printf("A[%d] ", i);
        print_array("under permutation", curObj.order, permutation);

        dim_t* sizeObj = FLA_Obj_size(curObj);
        dim_t* strideObj = FLA_Obj_stride(curObj);

        FLA_Obj P;
        FLA_Obj_create_tensor(FLA_DOUBLE, curObj.order, sizeObj, strideObj, &P);
        FLA_Permute(curObj, permutation, &P);

        print_array("P stored permutation", P.order, P.permutation);
        print_array("P stored size", P.order, P.size);
        FLA_Obj_print_matlab("P", P);
        FLA_free(sizeObj);
        FLA_free(strideObj);
    }
}

int main(int argc, char* argv[]){

	FLA_Init();

	dim_t order;
	TLA_sym sym;
	dim_t* n;
	dim_t* b;
	dim_t* permutation;

	//Parse inputs
	if(parse_input(argc, argv, &order, &sym, &n, &b, &permutation) == FLA_FAILURE){
	    Usage();
	    FLA_Finalize();
	    return 0;
	}

	//Error check
	if(check_errors(order, sym, n, b, permutation) == FLA_FAILURE){
	    FLA_free(n);
	    FLA_free(b);
	    FLA_free(permutation);
	    FLA_Finalize();
	    return 0;
	}

	//Perform test
	FLA_Obj T;

	create_psym_tensor(order, sym, n, b, &T);
	FLA_Obj_print_matlab("T", T);
	

	test_permute_tensor(permutation, T);


    FLA_Obj_blocked_psym_tensor_free_buffer(&T);
    FLA_Obj_free_without_buffer(&T);

	FLA_free(n);
    FLA_free(b);
	FLA_Finalize();
	return 0;
}
