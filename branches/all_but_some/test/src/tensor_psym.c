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

void create_psym_tensor(dim_t order, TLA_sym sym, dim_t n[], dim_t b[], FLA_Obj* A){
	dim_t i, j;
	dim_t flat_size[FLA_MAX_ORDER];
	dim_t blocked_size[FLA_MAX_ORDER];
	dim_t blk_size[FLA_MAX_ORDER];
	dim_t blocked_stride[FLA_MAX_ORDER];

	for(i = 0; i < sym.nSymGroups; i++){
	    dim_t groupOffset = TLA_sym_group_mode_offset(sym, i);
	    for(j = 0; j < sym.symGroupLens[i]; j++){
	        flat_size[sym.symModes[groupOffset + j]] = n[i];
	        blk_size[sym.symModes[groupOffset + j]] = b[i];
	        blocked_size[sym.symModes[groupOffset + j]] = n[i] / b[i];
	    }
	}

	blocked_stride[0] = 1;
	for(i = 1; i < order; i++)
		blocked_stride[i] = blocked_size[i-1] * blocked_stride[i-1];

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, blk_size, sym, A);
	FLA_Random_psym_tensor(*A);
}

FLA_Error check_errors(dim_t order, TLA_sym sym, dim_t n[], dim_t b[], dim_t permutation[]){
    dim_t i, j;
    dim_t count;

    if(order <= 0){
        printf("m must be greater than 0\n");
        return FLA_FAILURE;
    }

	if(sym.nSymGroups <= 0){
		printf("nSymGroups must be greater than 0\n");
		return FLA_FAILURE;
	}

	count = 0;
	for(i = 0; i < sym.nSymGroups; i++){
		count += sym.symGroupLens[i];
		if(sym.symGroupLens[i] > order){
			printf("symGroupLens[i] must not be greater than order\n");
			return FLA_FAILURE;
		}
	}

	if(count != order){
		printf("Sum of symGroupLens must equal order\n");
		return FLA_FAILURE;
	}

	for(i = 0; i < order; i++){
		if((i < order - 1) && (sym.symModes[i+1] < sym.symModes[i])){
			printf("symModes must be in ascending order\n");
			return FLA_FAILURE;
		}

		if(sym.symModes[i] >= order){
			printf("symModes[i] must be less than order\n");
			return FLA_FAILURE;
		}
		for(j = 0; j < order; j++){
			if((i != j) && (sym.symModes[i] == sym.symModes[j])){
				printf("symModes[i] must not contain duplicates\n");
				return FLA_FAILURE;
			}
		}
	}

    for(i = 0; i < sym.nSymGroups; i++){
        if(n[i] <= 0){
            printf("n[%d] must be greater than 0\n", i);
            return FLA_FAILURE;
        }
        if(b[i] <= 0){
            printf("b[%d] must be greater than 0\n", i);
            return FLA_FAILURE;
        }
        if(n[i] % b[i] != 0){
            printf("b[%d] must evenly divide n[%d]\n", i, i);
            return FLA_FAILURE;
        }
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

FLA_Error parse_input(int argc, char* argv[], dim_t* order, TLA_sym* sym, dim_t* n, dim_t* b, dim_t* permutation){
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


    for(i = 0; i < sym->nSymGroups; i++)
        n[i] = atoi(argv[++argNum]);

    for(i = 0; i < sym->nSymGroups; i++)
        b[i] = atoi(argv[++argNum]);

    for(i = 0; i < *order; i++)
        permutation[i] = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

void test_permute_tensor(dim_t permutation[], FLA_Obj A){
    dim_t i,j;
    dim_t num_elem_alloc_A = FLA_Obj_num_elem_alloc(A);

	printf("Testing permutation routine\n");
    print_array("permutation", A.order, permutation);
    printf("-----------------------------\n");

    for(i = 0; i < num_elem_alloc_A; i++){
        FLA_Obj curObj = ((FLA_Obj*)FLA_Obj_base_buffer(A))[i];
        dim_t n_obj_alloc = FLA_Obj_num_elem_alloc(curObj);
        double* obj_base_buffer = (double*)FLA_Obj_base_buffer(curObj);
        dim_t* sizeObj = FLA_Obj_size(curObj);
        dim_t* strideObj = FLA_Obj_stride(curObj);
        FLA_Obj P;

        printf("A[%d] ", i);
        print_array("stored permutation", curObj.order, curObj.permutation);
        print_array("stored size", curObj.order, curObj.size);
        FLA_Obj_print_matlab("data", curObj);

        printf("raw data: ");
        for(j = 0; j < n_obj_alloc; j++)
            printf("%.3f ", obj_base_buffer[j]);
        printf("\n");

        printf("A[%d] ", i);
        print_array("under permutation", curObj.order, permutation);


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
	dim_t order;
	TLA_sym sym;
	dim_t n[FLA_MAX_ORDER];
	dim_t b[FLA_MAX_ORDER];
	dim_t permutation[FLA_MAX_ORDER];
	FLA_Obj T;

	FLA_Init();


	//Parse inputs
	if(parse_input(argc, argv, &order, &sym, n, b, permutation) == FLA_FAILURE){
	    Usage();
	    FLA_Finalize();
	    return 0;
	}

	//Error check
	if(check_errors(order, sym, n, b, permutation) == FLA_FAILURE){
		Usage();
	    FLA_Finalize();
	    return 0;
	}

	//Perform test

	create_psym_tensor(order, sym, n, b, &T);
	FLA_Obj_print_matlab("T", T);
	
	test_permute_tensor(permutation, T);

    FLA_Obj_blocked_psym_tensor_free_buffer(&T);
    FLA_Obj_free_without_buffer(&T);

    FLA_Finalize();
	return 0;
}
