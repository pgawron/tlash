#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test psym tensor view operations.\n\n");
    printf("  tensor_psym <m> <nA> <bA> <nGroups> <group1Len> ... <groupNLen> <symMode1> ... <symModeN> <nPartModes> <partMode1> ... <partModeK>\n\n");
    printf("  m: order of symmetric tensor\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
	printf("  nGroups: number of symmetric groups\n");
	printf("  groupKLen: length of symmetric group K\n");
	printf("  symModeK: Kth mode in symmetric group ordering\n");
	printf("  nPartModes: number modes to partition\n");
	printf("  partModeK: Kth mode to partition\n");
}

void create_psym_tensor(dim_t order, dim_t nA, dim_t bA, TLA_sym sym, FLA_Obj* A){
	dim_t i;
	dim_t flat_size[FLA_MAX_ORDER];
	dim_t blocked_size[FLA_MAX_ORDER];
	dim_t block_size[FLA_MAX_ORDER];
	dim_t blocked_stride[FLA_MAX_ORDER];

	for(i = 0; i < order; i++){
		flat_size[i] = nA;
		block_size[i] = bA;
	}

	FLA_array_elemwise_quotient(order, flat_size, block_size, blocked_size);
	FLA_Set_tensor_stride(order, blocked_size, blocked_stride);

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, block_size, sym, A);
	FLA_Random_psym_tensor(*A);
}

void test_repart_routines(FLA_Obj A, dim_t nModes_repart, dim_t repart_modes[nModes_repart]){
    dim_t i;
    dim_t nParts = 1 << nModes_repart;
    dim_t nReparts;
    FLA_Obj** Apart;
    FLA_Obj** Arepart;
    dim_t* part_sizes;
    FLA_Side* part_sides;
    FLA_Side* repart_sides;
    dim_t* repart_sizes;

    nReparts = 1;
    for(i = 0; i < nModes_repart; i++)
        nReparts *= 3;

    Apart = (FLA_Obj**)FLA_malloc(nParts * sizeof(FLA_Obj*));
    Arepart = (FLA_Obj**)FLA_malloc(nReparts * sizeof(FLA_Obj*));
    part_sizes = (dim_t*)FLA_malloc(nModes_repart * sizeof(dim_t));
    part_sides = (FLA_Side*)FLA_malloc(nModes_repart * sizeof(FLA_Side));
    repart_sides = (FLA_Side*)FLA_malloc(nModes_repart * sizeof(FLA_Side));
    repart_sizes = (dim_t*)FLA_malloc(nModes_repart * sizeof(dim_t));

    TLA_create_part_obj(nParts, Apart);
    TLA_create_part_obj(nReparts, Arepart);

    memset(&(part_sizes[0]), 0, nModes_repart * sizeof(dim_t));

    for(i = 0; i < nModes_repart; i++){
        part_sides[i] = FLA_TOP;
        repart_sides[i] = FLA_BOTTOM;
        repart_sizes[i] = 1;
    }

    FLA_Part_2powm(A, Apart, nModes_repart, repart_modes, part_sizes, part_sides);

    while(FLA_Obj_dimsize(*(Apart[0]), repart_modes[0]) < FLA_Obj_dimsize(A, repart_modes[0])){
        FLA_Repart_2powm_to_3powm(Apart, Arepart,
                                  nModes_repart, repart_modes,
                                  repart_sizes, repart_sides);
        /*********************************/

        for(i = 0; i < nParts; i++){
            printf("A[%d] ", i);
            print_array("offset", Apart[i]->order, (Apart[i]->offset));
            FLA_Obj_print_matlab("", *(Apart[i]));
        }
        /*********************************/
        FLA_Cont_with_3powm_to_2powm(Apart, Arepart,
                                     nModes_repart, repart_modes,
                                     part_sides);
    }

    for(i = 0; i < nParts; i++){
        printf("A[%d]:\n", i);
        FLA_Obj_print_matlab("", *(Apart[i]));
    }

    FLA_Merge_2powm( Apart, &A, nModes_repart, repart_modes);
    FLA_Obj_print_matlab("A_final", A);

    TLA_destroy_part_obj(nParts, Apart);
    TLA_destroy_part_obj(nReparts, Arepart);

    FLA_free(Apart);
    FLA_free(Arepart);
    FLA_free(part_sizes);
    FLA_free(part_sides);
    FLA_free(repart_sides);
    FLA_free(repart_sizes);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* bA, dim_t* nPartModes, dim_t* partModes, TLA_sym* sym){
    dim_t i;
    int argNum = 0;

    if(argc < 5){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);

    sym->order = *order;
    sym->nSymGroups = atoi(argv[++argNum]);


    if(argc < 5 + sym->nSymGroups + *order){
    	return FLA_FAILURE;
    }

    for(i = 0; i < sym->nSymGroups; i++)
    	(sym->symGroupLens)[i] = atoi(argv[++argNum]);
    for(i = 0; i < *order; i++)
        	(sym->symModes)[i] = atoi(argv[++argNum]);

    *nPartModes = atoi(argv[++argNum]);

    if(argc != 5 + sym->nSymGroups + (*order) + 1 + (*nPartModes)){
        return FLA_FAILURE;
    }

    for(i = 0; i < *nPartModes; i++)
		partModes[i] = atoi(argv[++argNum]);


    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t bA, dim_t nPartModes, dim_t partModes[], TLA_sym sym){
	dim_t i;
	dim_t j;
	dim_t count;

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		return FLA_FAILURE;
	}
	if(order <= 0 || nA <= 0 || bA <= 0){
		printf("m, nA and bA must be greater than 0\n");
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

	if(nPartModes > order){
		printf("nPartModes must be less than order\n");
		return FLA_FAILURE;
	}

	for(i = 0; i < nPartModes; i++){
			if((i < nPartModes - 1) && (partModes[i+1] < partModes[i])){
				printf("partModes must be in ascending order\n");
				return FLA_FAILURE;
			}
			if(partModes[i] >= order){
				printf("partModes[i] must be less than order\n");
				return FLA_FAILURE;
			}
			for(j = 0; j < nPartModes; j++){
				if((i != j) && (partModes[i] == partModes[j])){
					printf("partModes[i] must not contain duplicates\n");
					return FLA_FAILURE;
				}
			}
		}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t bA;
	TLA_sym sym;
	dim_t nPartModes;
	dim_t partModes[FLA_MAX_ORDER];
	FLA_Obj T;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &bA, &nPartModes, partModes, &sym) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, bA, nPartModes, partModes, sym) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	create_psym_tensor(order, nA, bA, sym, &T);
	FLA_Obj_print_matlab("T", T);
	
	test_repart_routines(T, nPartModes, partModes);

	FLA_Obj_blocked_psym_tensor_free_buffer(&T);
	FLA_Obj_free_without_buffer(&T);

	FLA_Finalize();

	return 0;
}
