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
	dim_t block_size[order];
	for(i = 0; i < order; i++){
		flat_size[i] = nA;
		block_size[i] = bA;
	}

	dim_t blocked_stride[order];
	blocked_stride[0] = 1;

	for(i = 1; i < order; i++)
		blocked_stride[i] = flat_size[i-1]/bA * blocked_stride[i-1];

	FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, block_size, sym, A);
	FLA_Random_psym_tensor(*A);
}

void test_repart_routines(FLA_Obj A, dim_t nModes_repart, dim_t repart_modes[nModes_repart]){

    dim_t i;
    dim_t nParts = 1 << nModes_repart;
    dim_t nReparts = 1;
    for(i = 0; i < nModes_repart; i++)
        nReparts *= 3;

    FLA_Obj* Apart[nParts];
    for(i = 0; i < nParts; i++){
        Apart[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
    }
    FLA_Obj* Arepart[nReparts];
    for(i = 0; i < nReparts; i++){
        Arepart[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
    }

    dim_t part_sizes[nModes_repart];
    memset(&(part_sizes[0]), 0, nModes_repart * sizeof(dim_t));

    FLA_Side part_sides[nModes_repart];
    FLA_Side repart_sides[nModes_repart];
    for(i = 0; i < nModes_repart; i++){
        part_sides[i] = FLA_TOP;
        repart_sides[i] = FLA_BOTTOM;
    }

    FLA_Part_2powm(A, Apart, nModes_repart, repart_modes, part_sizes, part_sides);

    //Hack, need to make this better
    dim_t continueLooping = FALSE;
    for(i = 0; i < nModes_repart; i++)
        if(FLA_Obj_dimsize(*(Apart[0]),repart_modes[i]) < FLA_Obj_dimsize(A,repart_modes[i])){
            continueLooping = TRUE;
        }
    while(continueLooping){
        dim_t repart_sizes[nModes_repart];
        for(i = 0; i < nModes_repart; i++)
            repart_sizes[i] = 1;
        FLA_Repart_2powm_to_3powm(Apart, Arepart,
                                  nModes_repart, repart_modes,
                                  repart_sizes, repart_sides);
        /*********************************/

        for(i = 0; i < nParts; i++){
            printf("A[%d] ", i);
            print_array("offset", Apart[i]->order, (Apart[i]->offset));
            print_tensor("", *(Apart[i]));
            printf("\n");
        }
        /*********************************/
        FLA_Cont_with_3powm_to_2powm(Apart, Arepart,
                                     nModes_repart, repart_modes,
                                     part_sides);

        continueLooping = FALSE;
        for(i = 0; i < nModes_repart; i++)
            if(FLA_Obj_dimsize(*(Apart[0]),repart_modes[i]) < FLA_Obj_dimsize(A,repart_modes[i])){
                continueLooping = TRUE;
            }
    }

    for(i = 0; i < nParts; i++){
        printf("A[%d]:\n", i);
        print_tensor("", *(Apart[i]));
        printf("\n");
    }
    FLA_Merge_2powm( Apart, &A, nModes_repart, repart_modes);
    printf("A final:");
    print_tensor("", A);
    printf("\n");

    for(i = 0; i < nParts; i++)
        FLA_free(Apart[i]);
    for(i = 0; i < nReparts; i++)
        FLA_free(Arepart[i]);
}

int main(int argc, char* argv[]){
	dim_t i;

	FLA_Init();

	if(argc < 5){
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

	const dim_t nPartModes = atoi(argv[++argNum]);

    if(argc != 5 + sym.nSymGroups + m + 1 + nPartModes){
        Usage();
        FLA_Finalize();
        return 0;
    }

	dim_t partModes[nPartModes];
	for(i = 0; i < nPartModes; i++)
	    partModes[i] = atoi(argv[++argNum]);


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
	
	test_repart_routines(T, nPartModes, partModes);

	FLA_Obj_blocked_psym_tensor_free_buffer(&T);
	FLA_Obj_free_without_buffer(&T);

	FLA_Finalize();

	return 0;
}
