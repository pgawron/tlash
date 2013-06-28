#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test TLASH Part operations.\n\n");
    printf("  tlash_part <m> <nA> <bA>\n\n");
    printf("  m: order of symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
}

void initBlockedTensor(dim_t order, dim_t size[], dim_t bSize[], FLA_Obj* obj){
  dim_t stride[FLA_MAX_ORDER];
  dim_t blked_size[FLA_MAX_ORDER];

  FLA_array_elemwise_quotient(order, size, bSize, blked_size);
  FLA_Set_tensor_stride(order, blked_size, stride);

  FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size, stride, bSize, obj);
  FLA_Random_tensor(*obj);
}

void print_Part(FLA_Obj A, FLA_Obj* parts[]){
	dim_t i;
	dim_t order = FLA_Obj_order(A);
	dim_t curIndex[FLA_MAX_ORDER];
	dim_t update_ptr = 0;
	dim_t curPart = 0;

	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	while(TRUE){
		printf("A");
		for(i = 0; i < order; i++)
			printf("%d", curIndex[i]);
		printf(" = [");
		FLA_Obj_print_tensor(*(parts[curPart++]));
		printf("];\n");
		//Update
		curIndex[update_ptr]++;
		if(update_ptr < order && curIndex[update_ptr] == 2){
			update_ptr++;
			curIndex[update_ptr]++;
			while(update_ptr < order && curIndex[update_ptr] == 2){
				update_ptr++;
				curIndex[update_ptr]++;
			}
		}
		if(update_ptr >= order)
			break;
		for(i = 0; i < update_ptr; i++)
			curIndex[i] = 0;
		update_ptr = 0;
	}
}

void test_tlash_part_routines(dim_t m, dim_t nA, dim_t bA){
  dim_t i;
  FLA_Obj A;
  dim_t sizeA[FLA_MAX_ORDER];
  dim_t sizeBlock[FLA_MAX_ORDER];
  dim_t partSize[FLA_MAX_ORDER];
  FLA_Side sides[FLA_MAX_ORDER];
  dim_t nPartitions;
  FLA_Obj** Apart;
  dim_t nModes_repart = m;
  dim_t repart_modes[FLA_MAX_ORDER];

  for(i = 0; i < m; i++){
	sizeA[i] = nA;
	sizeBlock[i] = bA;
	partSize[i] = 1;
	sides[i] = FLA_TOP;
  }

  initBlockedTensor(m, sizeA, sizeBlock, &(A) );

  nPartitions = (1 << m);

  Apart = (FLA_Obj**)FLA_malloc(nPartitions * sizeof(FLA_Obj*));

  TLA_create_part_obj(nPartitions, Apart);

  printf("A = [\n");
  FLA_Obj_print_tensor(A);
  printf("];\n");

  printf("Partitioning A\n");

  for(i = 0; i < m; i++)
      repart_modes[i] = i;
  FLA_Part_2powm(A, Apart, nModes_repart, repart_modes, partSize, sides);

  print_Part(A, Apart);

  TLA_destroy_part_obj(nPartitions, Apart);

  FLA_free(Apart);
  FLA_Obj_blocked_tensor_free_buffer(&A);
  FLA_Obj_free_without_buffer(&A);
}

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* bA){
    int argNum = 0;

    if(argc != 4){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *bA = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t bA){

	if(order <= 0 || nA <= 0 || bA <= 0){
		printf("m, nA, and bA must be greater than 0\n");
		return FLA_FAILURE;
	}

	if(nA % bA != 0){
		printf("bA must evenly divide nA\n");
		return FLA_FAILURE;
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t bA;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &bA) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, bA) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_tlash_part_routines(order, nA, bA);

	FLA_Finalize();

	return 0;
}
