#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test TLASH Part, Repart, Cont_with, Merge operations.\n\n");
    printf("  tlash_part <m> <nA> <bA>\n\n");
    printf("  m: order of symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  bA: mode-length of block of A\n");
}

void initBlockedTensor(dim_t order, dim_t size[order], dim_t bSize[order], FLA_Obj* obj){
  dim_t i;
  dim_t stride[order];
  stride[0] = 1;
  for(i = 1; i < order; i++)
	stride[i] = size[i-1]*stride[i-1];

  FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size, stride, bSize, obj);
  FLA_Random_tensor(*obj);
}

void print_Part(FLA_Obj A, FLA_Obj parts[]){
	dim_t i;
	dim_t order = FLA_Obj_order(A);

	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	dim_t update_ptr = 0;
	dim_t curPart = 0;
	while(TRUE){
		printf("A");
		for(i = 0; i < order; i++)
			printf("%d", curIndex[i]);
		printf(" = [");
		FLA_Obj_print_tensor(parts[curPart++]);
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
  dim_t sizeA[m];
  dim_t sizeBlock[m];
  dim_t partSize[m];
  FLA_Side sides[m];
  for(i = 0; i < m; i++){
	sizeA[i] = nA;
	sizeBlock[i] = bA;
	partSize[i] = 1;
	sides[i] = FLA_TOP;
  }

  initBlockedTensor(m, sizeA, sizeBlock, &(A) );

  dim_t nPartitions = 1;
  nPartitions <<= m;

  FLA_Obj Apart[nPartitions];

  printf("A = [\n");
  FLA_Obj_print_tensor(A);
  printf("];\n");

  printf("Partitioning A\n");
  dim_t nModes_repart = m;
  dim_t repart_modes[m];
  for(i = 0; i < m; i++)
      repart_modes[i] = i;
  FLA_Part_2powm(A, Apart, nModes_repart, repart_modes, partSize, sides);

  print_Part(A, Apart);

}

int main(int argc, char* argv[]){
	FLA_Init();

	if(argc < 4){
		Usage();
		FLA_Finalize();
		return 0;
	}
		
	int argNum = 0;
	const int m = atoi(argv[++argNum]);
	const int nA = atoi(argv[++argNum]);
	const int bA = atoi(argv[++argNum]);

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

	test_tlash_part_routines(m, nA, bA);

	FLA_Finalize();

	return 0;
}
