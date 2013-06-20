#include "FLAME.h"

void print_array(const char* header, dim_t nElems, dim_t arr[nElems]){
 	dim_t i;

	printf("%s:", header);
	for(i = 0; i < nElems; i++){
		printf(" %d", arr[i]);
	}
	printf("\n");
}
