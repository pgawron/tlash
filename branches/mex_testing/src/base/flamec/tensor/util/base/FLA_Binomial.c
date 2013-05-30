#include "FLAME.h"

dim_t binomial(dim_t n, dim_t k){
	dim_t i;
	dim_t r;
	r = 1;
	for(i = 1; i <= k; i++)
		r = (r * (n-k+i)) / i;
	return r;
}
