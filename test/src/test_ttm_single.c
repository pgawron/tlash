#include "FLAME.h"
#include "stdio.h"

void* genRandomData(dim_t order, dim_t size[]){
	dim_t i;
	dim_t numel = 1;
	double* buffer;

	for(i = 0; i < order; i++)
		numel *= size[i];

	buffer = (double*)FLA_malloc(numel * sizeof(double));
	for(i = 0; i < numel; i++)
		buffer[i] = ((double) rand() / ((double)RAND_MAX / 2.0F) ) - 1.0F;

	return (void*)buffer;
}

void* genSequentialData(dim_t order, dim_t size[]){
	dim_t i;
	dim_t numel = 1;
	double* buffer;

	for(i = 0; i < order; i++)
		numel *= size[i];

	buffer = (double*)FLA_malloc(numel * sizeof(double));
	for(i = 0; i < numel; i++)
		buffer[i] = (double)i;

	return (void*)buffer;
}

void initObj_ttm(dim_t order, dim_t size[], dim_t isRandom, FLA_Obj* obj){
  dim_t i;
  dim_t stride[FLA_MAX_ORDER];
  dim_t nData;
  double* data;

  stride[0] = 1;
  for(i = 1; i < order; i++)
	stride[i] = stride[i-1]*size[i-1];

  nData = 1;
  for(i = 0; i < order; i++)
    nData *= size[i];

  if(isRandom == 1)
	data = genRandomData(order, size);
  else
	data = genSequentialData(order, size);

  	FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size, obj);
  	FLA_Obj_attach_buffer_to_tensor(data, order, stride, obj);
}

void initObjZero_ttm(dim_t order, dim_t size[], FLA_Obj* obj){
  dim_t i;
  dim_t stride[FLA_MAX_ORDER];
  dim_t nData = 1;
  double* data;

  stride[0] = 1;
  for(i = 1; i < order; i++)
	stride[i] = stride[i-1]*size[i-1];


  for(i = 0; i < order; i++)
    nData *= size[i];

  data = (double*)FLA_malloc(nData * sizeof(double));;
  memset(data, 0, nData * sizeof(double));

  	FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size, obj);
  	FLA_Obj_attach_buffer_to_tensor(data, order, stride, obj);
}

void setToZero(FLA_Obj obj){
	dim_t i;
	dim_t order = FLA_Obj_order(obj);
	dim_t numel = 1;

	for(i = 0; i < order; i++)
		numel *= FLA_Obj_dimsize(obj, i);

	memset(FLA_Obj_base_buffer(obj), 0, numel * sizeof(dim_t));
}

void test_ttm_single(){

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj t, m, c;

	//Start setup params
	dim_t mode_mult = 1;
	dim_t tOrder = 2;
	dim_t tSize[] = {2, 2};
	dim_t mOrder = 2;
	dim_t mSize[] = {4, 2};
	dim_t cOrder = 2;
	dim_t cSize[] = {2, 4};
	//End setup

  printf("got here\n");
  initObj_ttm(tOrder, tSize, 0, &t);
  initObj_ttm(mOrder, mSize, 0, &m);
  initObjZero_ttm(cOrder, cSize, &c);


  printf("calling ttm\n");
  FLA_Ttm_single_mode(alpha, t, mode_mult, beta, m, c);

	printf("t tensor\n");
	FLA_Obj_print_tensor(t);
	
	printf("m matrix\n");
	FLA_Obj_print_tensor(m);

	printf("c tensor\n");
	FLA_Obj_print_tensor(c);
}
