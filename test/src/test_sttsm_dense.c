#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sttsm operation.\n\n");
    printf("  test_sttsm <m> <nA> <nC><b>\n\n");
    printf("  m: order of hyper-symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  nC: mode-length of hyper-symmetric tensor C\n");
}

void initSymmTensor(dim_t order, dim_t size[order], FLA_Obj* obj){
  dim_t i;
  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size, obj);
  dim_t numel = 1;
  for(i = 0; i < order; i++)
	numel *= size[i];

  dim_t stride[order];
  stride[0] = 1;
  for(i = 1; i < order; i++)
	stride[i] = size[i-1] * stride[i-1];

  double* buf = (double*)FLA_malloc(numel * sizeof(double));
  FLA_Obj_attach_buffer_to_tensor(buf, order, stride, obj);

  dim_t** symmetries = (dim_t**)FLA_malloc(1 * sizeof(dim_t*));
  symmetries[0] = (dim_t*)FLA_malloc(order * sizeof(dim_t));
  for(i = 0; i < order; i++)
	symmetries[0][i] = i;

  dim_t size_symmetries[] = {order};
 /*
  dim_t sizeB[] = {size[0],1};
  dim_t stride_B[] = {1, size[0]};

  FLA_Obj B;
  FLA_Obj_create_tensor(FLA_DOUBLE, 2, size_B, stride_B, &B);
  FLA_Random_matrix(B);

  FLA_Obj Barr[order];
  for(i = 0; i < order;i++)
	Barr[i] = B; 

  FLA_Obj A;
  dim_t sizeA[order];
  for(i = 0; i < order; i++)
	sizeA[i] = 1;
  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order; sizeA, &A);
  double* dataA = (double*)FLA_malloc(1 * sizeof(double));
  dataA[0] = 1;
  dim_t strideA[order];
  for(i = 0; i < order; i++)
	strideA[i] = 1;
  FLA_Obj_attach_buffer_to_tensor(dataA, order, strideA, &A);

  FLA_Ttm(FLA_ONE, A, order, symmetries, FLA_ONE, Barr, *obj);
*/
  FLA_Random_dense_symm_tensor(1, size_symmetries, symmetries, obj);
  FLA_free(symmetries[0]);
  FLA_free(symmetries);
}

void initMatrix(dim_t size[2], FLA_Obj* obj){
  dim_t order = 2;
  dim_t sizeObj[] = {size[0], size[1]};
  dim_t strideObj[] = {1, sizeObj[0]};

  FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, sizeObj, obj);
  obj->base->elemtype = FLA_SCALAR;

  double* buf = (double*)FLA_malloc(sizeObj[0] * sizeObj[1] * sizeof(double));
  FLA_Obj_attach_buffer_to_tensor(buf, order, strideObj, obj);

  FLA_Adjust_2D_info(obj);

  FLA_Random_matrix(*obj);
}

void setSymmTensorToZero(FLA_Obj obj){
    dim_t i;
	dim_t order = FLA_Obj_order(obj);
	dim_t numel = 1;
	for(i = 0; i < order; i++)
		numel *= FLA_Obj_dimsize(obj, i);

	memset(&(((double*)FLA_Obj_base_buffer(obj))[0]), 0, numel * sizeof(double));
}

void test_sttsm(int m, int nA, int nC, double* elapsedTime){
	//Setup parameters
  dim_t i;
  dim_t aSize[m];
  for(i = 0; i < m; i++)
    aSize[i] = nA;
  dim_t bSize[] = {nC, nA};
  dim_t cSize[m];
  for(i = 0; i < m; i++)
  cSize[i] = nC;
  //End setup parameters

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj A, B, C;

  initSymmTensor(m, aSize, &A);

  initMatrix(bSize, &B);

  initSymmTensor(m, cSize, &C);
  setSymmTensorToZero(C);

  dim_t modes[m];
  for(i = 0; i < m; i++)
	modes[i] = i;

//printf("begin computation\n");
  FLA_Obj Barr[m];
  for(i = 0; i < m; i++)
	Barr[i] = B;

  double startTime = FLA_Clock();
  FLA_Ttm(alpha, A, m, modes, beta, Barr, C);
  double endTime = FLA_Clock();
  *elapsedTime = endTime - startTime;
//printf("end computation\n");
//	printf("c tensor\n");
//	FLA_Obj_print_flat_tensor(c);

  FLA_Obj_free_buffer(&B);
  FLA_Obj_free_without_buffer(&B);

  FLA_Obj_free_buffer(&A);
  FLA_Obj_free_without_buffer(&A);
  FLA_Obj_free_buffer(&C);
  FLA_Obj_free_without_buffer(&C);
}

double Sttsm_GFlops(dim_t m, dim_t nA, dim_t nC, double elapsedTime){
	dim_t i;
	dim_t KpowTerm = 1;
	dim_t NpowTerm = 1;

	for(i = 0; i < m; i++){
		KpowTerm *= nC;
		NpowTerm *= nA;
	}

	return nC*nA*(KpowTerm-NpowTerm)/(nC-nA)/(1.e9*elapsedTime);
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
	const int nC = atoi(argv[++argNum]);

	if(m <= 0){
		printf("m must be greater than 0\n");
		FLA_Finalize();
		return 0;
	}

	double elapsedTime;
	test_sttsm(m, nA, nC, &elapsedTime);

	FLA_Finalize();

	//Print out results
	double gflops = Sttsm_GFlops(m, nA, nC, elapsedTime);
	printf("ARGS DENSE %d %d %d %d %d\n", m, nA, nC, nA, nC);
	printf("TIME %.6f\n", elapsedTime);
	printf("GFLOPS %.6f\n", gflops);
	return 0;
}
