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
  dim_t blocked_stride[order];
  blocked_stride[0] = 1;
  for(i = 1; i < order; i++)
      blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / size[0]);

  TLA_sym sym;
  sym.order = FLA_Obj_order(*obj);
  sym.nSymGroups = 1;
  sym.symGroupLens[0] = sym.order;
  for(i = 0; i < sym.order; i++)
      (sym.symModes)[i] = i;
  FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size, blocked_stride, size, sym, obj);

  FLA_Random_psym_tensor(*obj);

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
	printf("A tensor\n");
	FLA_Obj_print_flat_tensor(A);

	printf("B tensor\n");
	FLA_Obj_print_flat_tensor(B);

	printf("c tensor\n");
	FLA_Obj_print_flat_tensor(C);


  FLA_Obj_free_buffer(&B);
  FLA_Obj_free_without_buffer(&B);

  FLA_Obj_blocked_psym_tensor_free_buffer(&A);
  FLA_Obj_free_without_buffer(&A);
  FLA_Obj_blocked_psym_tensor_free_buffer(&C);
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
//	double gflops = Sttsm_GFlops(m, nA, nC, elapsedTime);
	printf("ARGS DENSE %d %d %d %d %d\n", m, nA, nC, nA, nC);
	printf("TIME %.6f\n", elapsedTime);
	printf("GFLOPS %.6f\n", -1.0);
	return 0;
}
