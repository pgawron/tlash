#include "FLAME.h"
#include "stdio.h"

void Usage()
{
    printf("Test sttsm operation.\n\n");
    printf("  test_sttsm <m> <nA> <nC>\n\n");
    printf("  m: order of symmetric tensors\n");
    printf("  nA: mode-length of tensor A\n");
    printf("  nC: mode-length of symmetric tensor C\n");
}

void initSymmTensor(dim_t order, dim_t size[], FLA_Obj* obj){
  dim_t i;
  dim_t blocked_stride[FLA_MAX_ORDER];
  TLA_sym sym;

  blocked_stride[0] = 1;
  for(i = 1; i < order; i++)
      blocked_stride[i] = blocked_stride[i-1] * (size[i-1] / size[0]);


  sym.order = order;
  sym.nSymGroups = 1;
  sym.symGroupLens[0] = sym.order;
  for(i = 0; i < sym.order; i++)
      (sym.symModes)[i] = i;
  FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, size, blocked_stride, size, sym, obj);

  FLA_Random_psym_tensor(*obj);

}

void initMatrix(dim_t size[2], FLA_Obj* obj){
  dim_t order = 2;
  dim_t strideObj[] = {1, size[0]};

  FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size, strideObj, size, obj);
  FLA_Random_tensor(*obj);
}

void test_sttsm(int m, int nA, int nC, double* elapsedTime){
	//Setup parameters
  dim_t i;
  dim_t aSize[FLA_MAX_ORDER];
  dim_t bSize[] = {nC, nA};
  dim_t cSize[FLA_MAX_ORDER];
  double startTime;
  double endTime;

  for(i = 0; i < m; i++)
    aSize[i] = nA;
  for(i = 0; i < m; i++)
    cSize[i] = nC;
  //End setup parameters

  FLA_Obj alpha = FLA_ONE;
  FLA_Obj beta = FLA_ONE;

  FLA_Obj A, B, C;

  initSymmTensor(m, aSize, &A);
  initMatrix(bSize, &B);
  initSymmTensor(m, cSize, &C);

  FLA_Obj_print_matlab("A", A);
  FLA_Obj_print_matlab("B", B);
  FLA_Obj_print_matlab("C", C);


  startTime = FLA_Clock();
  FLA_Sttsm_without_psym_temps(alpha, A, beta, B, C);
  //FLA_Ttm(alpha, A, m, modes, beta, Barr, C);
  endTime = FLA_Clock();
  *elapsedTime = endTime - startTime;
//printf("end computation\n");
//	printf("c tensor\n");
//	FLA_Obj_print_flat_tensor(c);


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

FLA_Error parse_input(int argc, char* argv[], dim_t* order, dim_t* nA, dim_t* nC){
    int argNum = 0;

    if(argc != 4){
        return FLA_FAILURE;
    }

    *order = atoi(argv[++argNum]);
    *nA = atoi(argv[++argNum]);
    *nC = atoi(argv[++argNum]);

    return FLA_SUCCESS;
}

FLA_Error check_errors(dim_t order, dim_t nA, dim_t nC){

	if(order <= 0 || nA <= 0 || nC <= 0){
		printf("m, nA, and nC must be greater than 0\n");
		return FLA_FAILURE;
	}

	return FLA_SUCCESS;
}

int main(int argc, char* argv[]){
	dim_t order;
	dim_t nA;
	dim_t nC;
	double elapsedTime;

	FLA_Init();

	if(parse_input(argc, argv, &order, &nA, &nC) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	if(check_errors(order, nA, nC) == FLA_FAILURE){
		Usage();
		FLA_Finalize();
		return 0;
	}

	test_sttsm(order, nA, nC, &elapsedTime);

	FLA_Finalize();

	//Print out results
//	double gflops = Sttsm_GFlops(m, nA, nC, elapsedTime);
	printf("ARGS DENSE %d %d %d %d %d\n", order, nA, nC, nA, nC);
	printf("TIME %.6f\n", elapsedTime);
	printf("GFLOPS %.6f\n", -1.0);
	return 0;
}
