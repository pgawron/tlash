#include "FLAME.h"
#include "stdio.h"
#include "test_ttm_single.h"

int main(){
  FLA_Init();

//  test_ttm_single();
//  test_symm_tlash();

  test_sttsm();

  FLA_Finalize();
  return 0;
}
