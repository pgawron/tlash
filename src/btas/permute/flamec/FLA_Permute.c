/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2011, The University of Texas

   libflame is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   libflame is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with libflame; if you did not receive a copy, see
   http://www.gnu.org/licenses/.

   For more information, please contact us at flame@cs.utexas.edu or
   send mail to:

   Field G. Van Zee and/or
   Robert A. van de Geijn
   The University of Texas at Austin
   Department of Computer Sciences
   1 University Station C0500
   Austin TX 78712
*/

#include "FLA_Permute.h"


//FLAME way of permuting.... Very slow
FLA_Error FLA_Permute_helper(FLA_Obj A, const dim_t permutation[], FLA_Obj B, dim_t repart_mode_index){
	dim_t repartModeA = A.permutation[permutation[repart_mode_index]];
	dim_t repartModeB = repart_mode_index;

	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;

	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;

	//Down to a single column copy, call routine
	if(repart_mode_index == 0){
		TLA_Copy_col_mode(A, repartModeA, B, repartModeB);
		return FLA_SUCCESS;
	}



	FLA_Part_1xmode2(A, &AT,
						&AB, repartModeA, 0, FLA_TOP);

	FLA_Part_1xmode2(B, &BT,
						&BB, repartModeB, 0, FLA_TOP);

	while(FLA_Obj_dimsize(AT, repartModeA) < FLA_Obj_dimsize(A, repartModeA)){
		FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
										  &A1,
									  AB, &A2, repartModeA, 1, FLA_BOTTOM);
		FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
										  &B1,
									  BB, &B2, repartModeB, 1, FLA_BOTTOM);
		/***********************************/
			FLA_Permute_helper(A1, permutation, B1, repart_mode_index - 1);
		/***********************************/
		FLA_Cont_with_1xmode3_to_1xmode2(&AT, A0,
										  	  A1,
									  	 &AB, A2, repartModeA, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2(&BT, B0,
										  	  B1,
									  	 &BB, B2, repartModeB, FLA_TOP);
	}
	return FLA_SUCCESS;
}


//Permutes a tensor via permutation (loop-based)
FLA_Error FLA_Permute_single( FLA_Obj A, const dim_t permutation[], FLA_Obj* B){
	dim_t i;

	//FLA_Obj data
	dim_t order = FLA_Obj_order(A);
	dim_t stride_A[FLA_MAX_ORDER];
	dim_t stride_B[FLA_MAX_ORDER];
	dim_t size_A[FLA_MAX_ORDER];
	dim_t size_B[FLA_MAX_ORDER];
	double* buf_A;
	double* buf_B;
	
	//Inverse Permutation data
	dim_t ipermutation[FLA_MAX_ORDER];

	//Loop data
	dim_t curIndex[FLA_MAX_ORDER];
	dim_t updatePtr;
	dim_t linIndexFro;
	dim_t linIndexTo;

	//Init Obj data
	memcpy(&(size_A[0]), &(A.size[0]), order * sizeof(dim_t));
	memcpy(&(stride_A[0]), &((A.base->stride)[0]), order * sizeof(dim_t));
	memcpy(&(size_B[0]), &((B->size)[0]), order * sizeof(dim_t));
	memcpy(&(stride_B[0]), &((B->base->stride)[0]), order * sizeof(dim_t));


	//Init iperm data
	for(i = 0; i < A.order; i++)
		ipermutation[permutation[i]] = i;
	
	buf_A = (double*)FLA_Obj_base_buffer(A);
	buf_B = (double*)FLA_Obj_base_buffer(*B);
	
    //Init loop data
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	updatePtr = 0;
	linIndexFro = 0;
	linIndexTo = 0;
	
	//Loop over all entries, and explicitly permute the data
    while(TRUE){
		//Calculate linear index fro and to
		buf_B[linIndexTo] = buf_A[linIndexFro];
		
		//Update index pointer and linIndices Fro and To
		curIndex[updatePtr]++;
		linIndexFro += stride_A[updatePtr];
		linIndexTo += stride_B[ipermutation[updatePtr]];

		//Loop update
		while(updatePtr < order && curIndex[updatePtr] == size_A[updatePtr]){
			updatePtr++;
			if(updatePtr < order){
				curIndex[updatePtr]++;
				linIndexFro += stride_A[updatePtr];
				linIndexTo += stride_B[ipermutation[updatePtr]];
			}
		}
		if(updatePtr >= order)
			break;
		for(i = updatePtr-1; i < order; i--){
			curIndex[i] = 0;
			linIndexFro -= stride_A[i]*size_A[i];
			linIndexTo -= stride_B[ipermutation[i]]*size_B[ipermutation[i]];
		}
		updatePtr = 0;
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Permute(FLA_Obj A, const dim_t permutation[], FLA_Obj* B){
	//Account for A having a permutation of its own
	dim_t i;
	dim_t order = FLA_Obj_order(A);
	dim_t full_perm[FLA_MAX_ORDER];

	if(FLA_Obj_elemtype(A) != FLA_SCALAR){
		printf("FLA_Permute: Only defined for objects with FLA_SCALAR elements");
		return FLA_SUCCESS;
	}


	for(i = 0; i < order; i++)
		full_perm[i] = A.permutation[permutation[i]];

	//Update data of B to reflect result of permutation
	for(i = 0; i < order; i++){
	    (B->permutation)[i] = i;
	    (B->size)[i] = A.size[full_perm[i]];
		(B->base->size)[i] = A.size[full_perm[i]];
	}
	FLA_Set_tensor_stride(order, B->size, B->base->stride);

	return FLA_Permute_single(A, full_perm, B);
}	

