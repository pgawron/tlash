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


FLA_Error FLA_Permute_helper(FLA_Obj A, dim_t permutation[], FLA_Obj B, dim_t repart_mode_index){
	dim_t repartModeA = A.permutation[permutation[repart_mode_index]];
	dim_t repartModeB = repart_mode_index;

	//Down to a single column copy, call routine
	if(repart_mode_index == 0){
		TLA_Copy_col_mode(A, repartModeA, B, repartModeB);
		return FLA_SUCCESS;
	}

	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;

	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;

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

//TODO: Fix to account for initially permuted A objs
//NOTE: ONLY WORKS FOR FLA_SCALAR DATA
/*
FLA_Error FLA_Permute(FLA_Obj A, dim_t permutation[], FLA_Obj* B){
	if(FLA_Obj_elemtype(A) != FLA_SCALAR){
		printf("FLA_Permute: Only defined for objects with FLA_SCALAR elements");
		return FLA_SUCCESS;
	}
	
	dim_t i;
	dim_t order = FLA_Obj_order(A);
	for(i = 0; i < order; i++){
	    (B->permutation)[i] = i;
	    (B->size)[i] = A.size[A.permutation[permutation[i]]];
		(B->base->size)[i] = A.size[A.permutation[permutation[i]]];
	}
	(B->base->stride)[0] = 1;
	for(i = 1; i < order; i++)
	    (B->base->stride)[i] = (B->base->stride)[i-1] * ((B->size)[i-1]);

	return FLA_Permute_helper(A, permutation, *B, order - 1);
}
*/

FLA_Error FLA_Permute_single( FLA_Obj A, dim_t permutation[], FLA_Obj* B){
	
	dim_t i;
	dim_t order;
	dim_t* stride_A;
	dim_t* stride_B;
	dim_t* size_A;
	dim_t* size_B;
	
	order = FLA_Obj_order(A);
	stride_A = FLA_Obj_stride(A);
	stride_B = FLA_Obj_stride(*B);
	size_A = FLA_Obj_size(A);
	size_B = FLA_Obj_size(*B);
	dim_t ipermutation[A.order];
	
	for(i = 0; i < A.order; i++)
		ipermutation[permutation[i]] = i;
	
	double* buffer = (double*)FLA_Obj_base_buffer(*B);
	
    //Explicitly permute the data
	double* buf_A = (double*)FLA_Obj_base_buffer(A);
	dim_t curIndex[order];
	memset(&(curIndex[0]), 0, order * sizeof(dim_t));
	
	dim_t updatePtr = 0;
	dim_t linIndexFro = 0;
	dim_t linIndexTo = 0;
    while(TRUE){
		//Calculate linear index fro and to
		buffer[linIndexTo] = buf_A[linIndexFro];
		
		//Update
		//dim_t indexFro[A.order];
		//dim_t indexTo[A.order];
		//FLA_LinIndex_to_TIndex(A.order, stride_A, linIndexFro, indexFro);
		//FLA_LinIndex_to_TIndex(A.order, stride_B, linIndexTo, indexTo);

		/*
		printf("\n\nfro: %d to: %d\n", linIndexFro, linIndexTo);
		print_array("curIndex", A.order, curIndex);
		print_array("size A", A.order, size_A);
		print_array("size B", A.order, size_B);
		print_array("stride A", A.order, stride_A);
		print_array("stride B", A.order, stride_B);
		print_array("index fro", A.order, indexFro);
		print_array("index to ", A.order, indexTo);
		*/
		
		curIndex[updatePtr]++;
		linIndexFro += stride_A[updatePtr];
		linIndexTo += stride_B[ipermutation[updatePtr]];
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
		for(dim_t i = updatePtr-1; i < order; i--){
			curIndex[i] = 0;
			linIndexFro -= stride_A[i]*size_A[i];
			linIndexTo -= stride_B[ipermutation[i]]*size_B[ipermutation[i]];
		}
		updatePtr = 0;
	}
	
	FLA_free(stride_A);
	FLA_free(stride_B);
	FLA_free(size_A);
	FLA_free(size_B);
	
	return FLA_SUCCESS;
}

FLA_Error FLA_Permute(FLA_Obj A, dim_t permutation[], FLA_Obj* B){
	if(FLA_Obj_elemtype(A) != FLA_SCALAR){
		printf("FLA_Permute: Only defined for objects with FLA_SCALAR elements");
		return FLA_SUCCESS;
	}

	dim_t order = A.order;
	dim_t i;
	for(i = 0; i < order; i++){
	    (B->permutation)[i] = i;
	    (B->size)[i] = A.size[A.permutation[permutation[i]]];
		(B->base->size)[i] = A.size[A.permutation[permutation[i]]];
	}
	(B->base->stride)[0] = 1;
	for(i = 1; i < order; i++)
	    (B->base->stride)[i] = (B->base->stride)[i-1] * ((B->size)[i-1]);
	
	/*Check this*/
	dim_t full_perm[order];
	for(i = 0; i < order; i++)
		full_perm[i] = A.permutation[permutation[i]];
	return FLA_Permute_single(A, full_perm, B);
}	

