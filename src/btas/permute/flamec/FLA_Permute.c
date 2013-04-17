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
