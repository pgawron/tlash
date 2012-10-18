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

#include "FLAME.h"

FLA_Error FLA_Sttsm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t maxIndex )
{
	dim_t j;
	printf("size_A:");
	for(j = 0; j < FLA_Obj_order(A); j++)
		printf(" %d", FLA_Obj_dimsize(A, j));
	printf("\n");
	printf("size_B:");
	for(j = 0; j < FLA_Obj_order(B); j++)
		printf(" %d", FLA_Obj_dimsize(B, j));
	printf("\n");
	printf("size_C:");
	for(j = 0; j < FLA_Obj_order(C); j++)
		printf(" %d", FLA_Obj_dimsize(C, j));
	printf("\n");
	
	if(mode == 0){
		FLA_Ttm_single(alpha, A, mode, beta, B, C);
	}else{
		dim_t i;
		FLA_Obj X;
		dim_t order = FLA_Obj_order( A );

		
		//FOR
		FLA_Obj BT, BB;
		FLA_Obj B0, B1, B2;
		FLA_Obj CT, CB;
		FLA_Obj C0, C1, C2;
		
		FLA_Part_1xmode2(B, &BT,
							&BB, 0, 0, FLA_TOP);	
		FLA_Part_1xmode2(C, &CT,
							&CB, mode, 0, FLA_TOP);	
		//Only symmetric part touched
		//Ponder this
		i = 0;
		while(i < maxIndex){
			//Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
			//Mode-0 of B matches Mode-n of C
			dim_t b = 1;
			FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
										/**/ /**/
											  &B1,
										  BB, &B2, 0, b, FLA_BOTTOM); 
			FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
										/**/ /**/
											  &C1,
										  CB, &C2, mode, b, FLA_BOTTOM);

			//Set up X to be corect size
			dim_t* size_X = FLA_Obj_size( A );
			//HACK I KNOW, just need to figure out how to make X the correct blocked size
			dim_t blkSize[order];
			for(i = 0; i < order; i++)
				blkSize[i] = (((FLA_Obj*)FLA_Obj_base_buffer(B1))[0]).size[0];
			for(i = 0; i < order; i++)
				size_X[i] *= blkSize[i];
			size_X[mode] = FLA_Obj_dimsize( B1, 0) * blkSize[mode];
			
			dim_t nElem = size_X[0];
			dim_t stride_X[order];
			stride_X[0] = 1;
			for(i = 1; i < order; i++){
				nElem *= size_X[i];
				stride_X[i] = stride_X[i-1]*size_X[i-1];
			}
			FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, size_X, stride_X, blkSize, &X);
			FLA_Set_zero_tensor(X);
			//End X setup

			printf("mode mult: %d input:\n", mode);
			FLA_Obj_print_tensor(A);
			FLA_Obj_print_tensor(B1);
			FLA_Ttm_single(alpha, A, mode, beta, B1, X);
			printf("mode mult: %d res:\n", mode);
			FLA_Obj_print_tensor(X);
			FLA_Sttsm_single(alpha, X, mode-1, beta, B, C1, i);

			FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
												   C1,
											/********/
											  &CB, C2, mode, FLA_TOP);
			FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
												   B1,
											/********/
											  &BB, B2, 0, FLA_TOP);
			i++;
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sttsm( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	dim_t order = FLA_Obj_order(A);
	FLA_Sttsm_single( alpha, A, order-1, beta, B, C, order);

	return FLA_SUCCESS;
}
