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

FLA_Error FLA_Sttsm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex, FLA_Obj* temps[] )
{
	//FOR
	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;
	FLA_Obj CT, CB;
	FLA_Obj C0, C1, C2;
	dim_t loopCount;

	FLA_Part_1xmode2(B, &BT,
						&BB, 0, 0, FLA_TOP);
	FLA_Part_1xmode2(C, &CT,
						&CB, mode, 0, FLA_TOP);
	//Only symmetric part touched
	//Ponder this
	loopCount = 0;
	while(loopCount <= endIndex){
		dim_t b = 1;
		FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
									/**/ /**/
										  &B1,
									  BB, &B2, 0, b, FLA_BOTTOM);
		FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
									/**/ /**/
										  &C1,
									  CB, &C2, mode, b, FLA_BOTTOM);
		/*********************************/
		//Make sure that if we are at the bottom of the recursion,
		//We multiply into C (not X)
		if(mode == 0){
			FLA_Ttm_single_mode(alpha, A, mode, beta, B1, C1);
		}else{
			//Initialize X to 0
			FLA_Obj X = *(temps[mode]);
			FLA_Set_zero_tensor(X);

			//Compute X
			FLA_Ttm_single_mode(alpha, A, mode, beta, B1, X);
			//Use X for rest of computation
			FLA_Sttsm_single(alpha, X, mode-1, beta, B, C1, loopCount, temps);
		}
		/*********************************/
		FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
											   C1,
										/********/
										  &CB, C2, mode, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
											   B1,
										/********/
										  &BB, B2, 0, FLA_TOP);
		loopCount++;
	}


	return FLA_SUCCESS;
}

FLA_Error FLA_Sttsm_single_psttm( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex, FLA_Obj* temps[] )
{
    FLA_Obj BT, BB;
    FLA_Obj B0, B1, B2;
    FLA_Obj CT, CB;
    FLA_Obj C0, C1, C2;

    dim_t loopCount;

    FLA_Part_1xmode2(B, &BT,
                        &BB, 0, 0, FLA_TOP);
    FLA_Part_1xmode2(C, &CT,
                        &CB, mode, 0, FLA_TOP);
    //Only symmetric part touched
    //Ponder this
    loopCount = 0;
    while(loopCount <= endIndex){
        dim_t b = 1;
        FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
                                    /**/ /**/
                                          &B1,
                                      BB, &B2, 0, b, FLA_BOTTOM);
        FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
                                    /**/ /**/
                                          &C1,
                                      CB, &C2, mode, b, FLA_BOTTOM);
		/**************************************************/
        if(mode == 0){
            FLA_Ttm_single_mode(alpha, A, mode, beta, B1, C1);
        }else{
			//Set X to 0
			FLA_Obj X = *(temps[mode]);
            FLA_Set_zero_tensor(X);

            //Compute X
            FLA_Psttm(alpha, A, mode, beta, B1, X);

            //Use X for rest of computation
            FLA_Sttsm_single_psttm(alpha, X, mode-1, beta, B, C1, loopCount, temps);
        }
		/**************************************************/
        FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
                                               C1,
                                        /********/
                                          &CB, C2, mode, FLA_TOP);
        FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
                                               B1,
                                        /********/
                                          &BB, B2, 0, FLA_TOP);
        loopCount++;
    }

    return FLA_SUCCESS;
}

void initialize_psym_temporaries(FLA_Obj A, FLA_Obj C, FLA_Obj* temps[]){
	dim_t i, j;
	dim_t order = FLA_Obj_order(A);
	dim_t temp_blocked_size[FLA_MAX_ORDER];
	dim_t temp_blocked_stride[FLA_MAX_ORDER];
	dim_t temp_block_size[FLA_MAX_ORDER];
	dim_t temp_flat_size[FLA_MAX_ORDER];
	FLA_Obj* buf_A = (FLA_Obj*)FLA_Obj_base_buffer(A);
	FLA_Obj* buf_C = (FLA_Obj*)FLA_Obj_base_buffer(C);
	
	TLA_sym tmpSym = A.sym;
	for(i = order - 1; i > 0; i--){
		//Symmetry is different for each temporary
		TLA_sym Xsym;
		TLA_split_sym_group(tmpSym, 1, &i, &Xsym);

		//For simplicity, temp_block_size has all same size as A first
		memcpy(&(temp_block_size[0]), &(buf_A[0].size[0]), order * sizeof(dim_t));
		memcpy(&(temp_blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
		//Adjust output modes
		for(j = i; j < order; j++){
			//Each temporary is multiplied by a block vector
			temp_blocked_size[j] = 1;
			temp_block_size[j] = buf_C[0].size[j];
		}

		//Adjust the flat_size & stride
		FLA_array_elemwise_product(order, temp_blocked_size, temp_block_size, temp_flat_size);
		FLA_Set_tensor_stride(order, temp_blocked_size, temp_blocked_stride);

		//Create the temporary
		temps[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
		FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, temp_flat_size, temp_blocked_stride, temp_block_size, Xsym, temps[i]);
		
		tmpSym = Xsym;
	}
}

void destroy_psym_temporaries(dim_t order, FLA_Obj* temps[]){
	dim_t i;
	for(i = order - 1; i > 0; i--){
		FLA_Obj_blocked_psym_tensor_free_buffer(temps[i]);
        FLA_Obj_free_without_buffer(temps[i]);
		FLA_free(temps[i]);
	}
}

void initialize_temporaries(FLA_Obj A, FLA_Obj C, FLA_Obj* temps[]){
	dim_t i, j;
	dim_t order = FLA_Obj_order(A);
	dim_t temp_blocked_size[FLA_MAX_ORDER];
	dim_t temp_blocked_stride[FLA_MAX_ORDER];
	dim_t temp_block_size[FLA_MAX_ORDER];
	dim_t temp_flat_size[FLA_MAX_ORDER];
	FLA_Obj* buf_A = (FLA_Obj*)FLA_Obj_base_buffer(A);
	FLA_Obj* buf_C = (FLA_Obj*)FLA_Obj_base_buffer(C);

	//Set up each temporary
	for(i = order - 1; i > 0; i--){
		//For simplicity, temp_block_size has all same size as A first
		memcpy(&(temp_block_size[0]), &(buf_A[0].size[0]), order * sizeof(dim_t));
		memcpy(&(temp_blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
		//Adjust output modes
		for(j = i; j < order; j++){
			//Each temporary is multiplied by a block vector
			temp_blocked_size[j] = 1;
			temp_block_size[j] = buf_C[0].size[j];
		}

		//Adjust the flat_size & stride
		FLA_array_elemwise_product(order, temp_blocked_size, temp_block_size, temp_flat_size);
		FLA_Set_tensor_stride(order, temp_blocked_size, temp_blocked_stride);

		//Create the temporary
		temps[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
		FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, temp_flat_size, temp_blocked_stride, temp_block_size, temps[i]);
	}
}

void destroy_temporaries(dim_t order, FLA_Obj* temps[]){
	dim_t i;
	for(i = order - 1; i > 0; i--){
		FLA_Obj_blocked_tensor_free_buffer(temps[i]);
        FLA_Obj_free_without_buffer(temps[i]);
		FLA_free(temps[i]);
	}
}

//No psym temps
FLA_Error FLA_Sttsm_without_psym_temps( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	//Create temporaries used in sttsm
	FLA_Obj* temps[FLA_MAX_ORDER];
	initialize_temporaries(A, C, temps);
	
	//Compute
	FLA_Sttsm_single( alpha, A, FLA_Obj_order(C)-1, beta, B, C, FLA_Obj_dimsize(C,FLA_Obj_order(C)-1)-1, temps);
	
	//Cleanup
	destroy_temporaries(A.order, temps);

	return FLA_SUCCESS;
}


//Using psym temps

FLA_Error FLA_Sttsm_with_psym_temps( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	//Create temporaries used in sttsm
	FLA_Obj* temps[FLA_MAX_ORDER];
	initialize_psym_temporaries(A, C, temps);

	//Compute
    FLA_Sttsm_single_psttm( alpha, A, FLA_Obj_order(C)-1, beta, B, C, FLA_Obj_dimsize(C,FLA_Obj_order(C)-1)-1, temps);
	
    //Cleanup
	destroy_psym_temporaries(A.order, temps);

	return FLA_SUCCESS;
}
