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
	//dim_t order = FLA_Obj_order( A );

	if(mode == 0){
		FLA_Obj BT, BB;
		FLA_Obj B0, B1, B2;
		FLA_Obj CT, CB;
		FLA_Obj C0, C1, C2;

		FLA_Part_1xmode2(B, &BT,
						 &BB, 0, 0, FLA_TOP);
		FLA_Part_1xmode2(C, &CT,
						 &CB, mode, 0, FLA_TOP);
		dim_t loopCount = 0;
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

//			printf("ttm performed: %d\n", FLA_Ttm_Ops(order, A.size, B1.size, mode));
			FLA_Ttm_single_mode(alpha, A, mode, beta, B1, C1);

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
	}else{
		//dim_t i;
		//FLA_Obj X;


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
		dim_t loopCount = 0;
		while(loopCount <= endIndex){
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
			
			/*
			//Set up X to be corect size
			dim_t size_X[FLA_Obj_order( A )];
			//HACK I KNOW, just need to figure out how to make X the correct blocked size
			dim_t blkSize[order];
			for(i = 0; i < order; i++)
				blkSize[i] = (((FLA_Obj*)FLA_Obj_base_buffer(A))[0]).size[i];
			blkSize[mode] = (((FLA_Obj*)FLA_Obj_base_buffer(B1))[0]).size[0];

			for(i = 0; i < order; i++)
				size_X[i] = blkSize[i] * FLA_Obj_dimsize(A,i);
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
			 */
			FLA_Obj X = *(temps[mode]);
            FLA_Set_zero_tensor(X);
			
			//End X setup

			/*
			printf("\n\nMultiply happening\n");
			printf("a = tensor([");
			FLA_Obj_print_tensor(A);
			printf("],[");
			for(i = 0; i < FLA_Obj_order(A); i++)
				printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(A)))[0],i) * FLA_Obj_dimsize(A,i));
			printf("]);\n\n");

			printf("\ntimes %d\n", mode);
			printf("b1 = reshape([");
			FLA_Obj_print_tensor(B1);
			printf("],[");
			printf("%d ", 2);
			printf("%d ", 4);
			printf("]);\n\n");

			printf("preX = tensor([");
			FLA_Obj_print_tensor(X);
			printf("],[");
			for(i = 0; i < FLA_Obj_order(X); i++)
				printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(X)))[0],i) * FLA_Obj_dimsize(X,i));
			printf("]);\n\n");
*/
//			FLA_Obj_print_matlab("X", X);
			FLA_Ttm_single_mode(alpha, A, mode, beta, B1, X);
/*
			printf("postX = tensor([");
			FLA_Obj_print_tensor(X);
			printf("],[");
			for(i = 0; i < FLA_Obj_order(X); i++)
				printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(X)))[0],i) * FLA_Obj_dimsize(X,i));
			printf("]);\n\n");
*/

			FLA_Sttsm_single(alpha, X, mode-1, beta, B, C1, loopCount, temps);

//			FLA_Obj_blocked_free_buffer(&X);
//			FLA_Obj_free_without_buffer(&X);


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
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sttsm_single_psttm( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex, FLA_Obj* temps[] )
{
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
    dim_t loopCount = 0;
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
			//Set up X temporary
			FLA_Obj X = *(temps[mode]);
            FLA_Set_zero_tensor(X);
            //End X setup

			//print X
			//FLA_Obj_print_matlab("X", X);
            FLA_Psttm(alpha, A, mode, beta, B1, X);
            FLA_Sttsm_single_psttm(alpha, X, mode-1, beta, B, C1, loopCount, temps);
            
//            FLA_Obj_blocked_psym_tensor_free_buffer(&X);
//            FLA_Obj_free_without_buffer(&X);
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
	dim_t blocked_size[order];
	dim_t blocked_stride[order];
	dim_t input_block_size[order];
	dim_t output_block_size[order];
	dim_t temp_block_size[order];
	dim_t flat_size[order];
	
	memcpy(&(blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
	memcpy(&(blocked_stride[0]), &(((A.base)->stride)[0]), order * sizeof(dim_t));
	memcpy(&(input_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
	memcpy(&(output_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(C))[0].size[0]), order * sizeof(dim_t));
	memcpy(&(temp_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
	
	for(i = 0; i < order; i++)
		flat_size[i] = blocked_size[i] * input_block_size[i];
	
	TLA_sym tmpSym = A.sym;
	for(i = order - 1; i > 0; i--){
		TLA_sym Xsym;
		TLA_split_sym_group(tmpSym, 1, &i, &Xsym);

		temp_block_size[i] = output_block_size[i];
		flat_size[i] = output_block_size[i];
		blocked_size[i] = 1;
		for(j = i; j < order - 1; j++)
			blocked_stride[j] = blocked_stride[i];
		temps[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
		FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, temp_block_size, Xsym, temps[i]);
		
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
	dim_t blocked_size[order];
	dim_t blocked_stride[order];
	dim_t input_block_size[order];
	dim_t output_block_size[order];
	dim_t temp_block_size[order];
	dim_t flat_size[order];
	
	memcpy(&(blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
	memcpy(&(blocked_stride[0]), &(((A.base)->stride)[0]), order * sizeof(dim_t));
	memcpy(&(input_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
	memcpy(&(output_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(C))[0].size[0]), order * sizeof(dim_t));
	memcpy(&(temp_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
	
	for(i = 0; i < order; i++)
		flat_size[i] = blocked_size[i] * input_block_size[i];
	
	for(i = order - 1; i > 0; i--){
		temp_block_size[i] = output_block_size[i];
		flat_size[i] /= output_block_size[i];
		blocked_size[i] = 1;
		for(j = i; j < order - 1; j++)
			blocked_stride[j] = blocked_stride[i];
		temps[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
		FLA_Obj_create_blocked_tensor(FLA_DOUBLE, order, flat_size, blocked_stride, input_block_size, temps[i]);
	}
}

void destroy_temporaries(dim_t order, FLA_Obj* temps[]){
	dim_t i;
	for(i = order - 1; i > 0; i--){
		FLA_Obj_blocked_free_buffer(temps[i]);
        FLA_Obj_free_without_buffer(temps[i]);
		FLA_free(temps[i]);
	}
}

//No psym temps

FLA_Error FLA_Sttsm_without_psym_temps( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	FLA_Obj* temps[A.order];
	
	initialize_temporaries(A, C, temps);
    
	FLA_Sttsm_single( alpha, A, FLA_Obj_order(C)-1, beta, B, C, FLA_Obj_dimsize(C,FLA_Obj_order(C)-1)-1, temps);
	
	destroy_temporaries(A.order, temps);
	return FLA_SUCCESS;
}


//Using psym temps

FLA_Error FLA_Sttsm_with_psym_temps( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	FLA_Obj* temps[A.order];

	initialize_psym_temporaries(A, C, temps);
    FLA_Sttsm_single_psttm( alpha, A, FLA_Obj_order(C)-1, beta, B, C, FLA_Obj_dimsize(C,FLA_Obj_order(C)-1)-1, temps);
	

	destroy_psym_temporaries(A.order, temps);
	return FLA_SUCCESS;
}
