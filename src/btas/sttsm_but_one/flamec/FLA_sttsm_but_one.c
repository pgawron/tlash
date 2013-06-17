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

FLA_Error FLA_Sttsm_but_one_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, dim_t ignore_mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex, FLA_Obj* temps[] )
{

    if(mode == ignore_mode){
        if (mode != 0)
            FLA_Sttsm_but_one_single(alpha, A, mode-1, ignore_mode, beta, B, C, endIndex, temps);
        return FLA_SUCCESS;
    }
    
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

		if((mode == 0 && ignore_mode != 0) || (mode == 1 && ignore_mode == 0)){
			FLA_Ttm_single_mode(alpha, A, mode, beta, B1, C1);
		}else{
			FLA_Obj X = *(temps[mode]);
            FLA_Set_zero_tensor(*(temps[mode]));
            
			FLA_Psttm(alpha, A, mode, beta, B1, X);
            FLA_Sttsm_but_one_single(alpha, X, mode-1, ignore_mode, beta, B, C1, loopCount, temps);
		}
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

void initialize_psym_but_one_temporaries(FLA_Obj A, FLA_Obj C, dim_t ignore_mode, FLA_Obj* temps[]){
	dim_t i, j;
	dim_t order = FLA_Obj_order(A);
	dim_t temp_blocked_size[order];
	dim_t temp_blocked_stride[order];
	dim_t input_block_size[order];
	dim_t output_block_size[order];
	dim_t temp_block_size[order];
	dim_t temp_flat_size[order];
	


	memcpy(&(input_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
	memcpy(&(output_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(C))[0].size[0]), order * sizeof(dim_t));

	TLA_sym tmpSym = A.sym;
	for(i = order - 1; i > 0; i--){
		if(i != ignore_mode){
			TLA_sym Xsym;
			TLA_split_sym_group(tmpSym, 1, &i, &Xsym);

			//Set up the size and stride of the temporaries
			memcpy(&(temp_blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
			memcpy(&(temp_block_size[0]), &(((FLA_Obj*)FLA_Obj_base_buffer(A))[0].size[0]), order * sizeof(dim_t));
			for(j = i; j < order; j++){
				if(j != ignore_mode){
					temp_blocked_size[j] = 1;
					temp_block_size[j] = ((FLA_Obj*)FLA_Obj_base_buffer(C))[0].size[j];
				}
			}

			for(j = 0; j < order; j++)
				temp_flat_size[j] = temp_blocked_size[j] * temp_block_size[j];

			temp_blocked_stride[0] = 1;
			for(j = 1; j < order; j++)
				temp_blocked_stride[j] = temp_blocked_stride[j-1] * temp_blocked_size[j-1];
            
			temps[i] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
			FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, temp_flat_size, temp_blocked_stride, temp_block_size, Xsym, temps[i]);
			
			tmpSym = Xsym;
		}
	}
}

void destroy_psym_but_one_temporaries(dim_t order, dim_t ignore_mode, FLA_Obj* temps[]){
	dim_t i;
	for(i = order - 1; i > 0; i--){
		if (i != ignore_mode){
			FLA_Obj_blocked_psym_tensor_free_buffer(temps[i]);
			FLA_Obj_free_without_buffer(temps[i]);
		}
	}
}
//Using psym temps

FLA_Error FLA_Sttsm_but_one( FLA_Obj alpha, FLA_Obj A, dim_t ignore_mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	FLA_Obj* temps[A.order];

	initialize_psym_but_one_temporaries(A, C, ignore_mode, temps);
    if(ignore_mode == A.order - 1)
        FLA_Sttsm_but_one_single( alpha, A, C.order-2, ignore_mode, beta, B, C, FLA_Obj_dimsize(C,C.order-2)-1, temps);
    else
        FLA_Sttsm_but_one_single( alpha, A, C.order-1, ignore_mode, beta, B, C, FLA_Obj_dimsize(C,C.order-1)-1, temps);

	destroy_psym_but_one_temporaries(A.order, ignore_mode, temps);
	return FLA_SUCCESS;
}
