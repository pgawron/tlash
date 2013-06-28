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
FLA_Error FLA_Sttsm_but_some_compute( FLA_Obj alpha, FLA_Obj A, dim_t compute_ind, dim_t compute_modes[], dim_t ignore_ind, dim_t ignore_modes[], FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex_repart, dim_t endIndex_compute, FLA_Obj* temps[] )
{
	//FOR
	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;
	FLA_Obj CT, CB;
	FLA_Obj C0, C1, C2;
	dim_t mode = compute_modes[compute_ind];

	FLA_Part_1xmode2(B, &BT,
						&BB, 0, 0, FLA_TOP);
	FLA_Part_1xmode2(C, &CT,
						&CB, mode, 0, FLA_TOP);
	//Only symmetric part touched
	//Ponder this
	dim_t loopCount = 0;
	while(loopCount <= endIndex_compute){
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

		//Make sure that if we are at the bottom of the recursion,
		//We multiply into C (not X)
		if(compute_ind == 0){
            printf("Single mode compute\n");
            print_array("Aoffset", A.order, A.offset);
            FLA_Obj_print_matlab("A", A);
            print_array("B1offset", B1.order, B1.offset);
            FLA_Obj_print_matlab("B1", B1);
            printf("mode: %d\n", mode);
            print_array("preC1offset", C1.order, C1.offset);
            FLA_Obj_print_matlab("preC1", C1);
			FLA_Psttm(alpha, A, mode, beta, B1, C1);
			print_array("C1offset", C1.order, C1.offset);
            FLA_Obj_print_matlab("C1", C1);
		}else{
			//Set the temporary to zero
			FLA_Obj X = *(temps[mode]);
            FLA_Set_zero_tensor(*(temps[mode]));
            
            //Compute temporary
            FLA_Obj_print_matlab("A", A);
            FLA_Obj_print_matlab("B1", B1);
            printf("mode: %d\n", mode);
            FLA_Obj_print_matlab("preX", X);
			FLA_Psttm(alpha, A, mode, beta, B1, X);
            FLA_Obj_print_matlab("X", X);
            
			//Use temporary for recursion
            //Figure out how to fix guard
            FLA_Sttsm_but_some_compute(alpha, X, compute_ind - 1, compute_modes, ignore_ind, ignore_modes, beta, B, C1, endIndex_repart, loopCount, temps);
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

FLA_Error FLA_Sttsm_but_some_repart( FLA_Obj alpha, FLA_Obj A, dim_t compute_ind, dim_t compute_modes[], dim_t ignore_ind, dim_t ignore_modes[], FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t endIndex_repart, dim_t endIndex_compute, FLA_Obj* temps[] )
{
	//FOR
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;
	FLA_Obj CT, CB;
	FLA_Obj C0, C1, C2;
	dim_t mode = ignore_modes[ignore_ind];

	FLA_Part_1xmode2(A, &AT,
						&AB, mode, 0, FLA_TOP);
	FLA_Part_1xmode2(C, &CT,
						&CB, mode, 0, FLA_TOP);
	//Only symmetric part touched
	//Ponder this
	dim_t loopCount = 0;
	while(loopCount <= endIndex_repart){
		//Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
		//Mode-0 of B matches Mode-n of C
		dim_t b = 1;
		FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
									/**/ /**/
										  &A1,
									  AB, &A2, mode, b, FLA_BOTTOM);
		FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
									/**/ /**/
										  &C1,
									  CB, &C2, mode, b, FLA_BOTTOM);

		//Figure out how to fix guard
		if(ignore_ind == 0){
			FLA_Sttsm_but_some_compute(alpha, A1, compute_ind, compute_modes, ignore_ind, ignore_modes, beta, B, C1, loopCount, endIndex_compute, temps);
		}else{
			FLA_Sttsm_but_some_repart(alpha, A1, compute_ind, compute_modes, ignore_ind - 1, ignore_modes, beta, B, C1, loopCount, endIndex_compute, temps);
        }

		FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
											   C1,
										/********/
										  &CB, C2, mode, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
											   A1,
										/********/
										  &AB, A2, mode, FLA_TOP);
		loopCount++;
	}

	return FLA_SUCCESS;
}



void initialize_psym_but_some_temporaries(FLA_Obj A, FLA_Obj C, dim_t nModes_compute, dim_t compute_modes[], dim_t nModes_ignore, dim_t ignore_modes[], FLA_Obj* temps[]){
	dim_t i, j;
	dim_t order = FLA_Obj_order(A);
	dim_t temp_blocked_size[FLA_MAX_ORDER];
	dim_t temp_blocked_stride[FLA_MAX_ORDER];
	dim_t temp_block_size[FLA_MAX_ORDER];
	dim_t temp_flat_size[FLA_MAX_ORDER];
	FLA_Obj* buf_A = (FLA_Obj*)FLA_Obj_base_buffer(A);
	FLA_Obj* buf_C = (FLA_Obj*)FLA_Obj_base_buffer(C);
	TLA_sym tmpSym;
    TLA_sym A_post_part_sym = A.sym;
    
    for(i = 0; i < nModes_ignore; i++)
        TLA_split_sym_group(A_post_part_sym, 1, &(ignore_modes[i]), &A_post_part_sym);

	for(i = nModes_compute - 1; i <= nModes_compute; i--){
		TLA_sym Xsym;
		tmpSym = A_post_part_sym;
		for(j = nModes_compute - 1; j >= i && j < nModes_compute; j--)
			TLA_split_sym_group(tmpSym, 1, &(compute_modes[j]), &tmpSym);
		Xsym = tmpSym;

		//Set up the size and stride of the temporaries
		memcpy(&(temp_blocked_size[0]), &(A.size[0]), order * sizeof(dim_t));
		memcpy(&(temp_block_size[0]), &(buf_A[0].size[0]), order * sizeof(dim_t));

		//NOTE: Assume block partition size = 1
		for(j = i; j < nModes_compute; j++){
			temp_blocked_size[compute_modes[j]] = 1;
			temp_block_size[compute_modes[j]] = buf_C[0].size[compute_modes[j]];
		}

		for(j = 0; j < nModes_ignore; j++){
			temp_blocked_size[ignore_modes[j]] = 1;
			temp_block_size[ignore_modes[j]] = buf_A[0].size[ignore_modes[j]];
		}

		//Adjust size & stride
		FLA_array_elemwise_product(order, temp_blocked_size, temp_block_size, temp_flat_size);
		FLA_Set_tensor_stride(order, temp_blocked_size, temp_blocked_stride);

		//Create the temporary
		temps[compute_modes[i]] = (FLA_Obj*)FLA_malloc(sizeof(FLA_Obj));
		FLA_Obj_create_blocked_psym_tensor(FLA_DOUBLE, order, temp_flat_size, temp_blocked_stride, temp_block_size, Xsym, temps[compute_modes[i]]);
	}
}

void destroy_psym_but_some_temporaries(dim_t nModes_compute, dim_t compute_modes[], FLA_Obj* temps[]){
	dim_t i;

	for(i = 0; i < nModes_compute; i++){
		//Skip the modes we didn't do anything with
		FLA_Obj_blocked_psym_tensor_free_buffer(temps[compute_modes[i]]);
		FLA_Obj_free_without_buffer(temps[compute_modes[i]]);
		FLA_free(temps[compute_modes[i]]);
	}
}

//Using psym temps

FLA_Error FLA_Sttsm_but_some( FLA_Obj alpha, FLA_Obj A, dim_t nModes_ignore, dim_t ignore_modes[nModes_ignore], FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	dim_t i;
	FLA_Obj* temps[FLA_MAX_ORDER];
	dim_t compute_modes[FLA_MAX_ORDER];


    dim_t ignore_ind = 0;
    dim_t compute_ind = 0;
	for(i = 0; i < A.order; i++){
		if(i != ignore_modes[ignore_ind]){
			compute_modes[compute_ind] = i;
			compute_ind++;
		}else{
			ignore_ind++;
		}
	}

	initialize_psym_but_some_temporaries(A, C, compute_ind, compute_modes, nModes_ignore, ignore_modes, temps);

	FLA_Sttsm_but_some_repart(alpha, A, compute_ind - 1, compute_modes, nModes_ignore - 1, ignore_modes, beta, B, C, FLA_Obj_dimsize(C, ignore_modes[0]) - 1, FLA_Obj_dimsize(C, compute_modes[0]) - 1, temps);

    //Cleanup
    destroy_psym_but_some_temporaries(compute_ind, compute_modes, temps);

	return FLA_SUCCESS;
}
