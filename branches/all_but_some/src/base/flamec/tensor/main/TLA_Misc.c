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

FLA_Error FLA_Obj_print_scalar_tensor(FLA_Obj A, dim_t repart_mode_index){
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;
	double* buffer;
	
	//View could be under permutation.  Ensure print respects this
	FLA_Part_1xmode2(A, &AT,
						/**/
						&AB, A.permutation[repart_mode_index], 0, FLA_TOP);
	
	while(FLA_Obj_dimsize(AT,A.permutation[repart_mode_index]) < FLA_Obj_dimsize(A,A.permutation[repart_mode_index])){
		FLA_Repart_1xmode2_to_1xmode3(AT,   &A0,
										    &A1,
									  /**/  /**/
									  AB,   &A2, 
									  A.permutation[repart_mode_index], 1, FLA_BOTTOM);
		/************************/
		buffer = FLA_Obj_tensor_buffer_at_view(A1);
		if(repart_mode_index == 0){
			printf("%.3f", *buffer);
			printf(" ");	
		}else{
			FLA_Obj_print_scalar_tensor(A1, repart_mode_index - 1);
		}
		/************************/
		FLA_Cont_with_1xmode3_to_1xmode2(&AT, A0,
											  A1,
										 /**/ /**/
										 &AB, A2, 
										 A.permutation[repart_mode_index], FLA_TOP);
	}
	return FLA_SUCCESS;	
}

//Because we are printing a blocks of tensors
//We alternate between scalar columns and block columns
//This routine prints the element at that location
FLA_Error FLA_Obj_print_scalar_tensor_mode_at(FLA_Obj A, dim_t mode, dim_t index[]){
	FLA_Obj Atmp;
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;
	
	double* buffer;

	dim_t order = FLA_Obj_order(A);
	Atmp = A;
	memcpy(&(Atmp.offset[0]), &(index[0]), order * sizeof(dim_t));
	
	FLA_Part_1xmode2(Atmp, &AT,
					 /**/
					 &AB, mode, 0, FLA_TOP);
	
	while(FLA_Obj_dimsize(AT,mode) < FLA_Obj_dimsize(Atmp,mode)){
		FLA_Repart_1xmode2_to_1xmode3(AT,   &A0,
									  &A1,
									  /**/  /**/
									  AB,   &A2, mode, 1, FLA_BOTTOM);
		/************************/
		buffer = FLA_Obj_tensor_buffer_at_view(A1);
		printf("%.6f", *buffer);
		printf(" ");
		/************************/
		FLA_Cont_with_1xmode3_to_1xmode2(&AT, A0,
										 A1,
										 /**/ /**/
										 &AB, A2, mode, FLA_TOP);
	}
	
	return FLA_SUCCESS;
}

//Loop over mode in scalar object
FLA_Error FLA_Obj_print_hier_tensor_loop_scalar_mode(FLA_Obj A, dim_t mode, dim_t index[]){
    dim_t i;

    dim_t order = A.order;
    dim_t newIndex[FLA_MAX_ORDER];
    FLA_Obj* buffer;
    for(i = 0; i < FLA_Obj_order(A); i++)
        if(FLA_Obj_dimsize(A,i) == 0)
            return FLA_SUCCESS;

    buffer = FLA_Obj_tensor_buffer_at_view(A);
	if(mode == 0){
		const dim_t* permutation = buffer->permutation;
		dim_t ipermutation[FLA_MAX_ORDER];
		dim_t printMode = permutation[0];
		dim_t newIndex[FLA_MAX_ORDER];

		memset(&(ipermutation[0]), 0, order * sizeof(dim_t));
		for(i = 0; i < order; i++){
			ipermutation[permutation[i]] = i;
		}


		for(i = 0; i < order; i++)
			newIndex[i] = index[ipermutation[i]];
			
		
		return FLA_Obj_print_scalar_tensor_mode_at(*buffer, printMode, newIndex);
	}
	
	memcpy(&(newIndex[0]), &(index[0]), order * sizeof(dim_t));
	for(i = 0; i < FLA_Obj_dimsize(*buffer,mode); i++){
		newIndex[mode] = i;
		FLA_Obj_print_hier_tensor_repart_mode_at(A, mode - 1, newIndex);
	}
	
	return FLA_SUCCESS;
}

//Loop over the mode in blocked object
FLA_Error FLA_Obj_print_hier_tensor_repart_mode_at(FLA_Obj A, dim_t repart_mode, dim_t index[]){
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;



	FLA_Part_1xmode2(A, &AT,
						 /**/
						 &AB, repart_mode, 0, FLA_TOP);

	while(FLA_Obj_dimsize(AT,repart_mode) < FLA_Obj_dimsize(A,repart_mode)){
		FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
											  &A1,
										/**/  /**/
										AB,   &A2, repart_mode, 1, FLA_BOTTOM);
		/************************/
		FLA_Obj_print_hier_tensor_loop_scalar_mode(A1, repart_mode, index);
		/************************/
		FLA_Cont_with_1xmode3_to_1xmode2(&AT, A0,
												  A1,
											 /**/ /**/
											 &AB, A2, repart_mode, FLA_TOP);
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Obj_print_tensor(FLA_Obj A){
	FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
	dim_t order = FLA_Obj_order(A);

	if(elemtype == FLA_SCALAR){
		FLA_Obj_print_scalar_tensor(A, order - 1);
	}else{
		dim_t index[FLA_MAX_ORDER];
		memset(&(index[0]), 0, order * sizeof(dim_t));
		FLA_Obj_print_hier_tensor_repart_mode_at(A, order - 1, index);
	}
	
	return FLA_SUCCESS;
}


//Eliminate soon.
FLA_Error FLA_Obj_print(FLA_Obj obj){
    dim_t order = FLA_Obj_order(obj);
    const dim_t* permutation = obj.permutation;
    printf("FLA_Obj\n");
    printf("--------------------------\n");
    printf("  order: %d\n", order);
    print_array("  offset", order, obj.offset);
    print_array("  permutation", order, permutation);
    printf("  nSymGroups: %d\n", obj.sym.nSymGroups);
    print_array("  symGroupLens", obj.sym.nSymGroups, obj.sym.symGroupLens);
    print_array("  symModes", order, obj.sym.symModes);
    printf("--------------------------\n\n");
    return FLA_SUCCESS;
}

//Nice matlab printing
FLA_Error FLA_Obj_print_matlab(const char * varName, FLA_Obj obj){
    dim_t i;
    if(FLA_Obj_order(obj) <= 2)
        printf("%s = reshape([", varName);
    else
        printf("%s = tensor([", varName);
    FLA_Obj_print_tensor(obj);
    printf("],[");
    if(FLA_Obj_elemtype(obj) == FLA_SCALAR){
        for(i = 0; i < FLA_Obj_order(obj); i++)
            printf("%d ", FLA_Obj_dimsize(obj,obj.permutation[i]));
    }else{
        for(i = 0; i < FLA_Obj_order(obj); i++)
            printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(obj)))[0],((FLA_Obj*)(FLA_Obj_base_buffer(obj)))[0].permutation[i]) * FLA_Obj_dimsize(obj,i));
    }
    printf("]);\n\n");
    return FLA_SUCCESS;
}
