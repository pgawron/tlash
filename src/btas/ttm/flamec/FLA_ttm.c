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

//TODO: Handle permuted C
//Repartition until 2-D matrix (Very slow) (Ignore)
FLA_Error FLA_Ttm_single_mode_blis( FLA_Obj alpha, FLA_Obj A,
                                    dim_t mode,
                                    FLA_Obj beta, FLA_Obj B,
                                    FLA_Obj C, dim_t repart_mode )
{

	dim_t* permA = FLA_Obj_permutation(A);

	if((mode == 0 && repart_mode == 1) || repart_mode == 0){
		FLA_Obj A_2D, C_2D;
		dim_t size_A_2D[2];
		dim_t size_C_2D[2];
		dim_t stride_A_2D[2];
		dim_t stride_C_2D[2];
		dim_t firstModeA, firstModeC;

		void* A_buf;
		void* C_buf;

		if(permA[mode] < permA[repart_mode]){
			size_A_2D[0] = FLA_Obj_dimsize(A, permA[mode]);
			size_A_2D[1] = FLA_Obj_dimsize(A, permA[repart_mode]);
			stride_A_2D[0] = FLA_Obj_dimstride(A, permA[mode]);
			stride_A_2D[1] = FLA_Obj_dimstride(A, permA[repart_mode]);
			firstModeA = mode;
		}else{
			size_A_2D[0] = FLA_Obj_dimsize(A, permA[repart_mode]);
			size_A_2D[1] = FLA_Obj_dimsize(A, permA[mode]);
			stride_A_2D[0] = FLA_Obj_dimstride(A, permA[repart_mode]);
			stride_A_2D[1] = FLA_Obj_dimstride(A, permA[mode]);
			firstModeA = repart_mode;
		}
	
		if(mode < repart_mode){
			size_C_2D[0] = FLA_Obj_dimsize(C, mode);
			size_C_2D[1] = FLA_Obj_dimsize(C, repart_mode);
			stride_C_2D[0] = FLA_Obj_dimstride(C, mode);
			stride_C_2D[1] = FLA_Obj_dimstride(C, repart_mode);
			firstModeC = mode;
		}else{
			size_C_2D[0] = FLA_Obj_dimsize(C, repart_mode);
			size_C_2D[1] = FLA_Obj_dimsize(C, mode);
			stride_C_2D[0] = FLA_Obj_dimstride(C, repart_mode);
			stride_C_2D[1] = FLA_Obj_dimstride(C, mode);
			firstModeC = repart_mode;
		}	

		A_buf = FLA_Obj_tensor_buffer_at_view(A);
		C_buf = FLA_Obj_tensor_buffer_at_view(C);
		
		FLA_Obj_create_without_buffer(FLA_DOUBLE, size_A_2D[0], size_A_2D[1], &A_2D);
		FLA_Obj_attach_buffer(A_buf, stride_A_2D[0], stride_A_2D[1], &A_2D);
		
		FLA_Obj_create_without_buffer(FLA_DOUBLE, size_C_2D[0], size_C_2D[1], &C_2D);
		FLA_Obj_attach_buffer(C_buf, stride_C_2D[0], stride_C_2D[1], &C_2D);


		if (mode < repart_mode) {
			if(firstModeC == firstModeA){
				FLA_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, B, A_2D, beta, C_2D);				
			}else{
				FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, B, A_2D, beta, C_2D);				
			}
		}else{
			if(firstModeC == firstModeA){
				FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, A_2D, B, beta, C_2D);				
			}else{
				FLA_Gemm(FLA_TRANSPOSE, FLA_TRANSPOSE, alpha, A_2D, B, beta, C_2D);				
			}
		}
		FLA_Obj_free_without_buffer(&A_2D);
		FLA_Obj_free_without_buffer(&C_2D);
	}else{
	
		//FOR
		FLA_Obj AT, AB;
		FLA_Obj A0, A1, A2;
		FLA_Obj CT, CB;
		FLA_Obj C0, C1, C2;

		FLA_Part_1xmode2(A, &AT,
							&AB, permA[repart_mode], 0, FLA_TOP);	
		FLA_Part_1xmode2(C, &CT,
							&CB, repart_mode, 0, FLA_TOP);	
		//Only symmetric part touched
		//Ponder this
		while(FLA_Obj_dimsize(CT, repart_mode) < FLA_Obj_dimsize(C, repart_mode)){

			//Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
			//Mode-0 of B matches Mode-n of C
			dim_t b = 1;
			FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
										/**/ /**/
											  &A1,
										  AB, &A2, permA[repart_mode], b, FLA_BOTTOM); 
			FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
										/**/ /**/
											  &C1,
										  CB, &C2, repart_mode, b, FLA_BOTTOM);

	/********************/
			if(mode == repart_mode - 1)
				FLA_Ttm_single_mode_blis(alpha, A1, mode, beta, B, C1, repart_mode - 2);
			else
				FLA_Ttm_single_mode_blis(alpha, A1, mode, beta, B, C1, repart_mode - 1);
	/********************/

			FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
												   A1,
											/********/
											  &AB, A2, permA[repart_mode], FLA_TOP);
			FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
												   C1,
											/********/
											  &CB, C2, repart_mode, FLA_TOP);
		}
	}
	return FLA_SUCCESS;
}

/*****************/
 // End new code
/*****************/

//Original ttm implementation for no_permC variant (permutes A but assumes C is in permuted form)
FLA_Error FLA_Ttm_scalar_no_permC( FLA_Obj alpha, FLA_Obj A,
                                   dim_t mode,
                                   FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	FLA_Datatype datatype = FLA_Obj_datatype(A);
	FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
	dim_t order = FLA_Obj_order(A);
    dim_t permutation[FLA_MAX_ORDER];
    dim_t i;
	const dim_t* size_A = A.size;
	const dim_t* stride_A = (A.base)->stride;

	FLA_Obj P;

	if(elemtype != FLA_SCALAR){
		printf("NON-scalar detected in final ttm\n");
		return FLA_SUCCESS;
	}
	
	permutation[0] = mode;
    for(i = 0; i < mode; i++)
		permutation[i+1] = i;
    for(i = mode+1; i < order; i++)
		permutation[i] = i;

	
	FLA_Obj_create_tensor(datatype, order, size_A, stride_A, &P);
	FLA_Permute(A, permutation, &P);

	FLA_Adjust_2D_info(&B);
	FLA_Adjust_2D_info(&C);
	FLA_Adjust_2D_info(&P);

	/*********************************/
	FLASH_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, beta, B, P, alpha, C);
	/*********************************/

	FLA_Obj_free_buffer(&P);
	FLA_Obj_free_without_buffer(&P);

	return FLA_SUCCESS;
}

//Original form that permutes C
FLA_Error FLA_Ttm_scalar_permC( FLA_Obj alpha, FLA_Obj A,
                                dim_t mode,
                                FLA_Obj beta, FLA_Obj B,
                                FLA_Obj C )
{
    FLA_Datatype datatype = FLA_Obj_datatype(A);
    FLA_Elemtype elemtype = FLA_Obj_elemtype(A);
    dim_t order = FLA_Obj_order(A);
    dim_t ipermutation[FLA_MAX_ORDER];
    dim_t permutation[FLA_MAX_ORDER];
    dim_t i;
    const dim_t* size_A = A.size;
    const dim_t* stride_A = (A.base)->stride;
    const dim_t* size_C = C.size;
    const dim_t* stride_C = (C.base)->stride;

    FLA_Obj P, tmpC;

    if(elemtype != FLA_SCALAR){
        printf("NON-scalar detected in final ttm\n");
        return FLA_SUCCESS;
    }

    permutation[0] = mode;
    for(i = 0; i < mode; i++)
        permutation[i+1] = i;
    for(i = mode+1; i < order; i++)
        permutation[i] = i;

  FLA_Obj_create_tensor(datatype, order, size_A, stride_A, &P);
  FLA_Obj_create_tensor(datatype, order, size_C, stride_C, &tmpC);

  FLA_Permute(A, permutation, &P);
  FLA_Permute(C, permutation, &tmpC);

  FLA_Adjust_2D_info(&P);
  FLA_Adjust_2D_info(&tmpC);
	FLA_Adjust_2D_info(&B);
  FLASH_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, beta, B, P, alpha, tmpC);

  for(i = 0; i < order; i++)
      ipermutation[permutation[i]] = i;

  FLA_Permute(tmpC, ipermutation, &C);

  FLA_Obj_free_buffer(&tmpC);
  FLA_Obj_free_without_buffer(&tmpC);
  FLA_Obj_free_buffer(&P);
  FLA_Obj_free_without_buffer(&P);

    return FLA_SUCCESS;
}
/*
FLA_Error FLA_Ttm( FLA_Obj alpha, FLA_Obj A,
                   dim_t nModes, dim_t mode[nModes],
                   FLA_Obj beta, FLA_Obj B[nModes],
                   FLA_Obj C )
{
	dim_t i, j;
	FLA_Obj tmpA, tmpC;
	FLA_Datatype datatype = FLA_Obj_datatype( A );
	dim_t order;
	dim_t* size_tmpC;
	dim_t* size_tmpA;
	dim_t* stride_tmpC;
	dim_t* stride_tmpA;
	//When we mode-mult, we get B's second dim size
	dim_t  newdimsize[nModes];
	for(i = 0; i < nModes; i++)
		newdimsize[i] =  FLA_Obj_dimsize(B[i], 0);

	//Setup tmpA
	order = FLA_Obj_order(A);
	size_tmpA = FLA_Obj_size(A);
	stride_tmpA = FLA_Obj_stride(A);

	//WARNING: DOUBLE ATTACHMENT.  FIX THIS
	FLA_Obj_create_tensor_without_buffer(datatype, order, size_tmpA, &tmpA);
	double* tmpA_buf = (double*)FLA_malloc(FLA_Obj_num_elem_alloc(A) * sizeof(double*));
	memcpy(&(tmpA_buf[0]), &(((double*)FLA_Obj_base_buffer(A))[0]), FLA_Obj_num_elem_alloc(A) * sizeof(double));
    FLA_Obj_attach_buffer_to_tensor(tmpA_buf, order, stride_tmpA, &tmpA);

	for(i = 0; i < nModes; i++){
		//Setup tmpC
		size_tmpC = FLA_Obj_size(tmpA);
		size_tmpC[mode[i]] = newdimsize[i];
		stride_tmpC = FLA_Obj_stride(tmpA);
		dim_t nElem = size_tmpC[0];
		stride_tmpC[0] = 1;
		for(j = 1; j < order; j++){
			nElem *= size_tmpC[j];
			stride_tmpC[j] = stride_tmpC[j-1]*size_tmpC[j-1];
		}

		double* data = (double*)FLA_malloc(nElem * sizeof(double));
		memset(&(data[0]), 0, nElem * sizeof(double));

		FLA_Obj_create_tensor_without_buffer(datatype, order, size_tmpC, &tmpC);
	    FLA_Obj_attach_buffer_to_tensor(data, order, stride_tmpC, &tmpC);

		FLA_Ttm_single_mode(alpha, tmpA, mode[i], beta, B[i], tmpC);

		//Make tmpA = tmpC and clear tmpC for next iteration
		memcpy(&(tmpA.size[0]), &(tmpC.size[0]), order * sizeof(dim_t));
		memcpy(&((tmpA.base)->size[0]), &(tmpC.size[0]), order * sizeof(dim_t));
		memcpy(&(((tmpA.base)->stride)[0]), &(((tmpC.base)->stride)[0]), order * sizeof(dim_t));
		//Clear tmpA's buffer since we will point it elsewhere
		FLA_free((tmpA.base)->buffer);

		tmpA.base->buffer = FLA_malloc(tmpC.base->n_elem_alloc * sizeof(double));
		memcpy(&(((double*)tmpA.base->buffer)[0]), &(((double*)tmpC.base->buffer)[0]), tmpC.base->n_elem_alloc * sizeof(double));
		tmpA.base->n_elem_alloc = tmpC.base->n_elem_alloc;
		FLA_Adjust_2D_info(&tmpA);

		FLA_Obj_free(&tmpC);
		FLA_free(size_tmpC);
		FLA_free(stride_tmpC);
	}

	memcpy(&(((double*)C.base->buffer)[0]), &(((double*)tmpA.base->buffer)[0]), FLA_Obj_num_elem_alloc(C) * sizeof(double));

	FLA_Obj_free(&tmpA);
	FLA_free(size_tmpA);
	FLA_free(stride_tmpA);
	return FLA_SUCCESS;
}
*/

/*****************
*****    NOTE: ASSUMES THIS IS COMPUTING STATIONARY C ALGORITHM!!!!!
*****************/
/*
FLA_Error FLA_Ttm_hierAB_single_repart_mode( FLA_Obj alpha, FLA_Obj A,
                                             dim_t mode,
                                             FLA_Obj beta, FLA_Obj B,
                                             dim_t repart_mode, FLA_Obj C )
{
	//FOR
	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;
	
	FLA_Part_1xmode2(B, &BT,
						&BB, 1, 0, FLA_TOP);	
	FLA_Part_1xmode2(A, &AT,
						&AB, repart_mode, 0, FLA_TOP);	


	//Stationary-C specific code
	dim_t i;
	FLA_Datatype datatype = FLA_Obj_datatype( C );
	dim_t order = FLA_Obj_order(C);
	dim_t permutation[order];
	dim_t ipermutation[order];
	permutation[0] = mode;
	dim_t* size_C;
	dim_t* stride_C;
		
	size_C = FLA_Obj_size(C);
	stride_C = FLA_Obj_stride(C);

	for(i = 0; i < mode; i++)
		permutation[i+1] = i;
	for(i = mode+1; i < order; i++)
		permutation[i] = i;

	FLA_Obj tmpC;

	FLA_Obj_create_tensor(datatype, order, size_C, stride_C, &tmpC);
	FLA_Permute(C, permutation, &tmpC);
	//End Stationry-C code

	//Only symmetric part touched
	//Ponder this
	dim_t loopCount = 0;
	while(loopCount < FLA_Obj_dimsize(A, repart_mode)){

		//Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
		//Mode-0 of B matches Mode-n of C
		dim_t b = 1;
		FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
									//// ////
										  &B1,
									  BB, &B2, 1, b, FLA_BOTTOM); 
		FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
									//// ////
										  &A1,
									  AB, &A2, repart_mode, b, FLA_BOTTOM);

		FLA_Ttm_single_mode(alpha, A1, mode, beta, B1, tmpC);
		//Stationary C calls nopermC I think...


		FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
											   A1,
										//////////
										  &AB, A2, repart_mode, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
											   B1,
										//////////
										  &BB, B2, 1, FLA_TOP);
		loopCount++;
	}

	//Stationary-C code
	for(i = 0; i < order; i++)
		ipermutation[permutation[i]] = i;

	FLA_Permute(tmpC, ipermutation, &C);

	FLA_Obj_free_buffer(&tmpC);
	FLA_Obj_free_without_buffer(&tmpC);
	FLA_free(size_C);
	FLA_free(stride_C);
	//End Stationary-C code

	return FLA_SUCCESS;
}


FLA_Error FLA_Ttm_hierCB_single_repart_mode( FLA_Obj alpha, FLA_Obj A,
                                             dim_t mode,
                                             FLA_Obj beta, FLA_Obj B,
                                             dim_t repart_mode, FLA_Obj C )
{
	//FOR
	FLA_Obj BT, BB;
	FLA_Obj B0, B1, B2;
	FLA_Obj CT, CB;
	FLA_Obj C0, C1, C2;
	
	FLA_Part_1xmode2(B, &BT,
						&BB, 0, 0, FLA_TOP);	
	FLA_Part_1xmode2(C, &CT,
						&CB, repart_mode, 0, FLA_TOP);	
	//Only symmetric part touched
	//Ponder this
	dim_t loopCount = 0;
	while(loopCount < FLA_Obj_dimsize(C, repart_mode)){

		//Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
		//Mode-0 of B matches Mode-n of C
		dim_t b = 1;
		FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
									//// ////
										  &B1,
									  BB, &B2, 1, b, FLA_BOTTOM); 
		FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
									//// ////
										  &C1,
									  CB, &C2, repart_mode, b, FLA_BOTTOM);


		FLA_Ttm_single_mode(alpha, A, mode, beta, B1, C1);

		FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
											   C1,
										//////////
										  &CB, C2, repart_mode, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
											   B1,
										//////////
										  &BB, B2, 0, FLA_TOP);
		loopCount++;
	}

	return FLA_SUCCESS;
}
*/

FLA_Error FLA_Ttm_hierCA_single_repart_mode( FLA_Obj alpha, FLA_Obj A,
                                             dim_t mode,
                                             FLA_Obj beta, FLA_Obj B,
                                             dim_t repart_mode, FLA_Obj C )
{
	//FOR
	FLA_Obj AT, AB;
	FLA_Obj A0, A1, A2;
	FLA_Obj CT, CB;
	FLA_Obj C0, C1, C2;
	
	dim_t loopCount;

	FLA_Part_1xmode2(A, &AT,
						&AB, repart_mode, 0, FLA_TOP);	
	FLA_Part_1xmode2(C, &CT,
						&CB, repart_mode, 0, FLA_TOP);	
	//Only symmetric part touched
	//Ponder this
	loopCount = 0;
	while(loopCount < FLA_Obj_dimsize(C, repart_mode)){
		dim_t b = 1;
		FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
									/**/ /**/
										  &A1,
									  AB, &A2, repart_mode, b, FLA_BOTTOM); 
		FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
									/**/ /**/
										  &C1,
									  CB, &C2, repart_mode, b, FLA_BOTTOM);

		FLA_Ttm_single_mode(alpha, A1, mode, beta, B, C1);

		FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
											   C1,
										/********/
										  &CB, C2, repart_mode, FLA_TOP);
		FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
											   A1,
										/********/
										  &AB, A2, repart_mode, FLA_TOP);
		loopCount++;
	}

	return FLA_SUCCESS;
}

//Assumes C is already permuted (hence nopermC)
FLA_Error FLA_Tensor_innerprod_nopermC( FLA_Obj alpha, FLA_Obj A,
                                        dim_t mode,
                                        FLA_Obj beta, FLA_Obj B,
                                        FLA_Obj C )
{
    FLA_Obj BT, BB;
    FLA_Obj B0, B1, B2;
    FLA_Obj AT, AB;
    FLA_Obj A0, A1, A2;

    dim_t loopCount;

    FLA_Obj A1blk;
    FLA_Obj B1blk;

    FLA_Part_1xmode2(B, &BT,
                        &BB, 1, 0, FLA_TOP);
    FLA_Part_1xmode2(A, &AT,
                        &AB, mode, 0, FLA_TOP);
    //Only symmetric part touched
    //Ponder this
    loopCount = 0;
    while(loopCount < FLA_Obj_dimsize(A, mode)){

        //Check this mathc out.  I think it is correct, Mode-1 of B matches mode-n of A
        //Mode-0 of B matches Mode-n of C
        dim_t b = 1;
        FLA_Repart_1xmode2_to_1xmode3(BT, &B0,
                                    /**/ /**/
                                          &B1,
                                      BB, &B2, 1, b, FLA_BOTTOM);
        FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
                                    /**/ /**/
                                          &A1,
                                      AB, &A2, mode, b, FLA_BOTTOM);
        /*********************/
        A1blk = *((FLA_Obj*)FLA_Obj_tensor_buffer_at_view(A1));
        B1blk = *((FLA_Obj*)FLA_Obj_tensor_buffer_at_view(B1));
        FLA_Ttm_scalar_no_permC(alpha, A1blk, mode, beta, B1blk, C);
		/*********************/
        FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
                                               A1,
                                        /********/
                                          &AB, A2, mode, FLA_TOP);
        FLA_Cont_with_1xmode3_to_1xmode2( &BT, B0,
                                               B1,
                                        /********/
                                          &BB, B2, 1, FLA_TOP);
        loopCount++;
    }
    return FLA_SUCCESS;
}

//C is a vector, B a matrix, A a vector.
FLA_Error FLA_Tensor_mvmult_nopermC( FLA_Obj alpha, FLA_Obj A,
                                     dim_t mode,
                                     FLA_Obj beta, FLA_Obj B,
                                     FLA_Obj C )
{
    FLA_Obj BT, BB;
    FLA_Obj B0, B1, B2;
    FLA_Obj CT, CB;
    FLA_Obj C0, C1, C2;

    dim_t loopCount;

    //Normally each block of output vector would be re-permuted for each contribution.
    //To avoid this, we permute the output block once and compute assuming the data is
    //already permuted, then permute back.  This saves lots of mem-ops
	FLA_Obj tmpC1;

	
    FLA_Part_1xmode2(B, &BT,
                        &BB, 0, 0, FLA_TOP);
    FLA_Part_1xmode2(C, &CT,
                        &CB, mode, 0, FLA_TOP);
    //Only symmetric part touched
    //Ponder this
    loopCount = 0;
    while(loopCount < FLA_Obj_dimsize(C, mode)){

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

        //NOTE: Curly braces for C89 compliance.  Find a better way to fix this
        {
			//Permute this block of C
			FLA_Obj* C1blk = (FLA_Obj*)FLA_Obj_tensor_buffer_at_view(C1);
			dim_t i;
			dim_t order = FLA_Obj_order(*C1blk);
			const dim_t* size_C1blk = C1blk->size;
			const dim_t* stride_C1blk = C1blk->base->stride;
			dim_t permutation[FLA_MAX_ORDER];
			dim_t ipermutation[FLA_MAX_ORDER];
			permutation[0] = mode;
			for(i = 0; i < mode; i++)
				permutation[i+1] = i;
			for(i = mode+1; i < order; i++)
				permutation[i] = i;
			if(loopCount == 0)
				FLA_Obj_create_tensor(FLA_DOUBLE, C1blk->order, size_C1blk, stride_C1blk, &tmpC1);
			FLA_Permute(*C1blk, permutation, &tmpC1);
			//Done permuting

			/***********************************************/
			//Compute
			FLA_Tensor_innerprod_nopermC(alpha, A, mode, beta, B1, tmpC1);
			/***********************************************/

			//Permute the block back
			for(i = 0; i < order; i++)
				ipermutation[permutation[i]] = i;
			FLA_Permute(tmpC1, ipermutation, C1blk);
			//Done permuting
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
	
	FLA_Obj_free_buffer(&tmpC1);
	FLA_Obj_free_without_buffer(&tmpC1);
    return FLA_SUCCESS;
}

FLA_Error FLA_Ttm_single_mode( FLA_Obj alpha, FLA_Obj A,
                               dim_t mode,
                               FLA_Obj beta, FLA_Obj B,
                               FLA_Obj C )
{
	dim_t i;
	dim_t do_repart;
	dim_t repart_mode;

	//If A is a scalar, call the scalar version
	if(FLA_Obj_elemtype(A) == FLA_SCALAR){
		return FLA_Ttm_scalar_permC(alpha, A, mode, beta, B, C);
	}
	//Repartition C & A as much as you can before multiplying
	do_repart = FALSE;
	repart_mode = 0;
	for(i = 0; i < FLA_Obj_order(C); i++){
		if(i == mode)
			continue;
		else if(FLA_Obj_dimsize(C,i) > 1){
			repart_mode = i;
			do_repart = TRUE;
			break;
		}
	}

	//Multiply if can't repartition anymore
    if(do_repart == FALSE)
        return FLA_Tensor_mvmult_nopermC(alpha, A, mode, beta, B, C);
    else
        return FLA_Ttm_hierCA_single_repart_mode(alpha, A, mode, beta, B, repart_mode, C);
}

//Single ttm without permuting C
FLA_Error FLA_Ttm_single_mode_no_permC( FLA_Obj alpha, FLA_Obj A,
                                        dim_t mode,
                                        FLA_Obj beta, FLA_Obj B,
                                        FLA_Obj C )
{
	FLA_Elemtype elemtype_A = FLA_Obj_elemtype(A);
	FLA_Elemtype elemtype_B = FLA_Obj_elemtype(B);
	FLA_Elemtype elemtype_C = FLA_Obj_elemtype(C);

	if(elemtype_A == FLA_SCALAR && elemtype_B == FLA_SCALAR && elemtype_C == FLA_SCALAR){
		FLA_Ttm_scalar_no_permC(alpha, A, mode, beta, B, C);
	}else{
		dim_t i;
		dim_t singleElemA = TRUE;
		dim_t singleElemB = TRUE;
		dim_t singleElemC = TRUE;
		void* buf_A = FLA_Obj_base_buffer(A);
		void* buf_B = FLA_Obj_base_buffer(B);
		void* buf_C = FLA_Obj_base_buffer(C);

		for(i = 0; i < FLA_Obj_order(A); i++)
			if(FLA_Obj_dimsize(A, i) != 1){
				singleElemA = FALSE;
				break;
			}
		for(i = 0; i < FLA_Obj_order(B); i++)
			if(FLA_Obj_dimsize(B, i) != 1){
				singleElemB = FALSE;
				break;
			}
		for(i = 0; i < FLA_Obj_order(C); i++)
			if(FLA_Obj_dimsize(C, i) != 1){
				singleElemC = FALSE;
				break;
			}

		if(singleElemC && elemtype_C != FLA_SCALAR){
			dim_t linIndex = 0;
			dim_t order = FLA_Obj_order(C);
			linIndex = FLA_TIndex_to_LinIndex(order, &(((C.base)->stride)[0]), &(C.offset[0]));
			FLA_Ttm_single_mode_no_permC(alpha, A, mode, beta, B, ((FLA_Obj*)buf_C)[linIndex]);
		}else if(singleElemA && elemtype_A != FLA_SCALAR){
			dim_t linIndex = 0;
			dim_t order = FLA_Obj_order(A);
			linIndex = FLA_TIndex_to_LinIndex(order, &(((A.base)->stride)[0]), &(A.offset[0]));
			FLA_Ttm_single_mode_no_permC(alpha, ((FLA_Obj*)buf_A)[linIndex], mode, beta, B, C);
		}else if(singleElemB && elemtype_B != FLA_SCALAR){
			dim_t linIndex = 0;
			dim_t order = FLA_Obj_order(B);
			linIndex = FLA_TIndex_to_LinIndex(order, &(((B.base)->stride)[0]), &(B.offset[0]));
			FLA_Ttm_single_mode_no_permC(alpha, A, mode, beta, ((FLA_Obj*)buf_B)[linIndex], C);
		}
	}
	return FLA_SUCCESS;
}
