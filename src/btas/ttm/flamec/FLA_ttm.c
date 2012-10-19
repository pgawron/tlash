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

FLA_Error FLA_Ttm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
//	printf("res in:\n");
//	FLA_Obj_print_tensor(C);

	dim_t order = A.order;
    dim_t permutation[order];
	dim_t ipermutation[order];
    dim_t i;
	permutation[0] = mode;
    for(i = 0; i < mode; i++)
		permutation[i+1] = i;
    for(i = mode+1; i < order; i++)
		permutation[i] = i;

	FLA_Obj P, tmpC;

	FLA_Permute_hier(A, permutation, &P);
	FLA_Permute_hier(C, permutation, &tmpC);

	printf("A precomp:\n");
	FLA_Obj_print_tensor(A);
	
	printf("B precomp:\n");
	FLA_Obj_print_tensor(B);

	printf("P precomp:\n");
	FLA_Obj_print_tensor(P);

	printf("tmp precomp:\n");
	FLA_Obj_print_tensor(tmpC);

	//Maybe casting as flash works?
	FLASH_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, beta, B, P, alpha, tmpC);

	printf("tmp postcomp:\n");
	FLA_Obj_print_tensor(tmpC);
	
	for(i = 0; i < order; i++)
		ipermutation[permutation[i]] = i;

	FLA_Permute_hier(tmpC, ipermutation, &C);
/*
	printf("res out:\n");
	FLA_Obj_print_tensor(C);
*/
	return FLA_SUCCESS;
}

FLA_Error FLA_Ttm( FLA_Obj alpha, FLA_Obj A, dim_t nModes, dim_t mode[nModes], FLA_Obj beta, FLA_Obj B[nModes], FLA_Obj C )
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

	FLA_Obj_create_tensor_without_buffer(datatype, order, size_tmpA, &tmpA);
    FLA_Obj_attach_buffer_to_tensor((A.base)->buffer, order, stride_tmpA, &tmpA);

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

		//All objects setup, now perform computation
		//Probably should remove alpha + beta after first compute
		FLA_Ttm_single(alpha, tmpA, mode[i], beta, B[i], tmpC);

		//Make tmpA = tmpC and clear tmpC for next iteration
		memcpy(&(tmpA.size[0]), &(tmpC.size[0]), order * sizeof(dim_t));
		memcpy(&(((tmpA.base)->stride)[0]), &(((tmpC.base)->stride)[0]), order * sizeof(dim_t));
		(tmpA.base)->buffer = (tmpC.base)->buffer;
		FLA_Adjust_2D_info(&tmpA);
	}

	(C.base)->buffer = (tmpC.base)->buffer;
	return FLA_SUCCESS;
}
