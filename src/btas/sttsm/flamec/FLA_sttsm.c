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
	if(mode == 0){
		FLA_Ttm_single(alpha, A, mode, beta, B, C);
	}else{
		dim_t i;
		FLA_Obj X;
		dim_t order = FLA_Obj_order( A );
		dim_t* size_X = FLA_Obj_size( A );
		size_X[mode] = FLA_Obj_dimsize( B, 0);

		dim_t nElem = size_X[0];
		dim_t stride_X[order];
		stride_X[0] = 1;
		for(i = 1; i < order; i++){
			nElem *= size_X[i];
			stride_X[i] = stride_X[i-1]*size_X[i-1];
		}
		FLA_Obj_create_tensor_without_buffer(FLA_DOUBLE, order, size_X, &X);
		double* data = (double*)FLA_malloc(nElem * sizeof(dim_t));
		memset(data, 0, nElem * sizeof(dim_t));
		FLA_Obj_attach_buffer_to_tensor(data, order, stride_X, &X);

		//FOR
		FLA_Obj AT, AB;
		FLA_Obj A0, A1, A2;
		FLA_Obj CT, CB;
		FLA_Obj C0, C1, C2;
		
		FLA_Part_1xmode2(A, &AT,
							&AB, mode, 0, FLA_TOP);	
		FLA_Part_1xmode2(C, &CT,
							&CB, mode, 0, FLA_TOP);	
		//Only symmetric part touched
		//Ponder this
		while(i < maxIndex){
			dim_t b = 1;
			FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
										/**/ /**/
											  &A1,
										  AB, &A2, mode, b, FLA_BOTTOM); 
			FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
										/**/ /**/
											  &C1,
										  CB, &C2, mode, b, FLA_BOTTOM);

			memset((X.base)->buffer, 0, nElem * sizeof(dim_t));
			FLA_Ttm_single(alpha, A, mode, beta, B, X);

			FLA_Sttsm_single(alpha, X, mode-1, beta, B, C1, i);

			FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
												   C1,
											/********/
											  &CB, C2, mode, FLA_TOP);
			FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
												   A1,
											/********/
											  &AB, A2, mode, FLA_TOP);
			i++;
		}
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Sttsm( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	dim_t order = FLA_Obj_order(A);
	FLA_Sttsm_single( alpha, A, order, beta, B, C, order);

	return FLA_SUCCESS;
}
