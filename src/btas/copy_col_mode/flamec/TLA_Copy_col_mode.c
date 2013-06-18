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

#include "TLA_Copy_col_mode.h"

//Copies a single column from A to B
FLA_Error TLA_Copy_col_mode(FLA_Obj A, dim_t mode_A, FLA_Obj B, dim_t mode_B){
	dim_t i;
	dim_t stride_A = FLA_Obj_dimstride(A,mode_A);
	dim_t stride_B = FLA_Obj_dimstride(B,mode_B);

	/*
	 *	Replace the following with better routine
	 */
	dim_t dt_A = FLA_Obj_datatype(A);
	if(dt_A == FLA_DOUBLE){
		double* buf_A = FLA_Obj_tensor_buffer_at_view(A);
		double* buf_B = FLA_Obj_tensor_buffer_at_view(B);
		for(i = 0; i < FLA_Obj_dimsize(A,mode_A); i++){
			buf_B[i*stride_B] = buf_A[i*stride_A];
		}
	}else{
		FLA_Obj* buf_A = FLA_Obj_tensor_buffer_at_view(A);
		FLA_Obj* buf_B = FLA_Obj_tensor_buffer_at_view(B);
		for(i = 0; i < FLA_Obj_dimsize(A,mode_A); i++){
			buf_B[i*stride_B] = buf_A[i*stride_A];
		}
	}
	return FLA_SUCCESS;
}
