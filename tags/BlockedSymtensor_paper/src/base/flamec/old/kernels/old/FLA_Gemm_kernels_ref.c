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

void FLA_Gemm_pack_C( FLA_Trans transC, FLA_Obj C, FLA_Obj *packed_C )
{
  *packed_C = C;
}

void FLA_Gemm_unpack_andor_scale_C( FLA_Trans transC, FLA_Obj alpha, 
				    FLA_Obj C, FLA_Obj *packed_C )
{
}

void FLA_Gemm_pack_andor_scale_B( FLA_Trans transB, FLA_Obj alpha,
				  FLA_Obj B, FLA_Obj *packed_B )
{
  *packed_B = B;
}

void FLA_Gemm_release_pack_B( FLA_Trans transB, FLA_Obj *packed_B )
{
}

void FLA_Gemm_pack_andor_scale_A( FLA_Trans transA, FLA_Obj alpha,
				  FLA_Obj A, FLA_Obj *packed_A )
{
  *packed_A = A;
}

void FLA_Gemm_release_pack_A( FLA_Trans transA, FLA_Obj *packed_A )
{
}

void FLA_Gemm_kernel( FLA_Obj alpha, FLA_Obj packed_A, 
		      FLA_Obj packed_B, FLA_Obj packed_C )
{
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, packed_A, packed_B,
	    FLA_ONE, packed_C );
}
