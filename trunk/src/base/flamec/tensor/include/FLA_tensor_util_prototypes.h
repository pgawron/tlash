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

//---  Util functions ---------------

FLA_Error FLA_Set_zero_tensor( FLA_Obj A );
FLA_Error FLA_Adjust_2D_info( FLA_Obj *A );
FLA_Error FLA_Obj_create_Random_symm_tensor_data(dim_t b, FLA_Obj obj);
FLA_Error FLA_Random_dense_symm_tensor(dim_t nSymmGroups, dim_t symmGroupLens[nSymmGroups], dim_t** symmetries, FLA_Obj *obj);
FLA_Error FLA_Random_tensor(FLA_Obj A);

//---  Non-FLA utils ---------------
int compare_pairwise_sort(const void* a, const void* b);
dim_t binomial(dim_t n, dim_t k);
FLA_Error FLA_get_unique_info( dim_t order, dim_t index[order], dim_t* sortedIndex, dim_t* permutation);
FLA_Error FLA_Set_tensor_stride( dim_t order, dim_t size[order], dim_t* stride);
FLA_Error FLA_Set_tensor_permutation( dim_t order, dim_t permutation[order], FLA_Obj* A);
FLA_Error FLA_TIndex_to_LinIndex( dim_t order, dim_t stride[order], dim_t index[order], dim_t* linIndex);
FLA_Error FLA_Permute_array( dim_t order, dim_t arrfrom[order], dim_t permutation[order], dim_t* arrto);
dim_t FLA_Ttm_Ops( dim_t order, dim_t size_A[order], dim_t size_B[2], dim_t mode);

