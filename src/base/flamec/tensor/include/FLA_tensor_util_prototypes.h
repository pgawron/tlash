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
FLA_Error FLA_Random_tensor(FLA_Obj A);
FLA_Error FLA_Random_psym_tensor(FLA_Obj obj);

//---  Non-FLA utils ---------------
int compare_pairwise_sort(const void* a, const void* b);
int compare_dim_t(const void* a, const void* b);
dim_t binomial(dim_t n, dim_t k);
dim_t FLA_get_unique_info( TLA_sym sym, const dim_t index[], dim_t* sortedIndex, dim_t* permutation, dim_t* ipermutation);
FLA_Error FLA_Set_tensor_stride( dim_t order, const dim_t size[], dim_t* stride);
dim_t FLA_TIndex_to_LinIndex( dim_t order, dim_t const stride[], dim_t const index[]);
FLA_Error FLA_LinIndex_to_TIndex( dim_t order, dim_t const stride[], dim_t const linIndex, dim_t index[]);
dim_t FLA_Ttm_Ops( dim_t order, const dim_t size_A[], const dim_t size_B[2], dim_t mode);

//---  Array routines --------------
void print_array(const char* header, dim_t nElem, const dim_t arr[]);
dim_t FLA_array_product( dim_t order, const dim_t arr[]);
FLA_Error FLA_array_elemwise_product( dim_t order, const dim_t arr1[], const dim_t arr2[], dim_t* arrOut);
FLA_Error FLA_array_elemwise_quotient( dim_t order, const dim_t arr1[], const dim_t arr2[], dim_t* arrOut);

//--- Symmetry View routines ---------
FLA_Error TLA_create_part_obj( dim_t nPart, FLA_Obj* partitions[]);
FLA_Error TLA_destroy_part_obj( dim_t nPart, FLA_Obj* partitions[]);

//--- TLA_sym routines ------------
FLA_Error TLA_Sym_init_nonsymmetric(dim_t order, TLA_sym* sym);
