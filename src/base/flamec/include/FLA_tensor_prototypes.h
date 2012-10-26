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

//--------------------------------------------------------------------------
// TLASH stuff
//--------------------------------------------------------------------------

//--- Obj functions -----------------
FLA_Error FLA_Obj_create_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], FLA_Obj *obj);
FLA_Error FLA_Obj_create_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], FLA_Obj *obj );
FLA_Error FLA_Obj_create_blocked_tensor( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj);
FLA_Error FLA_Obj_create_blocked_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, dim_t size[order], dim_t size_inner[order], dim_t stride[order], dim_t blkSize[order], FLA_Obj *obj );
FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, dim_t stride[order], FLA_Obj *obj );
FLA_Error FLA_Obj_attach_buffer_to_symm_tensor( void *buffer[], dim_t order, dim_t stride[order], FLA_Obj *obj );
FLA_Error FLA_Obj_create_tensor_without_buffer( FLA_Datatype datatype, dim_t order, dim_t size[order], FLA_Obj *obj );
FLA_Error FLA_Obj_create_symm_tensor_without_buffer(FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t blkSize, FLA_Obj *obj);

//--- Query functions --------------

dim_t          FLA_Obj_order( FLA_Obj obj );
dim_t*         FLA_Obj_stride( FLA_Obj obj );
dim_t*         FLA_Obj_offset( FLA_Obj obj );
dim_t*         FLA_Obj_size( FLA_Obj obj );
dim_t*         FLA_Obj_permutation( FLA_Obj obj );
dim_t          FLA_Obj_dimsize( FLA_Obj obj, dim_t mode );
dim_t          FLA_Obj_base_dimstride( FLA_Obj obj, dim_t dim );
dim_t          FLA_Obj_dimstride( FLA_Obj obj, dim_t dim );


//--------------------------------------------------------------------------
FLA_Error     FLA_Part_1xmode2( FLA_Obj A,  FLA_Obj *A1,
                                            FLA_Obj *A2,
                                dim_t mode, dim_t  b,  FLA_Side side );

FLA_Error     FLA_Merge_1xmode2( FLA_Obj AT,
                                 FLA_Obj AB,  FLA_Obj *A, dim_t mode );

FLA_Error     FLA_Repart_1xmode2_to_1xmode3( FLA_Obj AT,  FLA_Obj *A0,
                                                          FLA_Obj *A1,
                                             FLA_Obj AB,  FLA_Obj *A2,
                                             dim_t mode,  dim_t    b,  FLA_Side side );

FLA_Error     FLA_Cont_with_1xmode3_to_1xmode2( FLA_Obj *AT,  FLA_Obj A0,
                                                              FLA_Obj A1,
                                                FLA_Obj *AB,  FLA_Obj A2,
                                                              dim_t mode, FLA_Side side );

FLA_Error     FLA_Part_1xmode2_check( FLA_Obj A,  FLA_Obj *A1,
                                                  FLA_Obj *A2,
                                      dim_t mode, dim_t  b,  FLA_Side side );

FLA_Error     FLA_Merge_1xmode2_check( FLA_Obj AT,
                                       FLA_Obj AB,  FLA_Obj *A, dim_t mode );

FLA_Error     FLA_Repart_1xmode2_to_1xmode3_check( FLA_Obj AT,  FLA_Obj *A0,
                                                                FLA_Obj *A1,
                                                   FLA_Obj AB,  FLA_Obj *A2,
                                                   dim_t mode,  dim_t    b,  FLA_Side side );

FLA_Error     FLA_Cont_with_1xmode3_to_1xmode2_check( FLA_Obj *AT,  FLA_Obj A0,
                                                                    FLA_Obj A1,
                                                      FLA_Obj *AB,  FLA_Obj A2,
                                                                    dim_t mode, FLA_Side side );

FLA_Error FLA_Obj_create_symm_tensor_without_buffer_check( FLA_Datatype datatype, dim_t order, dim_t size[order], dim_t b, FLA_Obj *obj );

//--------------------------------------------------------------------------

FLA_Error FLA_Check_attempted_repart_1xmode2( FLA_Obj A_side, dim_t mode, dim_t b );
FLA_Error FLA_Check_adjacent_objects_1xmode2( FLA_Obj AT,
                                              FLA_Obj AB, dim_t mode );

//--------------------------------------------------------------------------

FLA_Error FLA_Permute_single( FLA_Obj A, dim_t permtuation[], FLA_Obj* B);
FLA_Error FLA_Permute_hier( FLA_Obj A, dim_t permtuation[], FLA_Obj* B);
FLA_Error FLA_Permute_single_inplace( FLA_Obj* A, dim_t permutation[]);
FLA_Error FLA_Permute( dim_t permutation[], dim_t order, dim_t size[], FLA_Obj A[], dim_t* B_order, dim_t (*B_size)[], FLA_Obj (*B)[] );
FLA_Error FLA_Ttm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm_single_mode( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm( FLA_Obj alpha, FLA_Obj A, dim_t nModes, dim_t mode[nModes], FLA_Obj beta, FLA_Obj B[nModes], FLA_Obj C );
FLA_Error FLA_Ttm_hierC_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C);
FLA_Error FLA_Ttm_hierA_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C);
FLA_Error FLA_Sttsm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t maxIndex );
FLA_Error FLA_Sttsm( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C );

//---  Util functions ---------------
FLA_Error FLA_Set_zero_tensor( FLA_Obj A );
FLA_Error FLA_Adjust_2D_info( FLA_Obj *A );
FLA_Error FLA_Obj_create_Random_symm_tensor_data(dim_t b, FLA_Obj obj);

//---  Non-FLA utils ---------------
int compare_pairwise_sort(const void* a, const void* b);
FLA_Error FLA_get_unique_info( dim_t order, dim_t index[order], dim_t* sortedIndex, dim_t* permutation);
dim_t binomial(dim_t n, dim_t k);
FLA_Error FLA_Set_tensor_stride( dim_t order, dim_t size[order], dim_t* stride);
FLA_Error FLA_Set_tensor_permutation( dim_t order, dim_t permutation[order], FLA_Obj* A);
FLA_Error FLA_TIndex_to_LinIndex( dim_t order, dim_t stride[order], dim_t index[order], dim_t* linIndex);
FLA_Error FLA_Permute_array( dim_t order, dim_t arrfrom[order], dim_t permutation[order], dim_t* arrto);

//---  Needed sorting items ---------------
typedef struct pairwise_sort_struct{
  dim_t index;
  dim_t val;
} FLA_Paired_Sort;

//--- Misc functions ----------------------
FLA_Error FLA_Obj_print_tensor(FLA_Obj A);
FLA_Error FLA_Obj_print_tensor_in_mode(FLA_Obj A, dim_t mode);
