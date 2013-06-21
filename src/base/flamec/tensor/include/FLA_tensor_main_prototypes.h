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
FLA_Error FLA_Obj_create_tensor( FLA_Datatype datatype, dim_t order, const dim_t size[], const dim_t stride[], FLA_Obj *obj);
FLA_Error FLA_Obj_create_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, const dim_t size[], const dim_t size_inner[], const dim_t stride[], FLA_Obj *obj );
FLA_Error FLA_Obj_create_blocked_tensor_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t order, const dim_t size[], const dim_t size_inner[], const dim_t stride[], const dim_t blkSize[], FLA_Obj *obj );
FLA_Error FLA_Obj_create_psym_tensor( FLA_Datatype datatype, dim_t order, const dim_t size[], const dim_t stride[], TLA_sym sym, FLA_Obj *obj);
FLA_Error FLA_Obj_attach_buffer_to_tensor( void *buffer, dim_t order, const dim_t stride[], FLA_Obj *obj );
FLA_Error FLA_Obj_attach_buffer_to_blocked_tensor( void *buffer[], dim_t order, const dim_t stride[], FLA_Obj *obj );
FLA_Error FLA_Obj_attach_buffer_to_blocked_psym_tensor( void *buffer[], dim_t order, const dim_t stride[], FLA_Obj *obj );
FLA_Error FLA_Obj_create_tensor_without_buffer( FLA_Datatype datatype, dim_t order, const dim_t size[], FLA_Obj *obj );
FLA_Error FLA_Obj_create_blocked_tensor_without_buffer(FLA_Datatype datatype, dim_t order, const dim_t size[], const dim_t blkSize[], FLA_Obj *obj);
FLA_Error FLA_Obj_create_blocked_psym_tensor_without_buffer(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blk_size[], TLA_sym sym, FLA_Obj *obj);
FLA_Error FLA_Obj_create_blocked_tensor(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blocked_stride[], const dim_t blk_size[], FLA_Obj *obj);
FLA_Error FLA_Obj_create_blocked_psym_tensor(FLA_Datatype datatype, dim_t order, const dim_t flat_size[], const dim_t blocked_stride[], const dim_t blk_size[], TLA_sym sym, FLA_Obj *obj);

//--- Query functions --------------

dim_t		FLA_Obj_order( FLA_Obj obj );
dim_t*		FLA_Obj_stride( FLA_Obj obj );
dim_t*		FLA_Obj_offset( FLA_Obj obj );
dim_t*		FLA_Obj_size( FLA_Obj obj );
dim_t*		FLA_Obj_permutation( FLA_Obj obj );
dim_t		FLA_Obj_dimsize( FLA_Obj obj, dim_t mode );
dim_t		FLA_Obj_base_dimstride( FLA_Obj obj, dim_t dim );
dim_t		FLA_Obj_dimstride( FLA_Obj obj, dim_t dim );
FLA_Error	FLA_Obj_blocked_tensor_free_buffer(FLA_Obj* obj);
FLA_Error	FLA_Obj_blocked_psym_tensor_free_buffer( FLA_Obj *obj);
dim_t*		FLA_Obj_base_scalar_size(FLA_Obj A);
dim_t		FLA_Obj_base_scalar_dimsize(FLA_Obj A, dim_t mode);
void*		FLA_Obj_tensor_buffer_at_view( FLA_Obj obj );

//--- Symmetry related queries --------------------------
dim_t		TLA_mode_at_sym_pos( TLA_sym S, dim_t pos );
dim_t		TLA_sym_group_of_pos( TLA_sym S, dim_t pos );
dim_t		TLA_sym_group_size( TLA_sym S, dim_t symgroup);
dim_t		TLA_sym_pos_of_mode( TLA_sym S, dim_t mode);
dim_t       TLA_sym_group_of_mode( TLA_sym S, dim_t mode);
dim_t       TLA_sym_group_mode_offset( TLA_sym S, dim_t symGroup);
dim_t       TLA_sym_group_of_mode_size( TLA_sym S, dim_t mode);
FLA_Error   TLA_update_sym_based_offset( TLA_sym S, FLA_Obj* A);
FLA_Error   TLA_split_sym_group(TLA_sym S, dim_t nSplit_modes, const dim_t split_modes[], TLA_sym* S1);

//--------------------------------------------------------------------------
// --- FLA_View functions (non-symmetric) -------------------------------
FLA_Error   FLA_Part_1xmode2( FLA_Obj A, FLA_Obj *A1,
                                         /**/
                                         FLA_Obj *A2,
                              dim_t mode, dim_t b, FLA_Side side );


FLA_Error	FLA_Merge_1xmode2( FLA_Obj AT,
         	                   /**/
                               FLA_Obj AB,  FLA_Obj *A,
                               dim_t mode );


FLA_Error	FLA_Repart_1xmode2_to_1xmode3( FLA_Obj AT,  FLA_Obj *A0,
                                                        FLA_Obj *A1,
                                           FLA_Obj AB,  FLA_Obj *A2,
                                           dim_t mode,  dim_t    b,  FLA_Side side );


FLA_Error	FLA_Cont_with_1xmode3_to_1xmode2( FLA_Obj *AT,  FLA_Obj A0,
                                                            FLA_Obj A1,
                                              FLA_Obj *AB,  FLA_Obj A2,
                                                            dim_t mode, FLA_Side side );
	

//--------------------------------------------------------------------------
// --- FLA_View functions (symmetric) -------------------------------


FLA_Error   FLA_Part_2powm( FLA_Obj A,  FLA_Obj* Apart[],
                            dim_t nModes_repart, const dim_t repart_modes[],
                            const dim_t sizes[], const FLA_Side sides[]);

                                       //C won't allow
                                       //FLA_Obj const * const Apart[]
FLA_Error   FLA_Repart_2powm_to_3powm( FLA_Obj * Apart[], FLA_Obj* Arepart[],
                                       dim_t nModes_repart, const dim_t repart_modes[],
                                       const dim_t sizes[], const FLA_Side sides[]);

                          //C won't allow
                          //FLA_Obj const * const Apart[]
FLA_Error   FLA_Merge_2powm(FLA_Obj* Apart[], FLA_Obj* A,
                            dim_t nModes_repart, const dim_t repart_modes[]);

                                                          //C won't allow
                                                          //FLA_Obj const * const Arepart[]
FLA_Error   FLA_Cont_with_3powm_to_2powm( FLA_Obj* Apart[], FLA_Obj* Arepart[],
                                          dim_t nModes_repart, const dim_t repart_modes[],
                                          const FLA_Side sides[]);

// -------------------------------------------------------------------------
// ---  Check functions

FLA_Error	FLA_Obj_create_blocked_sym_tensor_without_buffer_check( FLA_Datatype datatype, dim_t order, const dim_t size[], dim_t b, FLA_Obj *obj );

FLA_Error   FLA_Repart_1xmode2_to_1xmode3_check( FLA_Obj AT,  FLA_Obj *A0,
                                                              FLA_Obj *A1,
                                                 FLA_Obj AB,  FLA_Obj *A2,
                                                 dim_t mode,  dim_t    b,  FLA_Side side );

FLA_Error   FLA_Cont_with_1xmode3_to_1xmode2_check( FLA_Obj *AT,  FLA_Obj A0,
                                                                  FLA_Obj A1,
                                                    FLA_Obj *AB,  FLA_Obj A2,
                                                    dim_t mode, FLA_Side side );


FLA_Error   FLA_Part_1xmode2_check( FLA_Obj A,  FLA_Obj *A1,
                                                FLA_Obj *A2,
                                    dim_t mode, dim_t  b,  FLA_Side side );

FLA_Error   FLA_Merge_1xmode2_check( FLA_Obj AT,
                                     FLA_Obj AB,  FLA_Obj *A, dim_t mode );

// -------------------------------------------------------------------------
FLA_Error	FLA_Check_attempted_repart_1xmode2( FLA_Obj A_side, dim_t mode, dim_t b );
FLA_Error	FLA_Check_adjacent_objects_1xmode2( FLA_Obj AT,
                                                FLA_Obj AB, dim_t mode );

//---  Needed sorting items ---------------
typedef struct pairwise_sort_struct{
  dim_t index;
  dim_t val;
} FLA_Paired_Sort;

//--- Misc functions ----------------------
FLA_Error	FLA_Obj_print_hier_tensor_repart_mode_at(FLA_Obj A, dim_t repart_mode, dim_t index[]);
FLA_Error	FLA_Obj_print_tensor(FLA_Obj A);
FLA_Error	FLA_Obj_print_flat_tensor(FLA_Obj A);
FLA_Error   FLA_Obj_print(FLA_Obj A);
FLA_Error   FLA_Obj_print_matlab(const char * varName, FLA_Obj obj);
