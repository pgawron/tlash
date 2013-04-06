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

//
// --- FLA_Part_2powm() --------------------------------------------------------
//
//TODO: Quadrant not yet generalized.  For now only handles TTTTTT... quadrant

FLA_Error FLA_Part_2powm( FLA_Obj A, FLA_Obj* Apart[],
                          dim_t nModes_part, dim_t part_modes[nModes_part],
                          dim_t sizes[], FLA_Side sides[]){

    dim_t i,j,k;
    dim_t num_part = 1 << nModes_part;

    const dim_t* mode_size = FLA_Obj_size(A);
    const dim_t* mode_offset = FLA_Obj_offset(A);
    dim_t part_mode_stride = num_part / 2;              //mode stride of part array

    //For each mode, we update all regions appropriately
    //Apart regions are laid out in column-major order, so each offset of the regions in linear array
    //can be viewed as groups of 0s followed by 1s
    //For instance, offset-values of an order-3 tensor looks like this in memory:
    //mode 0: [0 1 0 1 0 1 0 1]
    //mode 1: [0 0 1 1 0 0 1 1]
    //mode 2: [0 0 0 0 1 1 1 1]


    //All regions have base pointer of A.base

    for(i = 0; i < num_part; i++){
        Apart[i]->order = A.order;
        Apart[i]->base = A.base;

        //This will change when hierarchy comes in.  For now, we need to keep consistent sizes and offsets
        //for ALL modes (not just ones we are partitioning)
        memcpy(&(((Apart[i])->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
        memcpy(&(((Apart[i])->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
    }

    for(i = nModes_part - 1; i < nModes_part; i--){
        dim_t ind = 0;
        for(j = 0; j < num_part / (part_mode_stride * 2); j++){
            //If this region has a 0 in this mode of its offset
            for(k = 0; k < part_mode_stride; k++){
                ((Apart[ind])->offset)[part_modes[i]] = mode_offset[part_modes[i]];
                ((Apart[ind])->size)[part_modes[i]] = sizes[i];
                ind++;
            }
            //If this region has a 1 in this mode of its offset
            for(k = 0; k < part_mode_stride; k++){
                Apart[ind]->offset[part_modes[i]] = mode_offset[part_modes[i]] + sizes[i];
                Apart[ind]->size[part_modes[i]] = mode_size[part_modes[i]] - sizes[i];
                ind++;
            }
        }
        part_mode_stride /= 2;
    }

    //All regions need to update the symmetry of the objects
    for(i = 0; i < num_part; i++){
        memcpy(&(((Apart[i])->permutation)[0]), &(A.permutation[0]), FLA_Obj_order(A) * sizeof(dim_t));
        TLA_update_sym_based_offset(A.sym, Apart[i]);
    }

    return FLA_SUCCESS;
}
/*
FLA_Error FLA_Part_2powm_helper(FLA_Obj A, FLA_Obj* const Apart[],
                                dim_t nModes_repart, dim_t repart_modes[nModes_repart],
                                dim_t sizes[], FLA_Side sides[],
                                dim_t stride[], TLA_sym origSym, dim_t mode_index){
	dim_t i,j;
    dim_t order = FLA_Obj_order(A);
    dim_t size = FLA_Obj_dimsize(A,repart_modes[0]);

    *(Apart[0]) = A;
	dim_t stride_src = size;

    for(i = 0; i < order; i++){
        dim_t stride_part = stride_src / 2;
        for(j = 0; j < size; j += stride_src){
            FLA_Part_1xmode2(*(Apart[j]), Apart[j],
                                       Apart[j + stride_part],
                                       repart_modes[i], sizes[i], sides[i]);
        }
        stride_src /= 2;
    }

    dim_t num_elem_alloc = FLA_Obj_num_elem_alloc(A);
    for(i = 0; i < num_elem_alloc; i++){
        TLA_update_sym_based_offset(origSym, Apart[i]);
    }

	return FLA_SUCCESS;
}

FLA_Error FLA_Part_2powm( FLA_Obj A, FLA_Obj* const Apart[],
                          dim_t nModes_repart, dim_t repart_modes[nModes_repart],
                          dim_t sizes[], FLA_Side sides[]){

	dim_t order = FLA_Obj_order(A);
	dim_t* stride = FLA_Obj_stride(A);
	FLA_Part_2powm_helper( A, Apart,
	                      nModes_repart, repart_modes,
	                      sizes, sides, stride, A.sym, order - 1);

	FLA_free(stride);

	return FLA_SUCCESS;
}
*/

//
// --- FLA_Part_1xmode2() ----------------------------------------------------------
//

FLA_Error FLA_Part_1xmode2( FLA_Obj A,  FLA_Obj *A1,
                                        FLA_Obj *A2,
                            dim_t mode, dim_t  b,  FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Part_1xmode2_check( A,    A1,
                                  A2, mode, b, side );

  // Safeguard: if mb > m, reduce mb to m.
  if ( b > A.size[mode] ) b = A.size[mode];

  // Set mb to be the dimension of A1.
  if ( side == FLA_BOTTOM ) b = A.size[mode] - b;

  A1->order = A.order;

  memcpy(&((A1->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
  (A1->size)[mode] = b;
  memcpy(&((A1->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
  A1->base = A.base;
  memcpy(&((A1->permutation)[0]), &(A.permutation[0]), A.order * sizeof(dim_t));

  A2->order = A.order;
  memcpy(&((A2->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
  (A2->size)[mode] = A.size[mode] - b;
  memcpy(&((A2->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
  (A2->offset)[mode] += b;
  A2->base = A.base;
  memcpy(&((A2->permutation)[0]), &(A.permutation[0]), A.order * sizeof(dim_t));

  //Update symmetries

  dim_t split_mode_arr[1];
  split_mode_arr[0] = mode;
  TLA_split_sym_group(A.sym, 1, split_mode_arr, &(A1->sym));
  TLA_split_sym_group(A.sym, 1, split_mode_arr, &(A2->sym));

  FLA_Adjust_2D_info(A1);
  FLA_Adjust_2D_info(A2);

  return FLA_SUCCESS;
}

//
// --- FLA_Repart_2powm_to_3powm() -------------------------------------------------
//
//Note: Only handles BBBBBB... for now
//believed to be easy to change, each T corresponds to a 0, B a 1.
//repartition all regions with 2powm_to_3powm except for the region corresponding to the Side->index mapping
//Manually update that region
//TODO: See if a hierarchical method is better
//NOTE: Only works for BBBBBBBBBBB  Should be straightforward to generalize
//1. Setup hierarchy to turn object into 2x1 in mode working with
//2. Repartition 2x1 to 3x1 in that mode
//3. Repeat for all modes
                                   //C won't allow
                                   //FLA_Obj const * const Apart[]
FLA_Error FLA_Repart_2powm_to_3powm( FLA_Obj* Apart[],  FLA_Obj* Arepart[],
                                     dim_t nModes_repart, dim_t repart_modes[nModes_repart],
                                     dim_t sizes[nModes_repart], FLA_Side sides[nModes_repart])
{
    dim_t i,j,k;
    dim_t const num_parts = 1 << nModes_repart;
    dim_t num_reparts = 1;
    for(i = 0; i < nModes_repart; i++)
        num_reparts *= 3;

    dim_t part_mode_stride = num_parts / 2;
    dim_t repart_mode_stride = num_reparts / 3;

    //Update base pointer for all repart regions
    dim_t repart_ind = 0;
    for(i = 0; i < num_parts; i+=2){
        Arepart[repart_ind]->order = Apart[i]->order;
        Arepart[repart_ind]->base = Apart[i]->base;
        Arepart[repart_ind]->sym = Apart[i]->sym;

        //This will change when hierarchy comes in.  For now, we need to keep consistent sizes and offsets
        //for ALL modes (not just ones we are partitioning)
        memcpy(&(((Arepart[repart_ind])->offset)[0]), &(((Apart[i])->offset)[0]), (Apart[i])->order * sizeof(dim_t));
        memcpy(&(((Arepart[repart_ind])->size)[0]), &(((Apart[i])->size)[0]), (Apart[i])->order * sizeof(dim_t));
        memcpy(&(((Arepart[repart_ind])->permutation)[0]), &(((Apart[i])->permutation)[0]), (Apart[i])->order * sizeof(dim_t));
        repart_ind++;

        //Indices with 1 and 2 share the same base
        Arepart[repart_ind]->order = Apart[i+1]->order;
        Arepart[repart_ind+1]->order = Apart[i+1]->order;
        Arepart[repart_ind]->sym = Apart[i+1]->sym;
        Arepart[repart_ind]->base = Apart[i+1]->base;
        Arepart[repart_ind+1]->base = Apart[i+1]->base;
        Arepart[repart_ind]->sym = Apart[i+1]->sym;
        Arepart[repart_ind+1]->sym = Apart[i+1]->sym;
        //This will change when hierarchy comes in.  For now, we need to keep consistent sizes and offsets
        //for ALL modes (not just ones we are partitioning)
        memcpy(&(((Arepart[repart_ind])->offset)[0]), &(((Apart[i+1])->offset)[0]), (Apart[i+1])->order * sizeof(dim_t));
        memcpy(&(((Arepart[repart_ind])->size)[0]), &(((Apart[i+1])->size)[0]), (Apart[i+1])->order * sizeof(dim_t));
        memcpy(&(((Arepart[repart_ind+1])->offset)[0]), &(((Apart[i+1])->offset)[0]), (Apart[i+1])->order * sizeof(dim_t));
        memcpy(&(((Arepart[repart_ind+1])->size)[0]), &(((Apart[i+1])->size)[0]), (Apart[i+1])->order * sizeof(dim_t));
        repart_ind+=2;
    }

    for(i = nModes_repart - 1; i < nModes_repart; i--){
        dim_t repart_ind = 0;
        dim_t part_ind = 0;
        for(j = 0; j < num_reparts / (repart_mode_stride * 3); j++){
            //If this region has a 0 in this mode of its offset
            for(k = 0; k < repart_mode_stride; k++){
                ((Arepart[repart_ind])->offset)[repart_modes[i]] = ((Apart[part_ind])->offset)[repart_modes[i]];
                ((Arepart[repart_ind])->size)[repart_modes[i]] = ((Apart[part_ind])->size)[repart_modes[i]];
                repart_ind++;
            }
            part_ind += part_mode_stride;

            //If this region has a 1 in this mode of its offset
            for(k = 0; k < repart_mode_stride; k++){
                Arepart[repart_ind]->offset[repart_modes[i]] = Apart[part_ind]->offset[repart_modes[i]];
                Arepart[repart_ind]->size[repart_modes[i]] = sizes[i];
                repart_ind++;
            }

            //If this region has a 2 in this mode of its offset
            for(k = 0; k < repart_mode_stride; k++){
                Arepart[repart_ind]->offset[repart_modes[i]] = Apart[part_ind]->offset[repart_modes[i]] + sizes[i];
                Arepart[repart_ind]->size[repart_modes[i]] = Apart[part_ind]->size[repart_modes[i]] - sizes[i];
                repart_ind++;
            }
        }
        part_mode_stride /= 2;
        repart_mode_stride /= 3;
    }

    return FLA_SUCCESS;
}

/*
FLA_Error FLA_Repart_2powm_to_3powm( FLA_Obj const * const Apart[],  FLA_Obj* const Arepart[],
                                     dim_t nModes_repart, dim_t repart_modes[nModes_repart],
                                 	 dim_t sizes[nModes_repart], FLA_Side sides[nModes_repart])
{
    dim_t i, j;

    dim_t order = nModes_repart;
    dim_t num_Apart = 1 << order;
    dim_t num_Arepart = 1;

    dim_t stride_part[order];
    for(i = 0; i < order; i++){
        stride_part[i] = 2;
        stride_part[i] = 3;
        num_Arepart *= 3;
    }

    FLA_Obj * part_storage[num_Apart];
    dim_t nParts;
    dim_t part_sizes[nModes_repart];
    FLA_Side part_sides[nModes_repart];
    dim_t part_modes[nModes_repart];

    for(i = 0; i < num_Apart; i++){

        dim_t repart_lin_offset;
        FLA_Obj const * curObj;

        FLA_TIndex_to_LinIndex(order, curObj->offset, stride_part, &repart_lin_offset);

        dim_t size_storage_needed = 1;
        nParts = 0;
        for(j = 0; j < FLA_Obj_order(*curObj); j++){
            if((curObj->offset)[j] == 1){
                size_storage_needed *= 2;
                part_sizes[nParts] = sizes[repart_modes[j]];
                part_modes[nParts] = repart_modes[j];
                part_sides[nParts] = sides[repart_modes[j]];
                nParts++;
            }
        }

        dim_t stride_storage = size_storage_needed;
        dim_t next_stride_storage = size_storage_needed / 2;
        dim_t stride_repart = num_Arepart;
        dim_t next_stride_repart = num_Arepart / 3;

        part_storage[0] = Arepart[repart_lin_offset];

        for(j = 0; j < order; j++){
            if((curObj->offset)[j] == 1){
                part_storage[j*stride_storage + next_stride_storage] = Arepart[repart_lin_offset + j * stride_repart + next_stride_repart];
            }
            stride_storage /= 2;
            next_stride_storage /= 2;
            stride_repart /= 3;
            next_stride_repart /= 3;
        }


        FLA_Part_2powm(*curObj, part_storage, nParts, part_modes, part_sizes, part_sides);
    }

    return FLA_SUCCESS;
}
*/
/*
FLA_Error FLA_Repart_2xmodeBlank_to_3xmodeBlank( dim_t order, dim_t sizeApart[], dim_t strideApart[], FLA_Obj Apart[], 
												 			  dim_t strideArepart[], FLA_Obj Arepart[],
										 		 			  dim_t mode, dim_t size, FLA_Side side)
{
	dim_t baseOffset[order];

	memset(&(baseOffset[0]), 0, order * sizeof(dim_t));

	dim_t update_ptr = order - 1;
	while(TRUE){
		dim_t inObj0, inObj1;
		FLA_TIndex_to_LinIndex(order, baseOffset, strideApart, &inObj0);
		inObj1 = inObj0 + strideApart[mode];

		dim_t outObj0, outObj1, outObj2;
		FLA_TIndex_to_LinIndex(order, baseOffset, strideArepart, &outObj0);
		outObj1 = outObj0 + strideArepart[mode];
		outObj2 = outObj1 + strideArepart[mode];

		FLA_Repart_1xmode2_to_1xmode3( Apart[inObj0], &(Arepart[outObj0]),
													  &(Arepart[outObj1]),
									   Apart[inObj1], &(Arepart[outObj2]),
									   mode, size, side);

		//Update
		baseOffset[update_ptr]++;
		while(update_ptr < order && baseOffset[update_ptr] == sizeApart[update_ptr]){
			update_ptr--;
			if(update_ptr == mode){
				if(mode == 0)
					break;
				else
					update_ptr--;
			}
			if(update_ptr < order)
				baseOffset[update_ptr]++;
		}
		if(update_ptr >= order)
			break;
		for(dim_t i = update_ptr+1; i < order; i++)
			baseOffset[i] = 0;
		update_ptr = order - 1;
	}

	return FLA_SUCCESS;
}
*/

//
// --- FLA_Repart_1xmode2_to_1xmode3() -----------------------------------------
//

FLA_Error FLA_Repart_1xmode2_to_1xmode3( FLA_Obj AT,   FLA_Obj *A0,
                                                       FLA_Obj *A1,
                                         FLA_Obj AB,   FLA_Obj *A2,
                                         dim_t   mode, dim_t    b,
                                         FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Repart_1xmode2_to_1xmode3_check( AT,     A0,
                                                 A1,

                                         AB,     A2, mode, b, side );

  if ( side == FLA_TOP )
  {
    FLA_Part_1xmode2 ( AT,    A0,
                              A1,    mode, b, FLA_BOTTOM );

    A2->order = AB.order;
    memcpy(&((A2->size)[0]), &(AB.size[0]), AB.order * sizeof(dim_t));
    memcpy(&((A2->offset)[0]), &(AB.offset[0]), AB.order * sizeof(dim_t));
	memcpy(&((A2->permutation)[0]), &(AB.permutation[0]), AB.order * sizeof(dim_t));
    A2->base = AB.base;

    A2->sym = AB.sym;
    //A2->nSymGroups = AB.nSymGroups;
    //memcpy(&((A2->symGroupLens)[0]), &(AB.symGroupLens[0]), AB.nSymGroups * sizeof(dim_t));
    //memcpy(&((A2->symModes)[0]), &(AB.symModes[0]), AB.order * sizeof(dim_t));
  }
  else
  {
    A0->order = AT.order;
    memcpy(&((A0->size)[0]), &(AT.size[0]), AT.order * sizeof(dim_t));
    memcpy(&((A0->offset)[0]), &(AT.offset[0]), AT.order * sizeof(dim_t));
	memcpy(&((A0->permutation)[0]), &(AT.permutation[0]), AT.order * sizeof(dim_t));
    A0->base = AT.base;

    A0->sym = AT.sym;
    //A0->nSymGroups = AT.nSymGroups;
    //memcpy(&((A0->symGroupLens)[0]), &(AT.symGroupLens[0]), AT.nSymGroups * sizeof(dim_t));
    //memcpy(&((A0->symModes)[0]), &(AT.symModes[0]), AT.order * sizeof(dim_t));

    FLA_Part_1xmode2 ( AB,    A1,
                              A2,    mode, b, FLA_TOP );
  }

  FLA_Adjust_2D_info(A0);
  FLA_Adjust_2D_info(A1);
  FLA_Adjust_2D_info(A2);
	
  return FLA_SUCCESS;
}


//
// --- FLA_Cont_with_3powm_to_2powm() -------------------------------------------------
//
//TODO: See if a hierarchical method is better
//NOTE: ONLY FOR TTTTTTTTT.... Should be straightforward to switch
//1. Setup hierarchy to turn object into 2x1 in mode working with
//2. Repartition 2x1 to 3x1 in that mode
//3. Repeat for all modes
                                                          //C won't allow
                                                          //FLA_Obj const * const Arepart[]
FLA_Error FLA_Cont_with_3powm_to_2powm( FLA_Obj* Apart[], FLA_Obj * Arepart[],
                                        dim_t nModes_cont_with, dim_t cont_with_modes[nModes_cont_with],
                                        FLA_Side sides[nModes_cont_with]){
    dim_t i,j,k;
    dim_t const num_part = 1 << nModes_cont_with;
    dim_t num_repart = 1;
    for(i = 0; i < nModes_cont_with; i++)
        num_repart *= 3;

    dim_t part_mode_stride = num_part / 2;
    dim_t repart_mode_stride = num_repart / 3;

    //Update offsets and base pointers
    //This will change when hierarchy introduces, for now we need to keep all sizes and offsets consistent
    //regardless of what modes we are repartitioning
    dim_t repart_ind = 0;
    for(i = 0; i < num_part; i += 2){
        (Apart[i])->order = (Arepart[repart_ind])->order;
        memcpy(&(((Apart[i])->offset)[0]), &(((Arepart[repart_ind])->offset)[0]), (Arepart[repart_ind])->order * sizeof(dim_t));
        memcpy(&(((Apart[i])->permutation)[0]), &(((Arepart[repart_ind])->permutation)[0]), (Arepart[repart_ind])->order * sizeof(dim_t));
        Apart[i]->sym = Arepart[repart_ind]->sym;
        repart_ind += 2;
        (Apart[i+1])->order = (Arepart[repart_ind])->order;
        memcpy(&(((Apart[i+1])->offset)[0]), &(((Arepart[repart_ind])->offset)[0]), (Arepart[repart_ind])->order * sizeof(dim_t));
        memcpy(&(((Apart[i+1])->permutation)[0]), &(((Arepart[repart_ind])->permutation)[0]), (Arepart[repart_ind])->order * sizeof(dim_t));
        Apart[i+1]->sym = Arepart[repart_ind]->sym;
        repart_ind++;
    }


    for(i = nModes_cont_with - 1; i < nModes_cont_with; i--){
        dim_t part_ind = 0;
        for(j = 0; j < num_part / (part_mode_stride*2); j++){
            for(k = 0; k < part_mode_stride; k++){
                ((Apart[part_ind])->size)[cont_with_modes[i]] = ((Arepart[0])->size)[cont_with_modes[i]] + ((Arepart[repart_mode_stride])->size)[cont_with_modes[i]];
                part_ind++;
            }

            for(k = 0; k < part_mode_stride; k++){
                ((Apart[part_ind])->size)[cont_with_modes[i]] = ((Arepart[2*repart_mode_stride])->size)[cont_with_modes[i]];
                part_ind++;
            }
        }
        part_mode_stride /= 2;
        repart_mode_stride /= 3;
    }

    //Update offsets and base pointers
    repart_ind = 0;
    for(i = 0; i < num_part; i += 2){
        for(j = 0; j < nModes_cont_with; j++){
            ((Apart[i])->offset)[cont_with_modes[j]] = ((Arepart[repart_ind])->offset)[cont_with_modes[j]];
        }
        repart_ind += 2;
        for(j = 0; j < nModes_cont_with; j++){
            ((Apart[i+1])->offset)[cont_with_modes[j]] = ((Arepart[repart_ind])->offset)[cont_with_modes[j]];
        }
        repart_ind++;
    }
    return FLA_SUCCESS;
}

/*
FLA_Error FLA_Cont_with_3powm_to_2powm( FLA_Obj* const Apart[], FLA_Obj const * const Arepart[],
                                        dim_t nModes_repart, dim_t repart_modes[nModes_repart],
                                        FLA_Side sides[nModes_repart]){
    dim_t i;
    dim_t size_Arepart[nModes_repart];
    for(i = 0; i < nModes_repart; i++){
        size_Arepart[i] = 3;
    }

    dim_t stride_Apart[nModes_repart];
    dim_t stride_Arepart[nModes_repart];
    stride_Apart[0] = 1;
    stride_Arepart[0] = 1;
    for(i = 1; i < nModes_repart; i++){
        stride_Arepart[i] = stride_Arepart[i-1] * 3;
        stride_Apart[i] = stride_Apart[i-1] * 2;
    }

    FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( nModes_repart, stride_Apart, Apart,
                                                     size_Arepart, stride_Arepart, Arepart,
                                           repart_modes[0], sides[0] );
    for(i = 1; i < nModes_repart; i++){
        dim_t newSizeArepart[nModes_repart];
        memcpy(&(newSizeArepart[0]), &(size_Arepart[0]), nModes_repart * sizeof(dim_t));
        newSizeArepart[repart_modes[i]] = 2;
        FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( nModes_repart, stride_Apart, Apart,
                                                  newSizeArepart, stride_Apart, Apart,
                                                  repart_modes[i], sides[i] );
    }

    return FLA_SUCCESS;
}
*/
/*
FLA_Error FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( dim_t order, dim_t strideApart[], FLA_Obj Apart[],
												 			  dim_t sizeArepart[], dim_t strideArepart[], FLA_Obj Arepart[],
										 		 			  dim_t mode, FLA_Side side)
{
	dim_t baseOffset[order];

	memset(&(baseOffset[0]), 0, order * sizeof(dim_t));

	dim_t update_ptr = order - 1;
	while(TRUE){
		dim_t inObj0, inObj1, inObj2;
		FLA_TIndex_to_LinIndex(order, baseOffset, strideArepart, &inObj0);
		inObj1 = inObj0 + strideArepart[mode];
		inObj2 = inObj1 + strideArepart[mode];

		dim_t outObj0, outObj1;
		FLA_TIndex_to_LinIndex(order, baseOffset, strideApart, &outObj0);
		outObj1 = outObj0 + strideApart[mode];

		FLA_Cont_with_1xmode3_to_1xmode2( &(Apart[outObj0]), Arepart[inObj0],
													  Arepart[inObj1],
									   &(Apart[outObj1]), Arepart[inObj2],
									   mode, side);

		//Update
		baseOffset[update_ptr]++;
		while(update_ptr < order && baseOffset[update_ptr] == sizeArepart[update_ptr]){
			update_ptr--;
			if(update_ptr == mode){
				if(mode == 0)
					break;
				else
					update_ptr--;
			}
			if(update_ptr < order)
				baseOffset[update_ptr]++;
		}
		if(update_ptr >= order)
			break;
		for(dim_t i = update_ptr+1; i < order; i++)
			baseOffset[i] = 0;
		update_ptr = order - 1;
	}
	return FLA_SUCCESS;
}
*/

//
// --- FLA_Cont_with_1xmode3_to_1xmode2() ----------------------------------------------
//

FLA_Error FLA_Cont_with_1xmode3_to_1xmode2( FLA_Obj *AT,  FLA_Obj A0,
                                                          FLA_Obj A1,
                                            FLA_Obj *AB,  FLA_Obj A2,
                                                          dim_t mode, FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Cont_with_1xmode3_to_1xmode2_check( AT,     A0,
                                                    A1,
                                            AB,     A2, mode, side );

  AT->sym = A0.sym;

  if ( side == FLA_TOP )
  {
    AT->order = A0.order;
    memcpy(&((AT->size)[0]), &(A0.size[0]), A0.order * sizeof(dim_t));
    AT->size[mode] += A1.size[mode];
    memcpy(&((AT->offset)[0]), &(A0.offset[0]), A0.order * sizeof(dim_t));
    AT->base = A0.base;
    memcpy(&((AT->permutation)[0]), &(A0.permutation[0]), A0.order * sizeof(dim_t));

    AB->order = A2.order;
    memcpy(&((AB->size)[0]), &(A2.size[0]), A2.order * sizeof(dim_t));
    memcpy(&((AB->offset)[0]), &(A2.offset[0]), A2.order * sizeof(dim_t));
    AB->base = A2.base;
    memcpy(&((AB->permutation)[0]), &(A2.permutation[0]), A2.order * sizeof(dim_t));

    AB->sym = A2.sym;
  }
  else
  {
    AT->order = A0.order;
    memcpy(&((AT->size)[0]), &(A0.size[0]), A0.order * sizeof(dim_t));
    memcpy(&((AT->offset)[0]), &(A0.offset[0]), A0.order * sizeof(dim_t));
    AT->base = A0.base;
    memcpy(&((AT->permutation)[0]), &(A0.permutation[0]), A0.order * sizeof(dim_t));

    AB->order = A1.order;
    memcpy(&((AB->size)[0]), &(A1.size[0]), A1.order * sizeof(dim_t));
    AB->size[mode] += A2.size[mode];
    memcpy(&((AB->offset)[0]), &(A1.offset[0]), A1.order * sizeof(dim_t));
    AB->base = A1.base;
    memcpy(&((AB->permutation)[0]), &(A1.permutation[0]), A1.order * sizeof(dim_t));

    AB->sym = A1.sym;
  }
	FLA_Adjust_2D_info(AT);
	FLA_Adjust_2D_info(AB);
	return FLA_SUCCESS;
}


//
// --- FLA_Merge_1xmode2() ---------------------------------------------------------
//

FLA_Error FLA_Merge_1xmode2( FLA_Obj AT,
                             FLA_Obj AB,  FLA_Obj *A, dim_t mode )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Merge_1xmode2_check( AT,
                         AB,   A, mode );

  A->order = AT.order;
  memcpy(&((A->size)[0]), &(AT.size[0]), AT.order * sizeof(dim_t));
  (A->size)[mode] += AB.size[mode];
  memcpy(&((A->offset)[0]), &(AT.offset[0]), AT.order * sizeof(dim_t));
  A->base = AT.base;
  memcpy(&((A->permutation)[0]), &(AT.permutation[0]), AT.order * sizeof(dim_t));


  FLA_Adjust_2D_info(A);
  return FLA_SUCCESS;
}

//
// --- FLA_Merge_2powm() ---------------------------------------------------------
//
//NOTE: ONLY FOR TTTTTTTTTT...
                        //C won't allow
                        //FLA_Obj const * const Apart[]
FLA_Error FLA_Merge_2powm(FLA_Obj* Apart[], FLA_Obj* A,
                          dim_t nModes_merge, dim_t merge_modes[nModes_merge])
{
    dim_t i;
    dim_t const num_part = 1 << nModes_merge;
    dim_t stride = num_part / 2;

    A->order = Apart[0]->order;
    memcpy(&((A->size)[0]), &(((Apart[0])->size)[0]), ((Apart[0])->order) * sizeof(dim_t));
    for(i = nModes_merge; i < nModes_merge; i--){
        (A->size)[merge_modes[i]] = ((Apart[0])->size)[merge_modes[i]] + ((Apart[stride])->size)[merge_modes[i]];
        stride /= 2;
    }

    memcpy(&((A->offset)[0]), &(((Apart[0])->offset)[0]), ((Apart[0])->order) * sizeof(dim_t));
    memcpy(&((A->permutation)[0]), &(((Apart[0])->permutation)[0]), ((Apart[0])->order) * sizeof(dim_t));
    for(i = 0; i < nModes_merge; i++){
        (A->offset)[merge_modes[i]] = ((Apart[0])->offset)[merge_modes[i]];
    }
    A->sym = Apart[0]->sym;
    A->base = (Apart[0])->base;
    return FLA_SUCCESS;
}

/*
FLA_Error FLA_Merge_2powm(FLA_Obj const * const Aparts[], FLA_Obj* A,
                          dim_t nModes_repart, dim_t repart_modes[nModes_repart])
{
	dim_t i, j;

	dim_t sizeAparts = 1 << nModes_repart;
	FLA_Obj mergedObjs[sizeAparts / 2];

	memcpy(&(mergedObjs[0]), &(Aparts[0]), sizeAparts * sizeof(FLA_Obj));

	for(i = 0; i < nModes_repart; i++){
		for(j = 0; j < sizeAparts/2; j++){
			FLA_Merge_1xmode2(mergedObjs[j], 
							  mergedObjs[j+1], &(mergedObjs[j]), repart_modes[i]);
		}
		sizeAparts /= 2;
	}
	memcpy(A, &(mergedObjs[0]), sizeof(FLA_Obj));
	return FLA_SUCCESS;
}
*/
