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
                          dim_t nModes_part, const dim_t part_modes[],
                          const dim_t sizes[], const FLA_Side sides[]){

    dim_t i,j,k;
    dim_t num_part = 1 << nModes_part;

    //dim_t order = FLA_Obj_order(A);
    const dim_t* mode_size = A.size;
    const dim_t* mode_offset = A.offset;
    dim_t part_mode_stride = 1;              //mode stride of part array

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

        memcpy(&(((Apart[i])->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
        memcpy(&(((Apart[i])->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
    }

    for(i = 0; i < nModes_part; i++){
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
        part_mode_stride *= 2;
    }

    //All regions need to update the symmetry of the objects
    for(i = 0; i < num_part; i++){
        memcpy(&(((Apart[i])->permutation)[0]), &(A.permutation[0]), FLA_Obj_order(A) * sizeof(dim_t));
        TLA_update_sym_based_offset(A.sym, Apart[i]);
    }

    return FLA_SUCCESS;
}
//
// --- FLA_Part_1xmode2() ----------------------------------------------------------
//

FLA_Error FLA_Part_1xmode2( FLA_Obj A,  FLA_Obj *A1,
                                        FLA_Obj *A2,
                            dim_t mode, dim_t  b,  FLA_Side side )
{
  dim_t split_mode_arr[1];

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Part_1xmode2_check( A,    A1,
                                  A2, mode, b, side );

  // Safeguard: if mb > m, reduce mb to m.
  if ( b > A.size[mode] ) b = A.size[mode];

  // Set mb to be the dimension of A1.
  if ( side == FLA_BOTTOM ) b = A.size[mode] - b;

  //Adjust A1 (order, size, base, & permutation)
  A1->order = A.order;
  memcpy(&((A1->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
  (A1->size)[mode] = b;
  memcpy(&((A1->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
  A1->base = A.base;
  memcpy(&((A1->permutation)[0]), &(A.permutation[0]), A.order * sizeof(dim_t));

  //Adjust A2 (order, size, base, & permutation)
  A2->order = A.order;
  memcpy(&((A2->size)[0]), &(A.size[0]), A.order * sizeof(dim_t));
  (A2->size)[mode] = A.size[mode] - b;
  memcpy(&((A2->offset)[0]), &(A.offset[0]), A.order * sizeof(dim_t));
  (A2->offset)[mode] += b;
  A2->base = A.base;
  memcpy(&((A2->permutation)[0]), &(A.permutation[0]), A.order * sizeof(dim_t));

  //Update symmetries

  split_mode_arr[0] = mode;
  TLA_split_sym_group(A.sym, 1, split_mode_arr, &(A1->sym));
  TLA_split_sym_group(A.sym, 1, split_mode_arr, &(A2->sym));

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
FLA_Error FLA_Repart_2powm_to_3powm( FLA_Obj* Apart[], FLA_Obj* Arepart[],
                                     dim_t nModes_repart,
                                     const dim_t repart_modes[],
                                     const dim_t sizes[],
                                     const FLA_Side sides[] )
{
    dim_t i, j, k;
    dim_t nReparts = 1;
    dim_t* part_base;

    dim_t part_mode_stride;
    dim_t repart_mode_stride;
    for(i = 0; i < nModes_repart; i++)
        nReparts *= 3;

    part_base = (dim_t*)FLA_malloc(nReparts * sizeof(dim_t));
    memset(&(part_base[0]), 0, nReparts * sizeof(dim_t));

    part_mode_stride = 1;
    repart_mode_stride = 1;
    for(i = 0; i < nModes_repart; i++){
        for(j = 0; j < nReparts / (repart_mode_stride * 3); j++){
            //Region with 1 or 2 in current mode (share same base object)
            for(k = repart_mode_stride * (j*3+1); k < (j*3+3)*repart_mode_stride; k++){
                part_base[k] += part_mode_stride;
            }
        }
        repart_mode_stride *= 3;
        part_mode_stride *= 2;
    }

    //Update base info for all repart regions
    for(i = 0; i < nReparts; i++){
        Arepart[i]->order = Apart[part_base[i]]->order;
        Arepart[i]->base = Apart[part_base[i]]->base;
		Arepart[i]->sym = Apart[part_base[i]]->sym;
        memcpy(&((Arepart[i]->permutation)[0]), &((Apart[part_base[i]]->permutation)[0]), (Apart[part_base[i]]->order) * sizeof(dim_t));
        //Very well could be wrong
        memcpy(&((Arepart[i]->size)[0]), &((Apart[part_base[i]]->size)[0]), (Apart[part_base[i]]->order) * sizeof(dim_t));
        memcpy(&((Arepart[i]->offset)[0]), &((Apart[part_base[i]]->offset)[0]), (Apart[part_base[i]]->order) * sizeof(dim_t));
    }

    //Update size and offset arrays
    repart_mode_stride = 1;
    for(i = 0; i < nModes_repart; i++){
        for(j = 0; j < nReparts / (repart_mode_stride * 3); j++){
            for(k = repart_mode_stride * (j*3); k < (j*3+1) * repart_mode_stride; k++){
                (Arepart[k]->size)[repart_modes[i]] = ((Apart[part_base[k]])->size)[repart_modes[i]];
            }
            for(k = repart_mode_stride * (j*3+1); k < (j*3+2) * repart_mode_stride; k++){
                (Arepart[k]->size)[repart_modes[i]] = sizes[i];
            }
            for(k = repart_mode_stride * (j*3+2); k < (j*3+3)*repart_mode_stride; k++){
                (Arepart[k]->offset)[repart_modes[i]] = (Apart[part_base[k]]->offset)[repart_modes[i]] + sizes[i];
                (Arepart[k]->size)[repart_modes[i]] = (Apart[part_base[k]]->size)[repart_modes[i]] - sizes[i];
            }
        }
        repart_mode_stride *= 3;
    }

    //Based on offsets of parts, adjust symmetry
    for(i = 0; i < nReparts; i++){
        TLA_update_sym_based_offset(Arepart[i]->sym, Arepart[i]);
    }

    FLA_free(part_base);
    return FLA_SUCCESS;
}


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
  }
  else
  {
    A0->order = AT.order;
    memcpy(&((A0->size)[0]), &(AT.size[0]), AT.order * sizeof(dim_t));
    memcpy(&((A0->offset)[0]), &(AT.offset[0]), AT.order * sizeof(dim_t));
	memcpy(&((A0->permutation)[0]), &(AT.permutation[0]), AT.order * sizeof(dim_t));
    A0->base = AT.base;

    A0->sym = AT.sym;

    FLA_Part_1xmode2 ( AB,    A1,
                              A2,    mode, b, FLA_TOP );
  }

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
                                        dim_t nModes_cont_with, const dim_t cont_with_modes[],
                                        const FLA_Side sides[]){
    dim_t i,j,k;
    dim_t const num_part = 1 << nModes_cont_with;
    dim_t num_repart = 1;

    dim_t part_mode_stride = 1;
    dim_t repart_mode_stride = 1;

    dim_t* repart_base = (dim_t*)FLA_malloc(num_part * sizeof(dim_t));
    dim_t* repart_update = (dim_t*)FLA_malloc(num_part * sizeof(dim_t));

    for(i = 0; i < nModes_cont_with; i++)
        num_repart *= 3;

    memset(&(repart_base[0]), 0, num_part * sizeof(dim_t));
    memset(&(repart_update[0]), 0, num_part * sizeof(dim_t));

    for(i = 0; i < nModes_cont_with; i++){
        for(j = 0; j < num_part / (part_mode_stride * 2); j++){
            for(k = (j*2)*part_mode_stride; k < (j*2+1)*part_mode_stride; k++){
                repart_update[k] += repart_mode_stride;
            }
            for(k = (j*2+1)*part_mode_stride; k < (j*2+2)*part_mode_stride; k++){
                repart_update[k] += 2*repart_mode_stride;
                repart_base[k] += 2*repart_mode_stride;
            }
        }
        repart_mode_stride *= 3;
        part_mode_stride *= 2;
    }

    for(i = 0; i < num_part; i++){
        Apart[i]->order = Arepart[repart_base[i]]->order;
        Apart[i]->sym = Arepart[repart_base[i]]->sym;
        Apart[i]->base = Arepart[repart_base[i]]->base;
        memcpy(&((Apart[i]->permutation)[0]), &((Arepart[repart_base[i]]->permutation)[0]), (Arepart[repart_base[i]]->order) * sizeof(dim_t));
        memcpy(&((Apart[i]->offset)[0]), &((Arepart[repart_base[i]]->offset)[0]), (Arepart[repart_base[i]]->order) * sizeof(dim_t));
        for(j = 0; j < nModes_cont_with; j++){
            (Apart[i]->size)[cont_with_modes[j]] = (Arepart[repart_base[i]]->size)[cont_with_modes[j]];
        }
    }

    part_mode_stride = 1;
    for(i = 0; i < nModes_cont_with; i++){
        for(j = 0; j < num_part / (part_mode_stride * 2); j++){
            for(k = (2*j)*part_mode_stride; k < (2*j + 1)*part_mode_stride; k++){
                (Apart[k]->size)[cont_with_modes[i]] += (Arepart[repart_update[k]]->size)[cont_with_modes[i]];
            }
        }
        part_mode_stride *= 2;
    }

    FLA_free(repart_base);
    FLA_free(repart_update);

    return FLA_SUCCESS;
}

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

	if (side == FLA_TOP) {
		AT->order = A0.order;
		memcpy(&((AT->size)[0]), &(A0.size[0]), A0.order * sizeof(dim_t));
		AT->size[mode] += A1.size[mode];
		memcpy(&((AT->offset)[0]), &(A0.offset[0]), A0.order * sizeof(dim_t));
		AT->base = A0.base;
		memcpy(&((AT->permutation)[0]), &(A0.permutation[0]),
				A0.order * sizeof(dim_t));

		AB->order = A2.order;
		memcpy(&((AB->size)[0]), &(A2.size[0]), A2.order * sizeof(dim_t));
		memcpy(&((AB->offset)[0]), &(A2.offset[0]), A2.order * sizeof(dim_t));
		AB->base = A2.base;
		memcpy(&((AB->permutation)[0]), &(A2.permutation[0]),
				A2.order * sizeof(dim_t));

		AB->sym = A2.sym;
	} else {
		AT->order = A0.order;
		memcpy(&((AT->size)[0]), &(A0.size[0]), A0.order * sizeof(dim_t));
		memcpy(&((AT->offset)[0]), &(A0.offset[0]), A0.order * sizeof(dim_t));
		AT->base = A0.base;
		memcpy(&((AT->permutation)[0]), &(A0.permutation[0]),
				A0.order * sizeof(dim_t));

		AB->order = A1.order;
		memcpy(&((AB->size)[0]), &(A1.size[0]), A1.order * sizeof(dim_t));
		AB->size[mode] += A2.size[mode];
		memcpy(&((AB->offset)[0]), &(A1.offset[0]), A1.order * sizeof(dim_t));
		AB->base = A1.base;
		memcpy(&((AB->permutation)[0]), &(A1.permutation[0]),
				A1.order * sizeof(dim_t));

		AB->sym = A1.sym;
	}
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

  return FLA_SUCCESS;
}

//
// --- FLA_Merge_2powm() ---------------------------------------------------------
//
//NOTE: ONLY FOR TTTTTTTTTT...
                        //C won't allow
                        //FLA_Obj const * const Apart[]
FLA_Error FLA_Merge_2powm(FLA_Obj* Apart[], FLA_Obj* A,
                          dim_t nModes_merge, const dim_t merge_modes[])
{
    dim_t i;
    dim_t stride = 1;

    A->order = Apart[0]->order;
    memcpy(&((A->size)[0]), &(((Apart[0])->size)[0]), ((Apart[0])->order) * sizeof(dim_t));
    for(i = 0; i < nModes_merge; i++){
        (A->size)[merge_modes[i]] = ((Apart[0])->size)[merge_modes[i]] + ((Apart[stride])->size)[merge_modes[i]];
        stride *= 2;
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
