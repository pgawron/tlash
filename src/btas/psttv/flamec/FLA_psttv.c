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

//Note: Only retains symmetry that exists...
//Note: Mode multiplies MUST be INORDER (so that traverse stored pieces correctly).
//(This might could be relaxed since no matter the loop order, we will hit the unique part only once...Not sure...I think we would just have to handle the permutations)
//Step 1: Partition unaltered symmetric groups of A and C first!!
//Step 2: Deal with the symGroup that will be broken
FLA_Error FLA_Psttv( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
    dim_t i;

    //Determine which (if any) symmetric group can safely be repartitioned similarly
    //between A and C.
    TLA_sym symC = C.sym;
    dim_t symGroupToSplit = -1;
    dim_t nModes_part;
    for(i = 0; i < symC.nSymGroups; i++)
        if(symC.symGroupLens[i] > 1){
            symGroupToSplit = i;
            nModes_part = symC.symGroupLens[symGroupToSplit];
            break;
        }

    //No group can be split, meaning mode multiplied in is on own in both tensors.
    //Multiply
	if(symGroupToSplit == -1){
	    FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
	    return FLA_SUCCESS;
	}else{
		//This is the symmetric group to split
	    dim_t symGroupToSplitOffset = TLA_sym_group_mode_offset(symC, symGroupToSplit);

	    dim_t* part_modes;
	    dim_t* sizes;
	    dim_t* repart_sizes;
	    FLA_Side* sides;
	    FLA_Side* repart_sides;

	    dim_t isSingleBlock;

	    FLA_Obj** Apart;
	    FLA_Obj** Cpart;
	    FLA_Obj** Arepart;
	    FLA_Obj** Crepart;

	    dim_t update_region_stride;
	    dim_t update_region;

	    FLA_Obj Apass;
	    FLA_Obj Cpass;

	    //Initialize Views & data for loop
	    dim_t nPart = 1 << nModes_part;
	    dim_t nRepart = 1;
	    for(i = 0; i < nModes_part; i++)
	        nRepart *= 3;

	    //Check if we are dealing with a single block
	    //If so, we get to just multiply in a mode
	    isSingleBlock = TRUE;
	    for(i = 0; i < nModes_part; i++){
	        if(FLA_Obj_dimsize(C, symC.symModes[symGroupToSplitOffset + i]) == 0){
	            return FLA_SUCCESS;
	        }
	        if(FLA_Obj_dimsize(C, symC.symModes[symGroupToSplitOffset + i]) > 1){
	            isSingleBlock = FALSE;
	        }
	    }

	    if(isSingleBlock){
	        FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
	        return FLA_SUCCESS;
	    }

	    part_modes = (dim_t*)FLA_malloc(nModes_part * sizeof(dim_t));
	    sizes = (dim_t*)FLA_malloc(nModes_part * sizeof(dim_t));
	    repart_sizes = (dim_t*)FLA_malloc(nModes_part * sizeof(dim_t));
	    sides = (FLA_Side*)FLA_malloc(nModes_part * sizeof(dim_t));
	    repart_sides = (FLA_Side*)FLA_malloc(nModes_part * sizeof(dim_t));

	    for(i = 0; i < nModes_part; i++){
	    	part_modes[i] = symC.symModes[symGroupToSplitOffset + i];
	        sizes[i] = 0;
	        repart_sizes[i] = 1;
	        sides[i] = FLA_TOP;
	        repart_sides[i] = FLA_BOTTOM;
	    }

	    //Begin loop for general tensor case
	    Apart = (FLA_Obj**)FLA_malloc(nPart * sizeof(FLA_Obj*));
	    Cpart = (FLA_Obj**)FLA_malloc(nPart * sizeof(FLA_Obj*));

	    Arepart = (FLA_Obj**)FLA_malloc(nRepart * sizeof(FLA_Obj*));
	    Crepart = (FLA_Obj**)FLA_malloc(nRepart * sizeof(FLA_Obj*));

	    TLA_create_part_obj(nPart, Apart);
	    TLA_create_part_obj(nPart, Cpart);

	    TLA_create_part_obj(nRepart, Arepart);
	    TLA_create_part_obj(nRepart, Crepart);

	    FLA_Part_2powm(A, Apart,
	                       nModes_part, part_modes,
	                       sizes, sides);

	    FLA_Part_2powm(C, Cpart,
	                       nModes_part, part_modes,
	                       sizes, sides);

	    while(FLA_Obj_dimsize(*(Cpart[0]), part_modes[0]) < FLA_Obj_dimsize(C, part_modes[0])){
	        FLA_Repart_2powm_to_3powm(Apart, Arepart,
	                                      nModes_part, part_modes,
	                                      repart_sizes, repart_sides);
	        FLA_Repart_2powm_to_3powm(Cpart, Crepart,
	                                      nModes_part, part_modes,
	                                      repart_sizes, repart_sides);

			/******************************/
	        update_region_stride = 1;
	        for(i = 1; i < nModes_part; i++){
	            update_region_stride *= 3;
	        }

	        //Symmetric region being partitioned includes
	        //symmetric tensors of order 0->order-1
	        //Must update ALL of them
	        update_region = update_region_stride;
	        for(i = 0; i < nModes_part; i++){
                Apass = *(Arepart[update_region]);
                Cpass = *(Crepart[update_region]);
	            FLA_Psttv(alpha, Apass, mode, beta, B, Cpass);
	            update_region_stride /= 3;
	            update_region += update_region_stride;
	        }
	        /******************************/
	        FLA_Cont_with_3powm_to_2powm(Apart, Arepart,
	                                         nModes_part, part_modes,
	                                         repart_sides);
	        FLA_Cont_with_3powm_to_2powm(Cpart, Crepart,
	                                         nModes_part, part_modes,
                                             repart_sides);
	    }

	    //Tidy up alloc'd data
	    TLA_destroy_part_obj(nPart, Apart);
	    TLA_destroy_part_obj(nPart, Cpart);

	    TLA_destroy_part_obj(nRepart, Arepart);
	    TLA_destroy_part_obj(nRepart, Crepart);


	    FLA_free(part_modes);
	    FLA_free(sizes);
	    FLA_free(repart_sizes);
	    FLA_free(sides);
	    FLA_free(repart_sides);

	    FLA_free(Apart);
	    FLA_free(Cpart);
	    FLA_free(Arepart);
	    FLA_free(Crepart);
	}

	return FLA_SUCCESS;
}

