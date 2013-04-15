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
//Note: Mode multiplies MUST be INORDER (so that traverse stored pieces correctly). (This might could be relaxed since no matter the loop order, we will hit the unique part only once...Not sure...I think we would just have to handle the permutations)
FLA_Error FLA_Psttv( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
    dim_t i;

    TLA_sym symC = C.sym;
    dim_t symGroupToSplit = -1;
    dim_t nModes_part;
    for(i = 0; i < symC.nSymGroups; i++)
        if(symC.symGroupLens[i] > 1){
            symGroupToSplit = i;
            nModes_part = symC.symGroupLens[symGroupToSplit];
            break;
        }

	if(symGroupToSplit == -1){
	    printf("Calling ttm\n");
	    printf("-----------\n\n");
	    FLA_Obj_print_matlab("A", A);
	    FLA_Obj_print_matlab("B", B);
	    FLA_Obj_print_matlab("preC", C);
	    FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
	    FLA_Obj_print_matlab("C", C);
	    return FLA_SUCCESS;
	}else{
	    dim_t symGroupToSplitOffset = TLA_sym_group_mode_offset(symC, symGroupToSplit);

	    dim_t nPart = 1 << nModes_part;
	    dim_t nRepart = 1;
	    for(i = 0; i < nModes_part; i++)
	        nRepart *= 3;

	    dim_t part_modes[nModes_part];
	    dim_t sizes[nModes_part];
	    dim_t repart_sizes[nModes_part];
	    FLA_Side sides[nModes_part];
	    FLA_Side repart_sides[nModes_part];
	    for(i = 0; i < nModes_part; i++){
	        part_modes[i] = symC.symModes[symGroupToSplitOffset + i];
	        sizes[i] = 0;
	        repart_sizes[i] = 1;
	        sides[i] = FLA_TOP;
	        repart_sides[i] = FLA_BOTTOM;
	    }


	    dim_t isSingleBlock = TRUE;
	    for(i = 0; i < nModes_part; i++){
	        if(FLA_Obj_dimsize(C,part_modes[i]) == 0){
	            return FLA_SUCCESS;
	        }
	        if(FLA_Obj_dimsize(C, part_modes[i]) > 1){
	            isSingleBlock = FALSE;
	        }
	    }

	    if(isSingleBlock){
	        printf("Calling ttm\n");
	        printf("-----------\n\n");
	        FLA_Obj_print_matlab("A", A);
	        FLA_Obj_print_matlab("B", B);
	        FLA_Obj_print_matlab("preC", C);
	        FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
	        FLA_Obj_print_matlab("C", C);
	        return FLA_SUCCESS;
	    }

	    FLA_Obj* Apart[nPart];
	    FLA_Obj* Cpart[nPart];

	    FLA_Obj* Arepart[nRepart];
	    FLA_Obj* Crepart[nRepart];

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

        printf("---------------\n");
        printf("After Part\n");
        printf("---------------\n");
        for(i = 0; i < nPart; i++){
            printf("Cpart[%d]", i);
            print_array("", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
            print_array("  size", (Cpart[i])->order, &(((Cpart[i])->size)[0]));
            print_array("  offset", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
            printf("\n");
        }

	    while(FLA_Obj_dimsize(*(Cpart[0]), part_modes[0]) < FLA_Obj_dimsize(C, part_modes[0])){
	        FLA_Repart_2powm_to_3powm(Apart, Arepart,
	                                      nModes_part, part_modes,
	                                      repart_sizes, repart_sides);
	        FLA_Repart_2powm_to_3powm(Cpart, Crepart,
	                                      nModes_part, part_modes,
	                                      repart_sizes, repart_sides);
	        printf("---------------\n");
            printf("After repart\n");
            printf("---------------\n");

            for(i = 0; i < nRepart; i++){
                printf("Crepart[%d]", i);
                print_array("", (Crepart[i])->order, &(((Crepart[i])->offset)[0]));
                print_array("  size", (Crepart[i])->order, &(((Crepart[i])->size)[0]));
                print_array("  offset", (Crepart[i])->order, &(((Crepart[i])->offset)[0]));
                printf("\n");
            }
	        /******************************/
	        dim_t update_region_stride = 1;
	        for(i = 1; i < nModes_part; i++){
	            update_region_stride *= 3;
	        }
	        dim_t update_region = update_region_stride;
	        for(i = 0; i < nModes_part; i++){
	            printf("Recurring on region: %d\n", update_region);
	            printf("-------------------\n");
	            FLA_Psttv(alpha, *(Arepart[update_region]), mode, beta, B, *(Crepart[update_region]));
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
	        printf("---------------\n");
	        printf("After cont with\n");
            printf("---------------\n");
            for(i = 0; i < nPart; i++){
                printf("Cpart[%d]", i);
                print_array("", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
                print_array("  size", (Cpart[i])->order, &(((Cpart[i])->size)[0]));
                print_array("  offset", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
                printf("\n");
            }
	    }

	    TLA_destroy_part_obj(nPart, Apart);
	    TLA_destroy_part_obj(nPart, Cpart);

	    TLA_destroy_part_obj(nRepart, Arepart);
	    TLA_destroy_part_obj(nRepart, Crepart);

	}

	return FLA_SUCCESS;
}

