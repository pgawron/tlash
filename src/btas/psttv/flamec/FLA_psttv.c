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

/*
FLA_Error FLA_Repart_sym_group(FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t symGroupsToPartition[], dim_t symGroupsPartitioned, dim_t symModeToPartition, dim_t endIndex){
    if(symModeToPartition == -1){
       // return FLA_Sttv_helper(alpha, A, mode, beta, B, C, symGroupsToPartition, symGroupsPartitioned + 1);
    }

    FLA_Obj AT, AB;
    FLA_Obj A0, A1, A2;

    FLA_Obj CT, CB;
    FLA_Obj C0, C1, C2;

    dim_t repart_mode = TLA_sym_group_mode_offset(C.sym, symGroupsToPartition[symGroupsPartitioned]) + symModeToPartition;
    FLA_Part_1xmode2(A, &AT,
                        &AB, repart_mode, 0, FLA_TOP);
    FLA_Part_1xmode2(C, &CT,
                        &CB, repart_mode, 0, FLA_TOP);
    dim_t loopIndex = 0;
    while(FLA_Obj_dimsize(CT, repart_mode) <= endIndex){
        dim_t b = 1;
        FLA_Repart_1xmode2_to_1xmode3(AT, &A0,
                                      //
                                          &A1,
                                      AB, &A2, repart_mode, b, FLA_BOTTOM);
        FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
                                      //
                                          &C1,
                                      CB, &C2, repart_mode, b, FLA_BOTTOM);

        //
        FLA_Repart_sym_group(alpha, A1, mode, beta, B, C1, symGroupsToPartition, symGroupsPartitioned, symModeToPartition - 1, loopIndex);
        //

        FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
                                               C1,
                                        //
                                          &CB, C2, repart_mode, FLA_TOP);
        FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
                                               A1,
                                        //
                                          &AB, A2, repart_mode, FLA_TOP);
        loopIndex++;
    }
    return FLA_SUCCESS;
}

FLA_Error FLA_Sttv_helper(FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, TLA_sym origSym, dim_t origSymGroupIgnore, dim_t origSymGroupsPartitioned){
    if(origSymGroupsPartitioned == origSym.nSymGroups - 1){
        //Perform the multiply into the view
        FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
        return FLA_SUCCESS;
    }
    //We are just repartitioning symGroups
    dim_t mode_in_repart_symGroup = origSym.symModes[TLA_sym_group_mode_offset(origSym, origSymGroupsPartitioned)];
    dim_t cur_symGroup_repart = TLA_sym_group_of_mode(C.sym, mode_in_repart_symGroup);
    dim_t cur_mode_repart = C.sym.symModes[TLA_sym_group_mode_offset(C.sym, cur_symGroup_repart) + C.sym.symGroupLens[cur_symGroup_repart] - 1];

    //dim_t endIndex = FLA_Obj_dimsize(C, sym_group_mode_offset + symGroupModePartition);

    //return FLA_Repart_sym_group(alpha, A, mode, beta, B, C, symGroupsToPartition, symGroupsPartitioned, symGroupModePartition, endIndex);
    return FLA_SUCCESS;
}
*/

void print_matlab_tensor(const char* varName, FLA_Obj A){
    dim_t i;
    printf("%s tensor\n", varName);
    printf("%s = tensor([", varName);
    FLA_Obj_print_tensor(A);
    printf("],[");
    for(i = 0; i < FLA_Obj_order(A); i++)
        printf("%d ", FLA_Obj_dimsize(((FLA_Obj*)(FLA_Obj_base_buffer(A)))[0],i) * FLA_Obj_dimsize(A,i));
    printf("]);\n\n");
}
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

    dim_t isSingleBlock = TRUE;
    for(i = 0; i < FLA_Obj_order(A); i++){
        if(FLA_Obj_dimsize(A, i) > 1){
            isSingleBlock = FALSE;
            break;
        }

    }
	if(symGroupToSplit == -1 || isSingleBlock){
	    printf("Calling ttm\n");
	    printf("-----------\n\n");
	    print_matlab_tensor("A", A);
	    print_matlab_tensor("B", B);
	    print_matlab_tensor("preC", C);
	    FLA_Ttm_single_mode(alpha, A, mode, beta, B, C);
	    print_matlab_tensor("C", C);
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
        for(i = 0; i < (1 << nModes_part); i++){
            printf("Cpart[%d]", i);
            print_array("", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
            print_array("  size", (Cpart[i])->order, &(((Cpart[i])->size)[0]));
            print_array("  offset", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
            //printf("  nSymGroups: %d\n", (Arepart[i]->sym).nSymGroups);
            //print_array("  symGroupLens", (Arepart[i]->sym).nSymGroups, (Arepart[i]->sym).symGroupLens);
            //print_array("  symModes", (Arepart[i]->sym).order, (Arepart[i]->sym).symModes);
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
            dim_t n_repart = 1;
            for(i = 0; i < nModes_part; i++)
                n_repart *= 3;
            for(i = 0; i < n_repart; i++){
                printf("Crepart[%d]", i);
                print_array("", (Crepart[i])->order, &(((Crepart[i])->offset)[0]));
                print_array("  size", (Crepart[i])->order, &(((Crepart[i])->size)[0]));
                print_array("  offset", (Crepart[i])->order, &(((Crepart[i])->offset)[0]));
                //printf("  nSymGroups: %d\n", (Arepart[i]->sym).nSymGroups);
                //print_array("  symGroupLens", (Arepart[i]->sym).nSymGroups, (Arepart[i]->sym).symGroupLens);
                //print_array("  symModes", (Arepart[i]->sym).order, (Arepart[i]->sym).symModes);
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
            for(i = 0; i < (1 << nModes_part); i++){
                printf("Cpart[%d]", i);
                print_array("", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
                print_array("  size", (Cpart[i])->order, &(((Cpart[i])->size)[0]));
                print_array("  offset", (Cpart[i])->order, &(((Cpart[i])->offset)[0]));
                //printf("  nSymGroups: %d\n", (Arepart[i]->sym).nSymGroups);
                //print_array("  symGroupLens", (Arepart[i]->sym).nSymGroups, (Arepart[i]->sym).symGroupLens);
                //print_array("  symModes", (Arepart[i]->sym).order, (Arepart[i]->sym).symModes);
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

