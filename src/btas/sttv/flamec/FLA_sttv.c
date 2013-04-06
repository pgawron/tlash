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
                                      /**/ /**/
                                          &A1,
                                      AB, &A2, repart_mode, b, FLA_BOTTOM);
        FLA_Repart_1xmode2_to_1xmode3(CT, &C0,
                                      /**/ /**/
                                          &C1,
                                      CB, &C2, repart_mode, b, FLA_BOTTOM);

        /*********************************/
        FLA_Repart_sym_group(alpha, A1, mode, beta, B, C1, symGroupsToPartition, symGroupsPartitioned, symModeToPartition - 1, loopIndex);
        /*********************************/

        FLA_Cont_with_1xmode3_to_1xmode2( &CT, C0,
                                               C1,
                                        /********/
                                          &CB, C2, repart_mode, FLA_TOP);
        FLA_Cont_with_1xmode3_to_1xmode2( &AT, A0,
                                               A1,
                                        /********/
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

//Note: Only retains symmetry that exists...
//Note: Mode multiplies MUST be INORDER (so that traverse stored pieces correctly). (This might could be relaxed since no matter the loop order, we will hit the unique part only once...Not sure...I think we would just have to handle the permutations)
FLA_Error FLA_Sttv( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C )
{
	TLA_sym symC = C.sym;
	dim_t i;
	dim_t symModePos = TLA_sym_pos_of_mode(C.sym, mode);
	dim_t symGroupToIgnore;
	for(i = 0; i < symC.nSymGroups; i++){
	    if(symC.symGroupLens[i] == 1 && symC.symModes[symModePos] == mode){
	        symGroupToIgnore = i;
	        break;
	    }
	}

	FLA_Sttv_helper(alpha, A, mode, beta, B, C, symC, symGroupToIgnore, 0);

	return FLA_SUCCESS;
}

