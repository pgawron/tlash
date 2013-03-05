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
//TODO: Quadrant not yet generalized.  For now only handles TTTTTT... 
FLA_Error FLA_Part_2powm_helper(FLA_Obj A, FLA_Obj Apart[], dim_t sizes[], FLA_Side sides[], dim_t stride[], dim_t mode){
	dim_t order = FLA_Obj_order(A);
	dim_t createdObj0, createdObj1;

	/* TODO: ERROR CHECKS NEEDED */

	dim_t baseOffset0[order];
	dim_t baseOffset1[order];
	memset(&(baseOffset0[0]), 0, order * sizeof(dim_t));
	memset(&(baseOffset1[0]), 0, order * sizeof(dim_t));

	memcpy(&(baseOffset0[mode+1]), &(A.offset[mode+1]), (order - (mode + 1)) * sizeof(dim_t));	
	memcpy(&(baseOffset1[mode+1]), &(A.offset[mode+1]), (order - (mode + 1)) * sizeof(dim_t));

	baseOffset0[mode] = 0;
	baseOffset1[mode] = 1;

	FLA_TIndex_to_LinIndex(order, stride, baseOffset0, &createdObj0);
	FLA_TIndex_to_LinIndex(order, stride, baseOffset1, &createdObj1);

	FLA_Part_1xmode2(A, &(Apart[createdObj0]),
						&(Apart[createdObj1]),
						mode, sizes[mode], sides[mode]);

	if(mode > 0){
		FLA_Part_2powm_helper(Apart[createdObj0], Apart, sizes, sides, stride, mode - 1);
		FLA_Part_2powm_helper(Apart[createdObj1], Apart, sizes, sides, stride, mode - 1);
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Part_2powm( FLA_Obj A,  FLA_Obj Apart[], dim_t sizes[], FLA_Side sides[])
{
	dim_t order = FLA_Obj_order(A);
	dim_t* stride = FLA_Obj_stride(A);
	FLA_Part_2powm_helper(A, Apart, sizes, sides, stride, order - 1);

	FLA_free(stride);

	return FLA_SUCCESS;
}


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

  FLA_Adjust_2D_info(A1);
  FLA_Adjust_2D_info(A2);

  return FLA_SUCCESS;
}

//
// --- FLA_Repart_2powm_to_3powm() -------------------------------------------------
//
//TODO: See if a hierarchical method is better
//1. Setup hierarchy to turn object into 2x1 in mode working with
//2. Repartition 2x1 to 3x1 in that mode
//3. Repeat for all modes
FLA_Error FLA_Repart_2powm_to_3powm_helper( dim_t order, dim_t sizeApart[], dim_t strideApart[], FLA_Obj Apart[], 
												 			  dim_t strideArepart[], FLA_Obj Arepart[],
										 		 			  dim_t sizes[order], FLA_Side sides[order])
{
	dim_t i;
	FLA_Repart_2xmodeBlank_to_3xmodeBlank( order, sizeApart, strideApart, Apart,
												  strideArepart, Arepart,
										   0, sizes[0], sides[0] );
	for(i = 1; i < order; i++){
		dim_t newSizeApart[order];
		memcpy(&(newSizeApart[0]), &(sizeApart[0]), order * sizeof(dim_t));
		newSizeApart[i] = 3;
		FLA_Repart_2xmodeBlank_to_3xmodeBlank( order, newSizeApart, strideArepart, Arepart,
													  strideArepart, Arepart,
											   i, sizes[i], sides[i] );
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Repart_2powm_to_3powm( dim_t order, FLA_Obj Apart[],  FLA_Obj Arepart[],
                                 	 dim_t sizes[], FLA_Side sides[])
{
	dim_t sizeApart[order];
	dim_t i;
	for(i = 0; i < order;i++){
		sizeApart[i] = 2;
	}

	dim_t strideApart[order];
	dim_t strideArepart[order];
	strideApart[0] = 1;
	strideArepart[0] = 1;

	for(i = 1; i < order; i++){
		strideApart[i] = strideApart[i-1]*sizeApart[i-1];
		strideArepart[i] = strideArepart[i-1]*3;
	}

	return FLA_Repart_2powm_to_3powm_helper( order, sizeApart, strideApart, Apart, 
													strideArepart, Arepart, 
													sizes, sides);
}

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
  }
  else
  {
    A0->order = AT.order;
    memcpy(&((A0->size)[0]), &(AT.size[0]), AT.order * sizeof(dim_t));
    memcpy(&((A0->offset)[0]), &(AT.offset[0]), AT.order * sizeof(dim_t));
	memcpy(&((A0->permutation)[0]), &(AT.permutation[0]), AT.order * sizeof(dim_t));
    A0->base = AT.base;

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
//1. Setup hierarchy to turn object into 2x1 in mode working with
//2. Repartition 2x1 to 3x1 in that mode
//3. Repeat for all modes
FLA_Error FLA_Cont_with_3powm_to_2powm_helper( dim_t order, dim_t strideApart[], FLA_Obj Apart[], 
												 			  dim_t sizeArepart[], dim_t strideArepart[], FLA_Obj Arepart[],
										 		 			  dim_t mode, dim_t sizes[order], FLA_Side sides[order])
{
	dim_t i;
	FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( order, strideApart, Apart,
													 sizeArepart, strideArepart, Arepart,
										   0, sizes[0], sides[0] );
	for(i = 1; i < order; i++){
		dim_t newSizeArepart[order];
		memcpy(&(newSizeArepart[0]), &(sizeArepart[0]), order * sizeof(dim_t));
		newSizeArepart[mode] = 2;
		FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( order, strideApart, Apart,
													  	 newSizeArepart, strideApart, Apart,
											   	  i, sizes[i], sides[i] );
	}
	return FLA_SUCCESS;
}

FLA_Error FLA_Cont_with_3xmodeBlank_to_2xmodeBlank( dim_t order, dim_t strideApart[], FLA_Obj Apart[], 
												 			  dim_t sizeArepart[], dim_t strideArepart[], FLA_Obj Arepart[],
										 		 			  dim_t mode, dim_t size, FLA_Side side)
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

FLA_Error FLA_Merge_2powm(FLA_Obj Aparts[], FLA_Obj* A)
{
	dim_t i, j;
	dim_t order = FLA_Obj_order(*A);
	dim_t sizeAparts = 1 << FLA_Obj_order(*A);
	FLA_Obj mergedObjs[sizeAparts / 2];

	memcpy(&(mergedObjs[0]), &(Aparts[0]), sizeAparts * sizeof(FLA_Obj));

	for(i = 0; i < order; i++){
		for(j = 0; j < sizeAparts/2; j++){
			FLA_Merge_1xmode2(mergedObjs[j], 
							  mergedObjs[j+1], &(mergedObjs[j]), i);
		}
		sizeAparts /= 2;
	}
	memcpy(A, &(mergedObjs[0]), sizeof(FLA_Obj));
	return FLA_SUCCESS;
}
