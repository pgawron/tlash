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

dim_t FLA_get_unique_info( TLA_sym sym, const dim_t index[], dim_t* sortedIndex, dim_t* permutation, dim_t* ipermutation)
{
	dim_t i, j;

	FLA_Paired_Sort index_pairs[FLA_MAX_ORDER];
	dim_t nSymGroups = sym.nSymGroups;
	dim_t* symGroupLens = &(sym.symGroupLens[0]);
	dim_t* symModes = &(sym.symModes[0]);
	dim_t orderedSymModes[FLA_MAX_ORDER];
	dim_t modeOffset = 0;
	dim_t uniqueIndex = TRUE;

	for(i = 0; i < nSymGroups; i++){
		for(j = 0; j < symGroupLens[i]; j++){
			orderedSymModes[j+modeOffset] = symModes[j+modeOffset];
		}
		qsort(&(orderedSymModes[modeOffset]), symGroupLens[i], sizeof(dim_t), compare_dim_t);

		for(j = 0; j < symGroupLens[i]; j++){
			index_pairs[j].index = orderedSymModes[j+modeOffset];
			index_pairs[j].val = index[orderedSymModes[j+modeOffset]];
		}
		qsort(index_pairs, symGroupLens[i], sizeof(FLA_Paired_Sort), compare_pairwise_sort);


		for(j = 0; j < symGroupLens[i]; j++){
			permutation[orderedSymModes[j+modeOffset]] = index_pairs[j].index;
			sortedIndex[orderedSymModes[j+modeOffset]] = index_pairs[j].val;
		}

		//Check if index is unique or not
		if (symGroupLens[i] > 1) {
			for (j = 0; j < symGroupLens[i] - 1; j++) {
				if (index[symModes[j + modeOffset]]
						> index[symModes[j + modeOffset + 1]]) {
					uniqueIndex = FALSE;
					break;
				}
			}
		}

		modeOffset += symGroupLens[i];
	}

	//Set the correct ipermutation
	memcpy(&(ipermutation[0]), &(permutation[0]), sym.order * sizeof(dim_t));
	modeOffset = 0;
	for(i = 0; i < nSymGroups; i++){
		for(j = 0; j < symGroupLens[i]; j++){
			ipermutation[permutation[symModes[j+modeOffset]]] = orderedSymModes[j+modeOffset];
		}
		modeOffset += symGroupLens[i];
	}

  return uniqueIndex;
}

