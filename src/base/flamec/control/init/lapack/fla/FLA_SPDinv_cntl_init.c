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

extern fla_chol_t*  fla_chol_cntl;
extern fla_trinv_t* fla_trinv_cntl;
extern fla_ttmm_t*  fla_ttmm_cntl;

fla_spdinv_t*       fla_spdinv_cntl;
fla_blocksize_t*    fla_spdinv_size_cutoff;

void FLA_SPDinv_cntl_init()
{
	// Rather than embed a blocksize, we store the cutoff matrix size for
	// switching from external routines to internal FLAME variants.
	fla_spdinv_size_cutoff = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	// Initialize a control tree node that calls the top-level Cholesky
	// factorization, Trinagular inversion, and Triangular-transpose matrix
	// multiply control trees. 
	fla_spdinv_cntl        = FLA_Cntl_spdinv_obj_create( FLA_FLAT,
	                                                     FLA_BLOCKED_VARIANT1, 
	                                                     fla_spdinv_size_cutoff,
	                                                     fla_chol_cntl,
	                                                     fla_trinv_cntl,
	                                                     fla_ttmm_cntl );
}

void FLA_SPDinv_cntl_finalize()
{
	FLA_Cntl_obj_free( fla_spdinv_cntl );

	FLA_Blocksize_free( fla_spdinv_size_cutoff );
}

