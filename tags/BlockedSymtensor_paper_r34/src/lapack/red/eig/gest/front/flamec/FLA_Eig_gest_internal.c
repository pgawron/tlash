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

extern fla_eig_gest_t* flash_eig_gest_cntl;
extern fla_eig_gest_t* fla_eig_gest_ix_cntl_leaf;
extern fla_eig_gest_t* fla_eig_gest_nx_cntl_leaf;

FLA_Error FLA_Eig_gest_internal( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Eig_gest_internal_check( inv, uplo, A, Y, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Eig_gest_internal( inv,
		                               uplo,
		                               *FLASH_OBJ_PTR_AT( A ),
		                               *FLASH_OBJ_PTR_AT( Y ),
		                               *FLASH_OBJ_PTR_AT( B ),
		                               flash_eig_gest_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Eig_gest( inv, uplo, A, Y, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			if ( inv == FLA_INVERSE )
  			  cntl = fla_eig_gest_ix_cntl_leaf;
			else
			  cntl = fla_eig_gest_nx_cntl_leaf;
		}

		// Parameter combinations
		if      ( inv == FLA_INVERSE )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				r_val = FLA_Eig_gest_il( A, Y, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				r_val = FLA_Eig_gest_iu( A, Y, B, cntl );
			}
		}
		else if ( inv == FLA_NO_INVERSE )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				r_val = FLA_Eig_gest_nl( A, Y, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				r_val = FLA_Eig_gest_nu( A, Y, B, cntl );
			}
		}
	}

	return r_val;
}
