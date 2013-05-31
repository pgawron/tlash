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

FLA_Error FLA_SPDinv_internal( FLA_Uplo uplo, FLA_Obj A, fla_spdinv_t* cntl )
{
	FLA_Error r_val;
	FLA_Error e_val;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_SPDinv_internal_check( uplo, A, cntl );

	r_val = FLA_Chol_internal( uplo, A,
	                           FLA_Cntl_sub_chol( cntl ) );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_chol_failure( r_val );
		FLA_Check_error_code( e_val );
	}

	FLA_Trinv_internal( uplo, FLA_NONUNIT_DIAG, A, 
	                    FLA_Cntl_sub_trinv( cntl ) );

	FLA_Ttmm_internal( uplo, A, 
	                   FLA_Cntl_sub_ttmm( cntl ) );

	return FLA_SUCCESS;
}
