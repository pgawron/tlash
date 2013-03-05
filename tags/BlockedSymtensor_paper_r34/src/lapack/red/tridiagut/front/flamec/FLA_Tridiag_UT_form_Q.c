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

FLA_Error FLA_Tridiag_UT_form_Q( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T )
{
	FLA_Error r_val = FLA_SUCCESS;
	FLA_Obj   ATL, ATR,
	          ABL, ABR;
	FLA_Obj   TL,  TR;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Tridiag_UT_form_Q_check( uplo, A, T );

	// Shift the Householder vectors one row/column towards the diagonal.
	FLA_Tridiag_UT_shift_U( uplo, A );

	FLA_Part_2x2( A,    &ATL, &ATR,
	                    &ABL, &ABR,    1, 1, FLA_TL );
	FLA_Part_1x2( T,    &TL, &TR,     1, FLA_RIGHT );

	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		FLA_QR_UT_form_Q( ABR, TL, ABR );
	}
	else // if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		//FLA_LQ_UT_form_Q( A, T, A );
	}

	return r_val;
}

