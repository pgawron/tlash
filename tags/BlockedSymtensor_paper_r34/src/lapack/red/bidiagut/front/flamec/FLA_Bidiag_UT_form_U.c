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

FLA_Error FLA_Bidiag_UT_form_U( FLA_Obj A, FLA_Obj T, FLA_Obj U )
{
	FLA_Uplo uplo;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Bidiag_UT_form_U_check( A, T, U );

	uplo = ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) ?
                                  FLA_UPPER_TRIANGULAR :
                                  FLA_LOWER_TRIANGULAR );

	if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		FLA_QR_UT_form_Q( A, T, U );
		//FLA_Copyr( FLA_LOWER_TRIANGULAR, A, U );
		//FLA_QR_UT_form_Q( U, T, U );
/*
		{
			FLA_Obj W;
			FLA_Set_to_identity( U );
			//FLA_Set( FLA_ZERO, U );
            FLA_Apply_Q_UT_create_workspace( T, U, &W );
			FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
			                A, T, W, U );
            FLA_Obj_free( &W );
		}
*/
	}
	else // if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		FLA_Obj ATL, ATR,
		        ABL, ABR;

		FLA_Obj UTL, UTR,
		        UBL, UBR;

		FLA_Obj TL,  TR;

		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    1, 1, FLA_TR );

		FLA_Part_2x2( U,    &UTL, &UTR,
		                    &UBL, &UBR,    1, 1, FLA_TL );

		FLA_Part_1x2( T,    &TL,  &TR,     1, FLA_RIGHT );

		FLA_Set( FLA_ONE,  UTL );
		FLA_Set( FLA_ZERO, UBL );
		FLA_Set( FLA_ZERO, UTR );

		FLA_QR_UT_form_Q( ABL, TL, UBR );
	}

	return FLA_SUCCESS;
}

