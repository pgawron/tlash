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

FLA_Error FLA_Bidiag_UT_form_V( FLA_Obj A, FLA_Obj S, FLA_Obj V )
{
	FLA_Uplo uplo;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Bidiag_UT_form_V_check( A, S, V );

	uplo = ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) ?
                                  FLA_UPPER_TRIANGULAR :
                                  FLA_LOWER_TRIANGULAR );

	if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		FLA_Obj ATL, ATR,
		        ABL, ABR;

		FLA_Obj VTL, VTR,
		        VBL, VBR;

		FLA_Obj SL,  SR;

		dim_t   n_A = FLA_Obj_width( A );

		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    n_A-1, n_A-1, FLA_TR );

		FLA_Part_2x2( V,    &VTL, &VTR,
		                    &VBL, &VBR,    1, 1, FLA_TL );

		FLA_Part_1x2( S,    &SL,  &SR,     1, FLA_RIGHT );

		FLA_Set( FLA_ONE,  VTL );
		FLA_Set( FLA_ZERO, VBL );
		FLA_Set( FLA_ZERO, VTR );

		//FLA_LQ_UT_form_Q( ATR, SL, VBR );
		FLA_Copyrt( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, ATR, VBR );
		FLA_QR_UT_form_Q( VBR, SL, VBR );
	}
	else // if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		//FLA_LQ_UT_form_Q( A, S, V );
		FLA_Copyrt( FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, A, V );
		FLA_QR_UT_form_Q( V, S, V );
	}

	return FLA_SUCCESS;
}

