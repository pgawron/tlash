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

FLA_Error FLA_Apply_QUD_UT_inc_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                         FLA_Obj T, FLA_Obj W,
                                                    FLA_Obj R,
                                         FLA_Obj U, FLA_Obj C,
                                         FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Apply_QUD_UT_inc_internal_check( side, trans, direct, storev, T, W, R, U, C, V, D, cntl );

	if      ( side == FLA_LEFT )
	{
		if      ( trans == FLA_NO_TRANSPOSE )
		{
			if      ( direct == FLA_FORWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
			else if ( direct == FLA_BACKWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
		}
		else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
		{
			if      ( direct == FLA_FORWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
					r_val = FLA_Apply_QUD_UT_inc_lhfc( T, W, R, U, C, V, D, cntl );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
			else if ( direct == FLA_BACKWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
		}
	}
	else if ( side == FLA_RIGHT )
	{
		if      ( trans == FLA_NO_TRANSPOSE )
		{
			if      ( direct == FLA_FORWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
			else if ( direct == FLA_BACKWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
		}
		else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
		{
			if      ( direct == FLA_FORWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
			else if ( direct == FLA_BACKWARD )
			{
				if      ( storev == FLA_COLUMNWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				else if ( storev == FLA_ROWWISE )
  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			}
		}
	}

	return r_val;
}
