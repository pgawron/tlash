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

FLA_Error FLA_Apply_G_internal( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A )
{
	FLA_Error r_val = FLA_SUCCESS;

	if      ( side == FLA_LEFT )
	{
		if      ( direct == FLA_FORWARD )
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			//r_val = FLA_Apply_G_lf_opt_var1( G, A );
		}
		else if ( direct == FLA_BACKWARD )
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			//r_val = FLA_Apply_G_lb_opt_var1( G, A );
		}
	}
	else if ( side == FLA_RIGHT )
	{
		if      ( direct == FLA_FORWARD )
		{
			r_val = FLA_Apply_G_rf_opt_var1( G, A );
		}
		else if ( direct == FLA_BACKWARD )
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
			//r_val = FLA_Apply_G_rb_opt_var1( G, A );
		}
	}

	return r_val;
}

