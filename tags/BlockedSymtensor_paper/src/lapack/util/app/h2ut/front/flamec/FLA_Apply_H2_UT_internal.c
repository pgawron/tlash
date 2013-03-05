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

FLA_Error FLA_Apply_H2_UT_internal( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
	FLA_Error r_val = FLA_SUCCESS;

	if      ( side == FLA_LEFT )
	{
		//r_val = FLA_Apply_H2_UT_l_unb_var1( tau, u2, a1, A2 );
		r_val = FLA_Apply_H2_UT_l_opt_var1( tau, u2, a1, A2 );
	}
	else if ( side == FLA_RIGHT )
	{
		//r_val = FLA_Apply_H2_UT_r_unb_var1( tau, u2, a1, A2 );
		r_val = FLA_Apply_H2_UT_r_opt_var1( tau, u2, a1, A2 );
	}

	return r_val;
}

