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

FLA_Error FLA_Pythag2( FLA_Obj chi, FLA_Obj psi, FLA_Obj rho )
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( chi );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  buff_chi = FLA_FLOAT_PTR( chi );
			float*  buff_psi = FLA_FLOAT_PTR( psi );
			float*  buff_rho = FLA_FLOAT_PTR( rho );

			FLA_Pythag2_ops( buff_chi,
			                 buff_psi,
			                 buff_rho );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_chi = FLA_DOUBLE_PTR( chi );
			double* buff_psi = FLA_DOUBLE_PTR( psi );
			double* buff_rho = FLA_DOUBLE_PTR( rho );

			FLA_Pythag2_opd( buff_chi,
			                 buff_psi,
			                 buff_rho );

			break;
		}

		case FLA_COMPLEX:
		{
			FLA_Check_error_code( FLA_OBJECT_NOT_REAL );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			FLA_Check_error_code( FLA_OBJECT_NOT_REAL );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Pythag2_ops( float*    chi,
                           float*    psi,
                           float*    rho )
{
	float  zero = bli_s0();
	float  one  = bli_s1();

	float  xabs, yabs;
	float  w, z;
	float  zdivw;

	xabs = fabsf( *chi );
	yabs = fabsf( *psi );
	w    = max( xabs, yabs );
	z    = min( xabs, yabs );

	if ( z == zero )
	{
		*rho = w;
	}
	else
	{
		zdivw = z / w;

		*rho = w * sqrt( one + zdivw * zdivw );
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Pythag2_opd( double*   chi,
                           double*   psi,
                           double*   rho )
{
	double zero = bli_d0();
	double one  = bli_d1();

	double xabs, yabs;
	double w, z;
	double zdivw;

	xabs = fabs( *chi );
	yabs = fabs( *psi );
	w    = max( xabs, yabs );
	z    = min( xabs, yabs );

	if ( z == zero )
	{
		*rho = w;
	}
	else
	{
		zdivw = z / w;

		*rho = w * sqrt( one + zdivw * zdivw );
	}

	return FLA_SUCCESS;
}

