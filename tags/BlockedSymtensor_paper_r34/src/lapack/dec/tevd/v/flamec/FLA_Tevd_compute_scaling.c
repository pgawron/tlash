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

FLA_Error FLA_Tevd_compute_scaling_ops( int       m_A,
                                        float*    buff_d, int inc_d, 
                                        float*    buff_e, int inc_e,
                                        float*    sigma )
{
	float  one   = bli_s1();
	float  three = 3.0F;
	float  norm;
	float  eps2;
	float  safmin;
	float  safmax;
	float  ssfmin;
	float  ssfmax;

	// Query some constants.
	eps2   = FLA_Mach_params_ops( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_ops( FLA_MACH_SFMIN );
	safmax = one / safmin;

	// Compute the acceptable range for the 1-norm;
	ssfmax = sqrt( safmax ) / three;
	ssfmin = sqrt( safmin ) / eps2;

	// Compute the 1-norm of the tridiagonal matrix.
	FLA_Norm1_tridiag_ops( m_A,
	                       buff_d, inc_d,
	                       buff_e, inc_e,
	                       &norm );

	// Compute sigma accordingly if norm is outside of the range.
	if ( norm > ssfmax )
	{
		*sigma = ssfmax / norm;
	}
	else if ( norm < ssfmin )
	{
		*sigma = ssfmin / norm;
	}
	else
	{
		*sigma = one;
	}

	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_compute_scaling_opd( int       m_A,
                                        double*   buff_d, int inc_d, 
                                        double*   buff_e, int inc_e,
                                        double*   sigma )
{
	double one   = bli_d1();
	double three = 3.0;
	double norm;
	double eps2;
	double safmin;
	double safmax;
	double ssfmin;
	double ssfmax;

	// Query some constants.
	eps2   = FLA_Mach_params_opd( FLA_MACH_EPS2 );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );
	safmax = one / safmin;

	// Compute the acceptable range for the 1-norm;
	ssfmax = sqrt( safmax ) / three;
	ssfmin = sqrt( safmin ) / eps2;

	// Compute the 1-norm of the tridiagonal matrix.
	FLA_Norm1_tridiag_opd( m_A,
	                       buff_d, inc_d,
	                       buff_e, inc_e,
	                       &norm );

	// Compute sigma accordingly if norm is outside of the range.
	if ( norm > ssfmax )
	{
		*sigma = ssfmax / norm;
	}
	else if ( norm < ssfmin )
	{
		*sigma = ssfmin / norm;
	}
	else
	{
		*sigma = one;
	}

	return FLA_SUCCESS;
}

