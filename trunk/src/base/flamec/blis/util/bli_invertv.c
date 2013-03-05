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

#include "blis.h"

void bli_sinvertv( conj_t conj, int n, float* x, int incx )
{
	float  one = 1.0F;
	float* chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = one / *chi;
	}
}

void bli_dinvertv( conj_t conj, int n, double* x, int incx )
{
	double  one = 1.0;
	double* chi;
	int     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = one / *chi;
	}
}

void bli_cinvertv( conj_t conj, int n, scomplex* x, int incx )
{
	float     one = 1.0F;
	float     temp;
	float     conjsign;
	scomplex* chi;
	int       i;

	if ( bli_is_conj( conj ) ) conjsign =  one;
	else                       conjsign = -one;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

	    temp = one / ( chi->real * chi->real +
	                   chi->imag * chi->imag );
	    chi->real = chi->real *            temp;
	    chi->imag = chi->imag * conjsign * temp;
	}
}

void bli_zinvertv( conj_t conj, int n, dcomplex* x, int incx )
{
	double    one = 1.0;
	double    temp;
	double    conjsign;
	dcomplex* chi;
	int       i;

	if ( bli_is_conj( conj ) ) conjsign =  one;
	else                       conjsign = -one;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

	    temp = one / ( chi->real * chi->real +
	                   chi->imag * chi->imag );
	    chi->real = chi->real *            temp;
	    chi->imag = chi->imag * conjsign * temp;
	}
}

