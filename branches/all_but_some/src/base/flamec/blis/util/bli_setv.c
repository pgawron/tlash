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

void bli_isetv( int n, int* sigma, int* x, int incx )
{
	int*   chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bli_ssetv( int n, float* sigma, float* x, int incx )
{
	float* chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bli_dsetv( int n, double* sigma, double* x, int incx )
{
	double* chi;
	int     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bli_csetv( int n, scomplex* sigma, scomplex* x, int incx )
{
	scomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

void bli_zsetv( int n, dcomplex* sigma, dcomplex* x, int incx )
{
	dcomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

