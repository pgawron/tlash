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

void bli_sconjv( int m, float* x, int incx )
{
	return;
}

void bli_dconjv( int m, double* x, int incx )
{
	return;
}

void bli_cconjv( int m, scomplex* x, int incx )
{
	float  m1        = bli_sm1();
	float* x_conj    = ( float* ) x + 1;
	int    incx_conj = 2 * incx;

	bli_sscal( m,
	           &m1,
	           x_conj, incx_conj );
}

void bli_zconjv( int m, dcomplex* x, int incx )
{
	double  m1        = bli_dm1();
	double* x_conj    = ( double* ) x + 1;
	int     incx_conj = 2 * incx;

	bli_dscal( m,
	           &m1,
	           x_conj, incx_conj );
}

