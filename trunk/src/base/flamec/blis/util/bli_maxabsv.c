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

void bli_smaxabsv( int n, float* x, int incx, float* maxabs )
{
	float*    chi;
	float     maxabs_cand;
	float     maxabs_temp;
	int       i;

	bli_sabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bli_sabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bli_dmaxabsv( int n, double* x, int incx, double* maxabs )
{
	double*   chi;
	double    maxabs_cand;
	double    maxabs_temp;
	int       i;

	bli_dabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bli_dabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bli_cmaxabsv( int n, scomplex* x, int incx, float* maxabs )
{
	scomplex* chi;
	float     maxabs_cand;
	float     maxabs_temp;
	int       i;

	bli_csabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bli_csabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}

void bli_zmaxabsv( int n, dcomplex* x, int incx, double* maxabs )
{
	dcomplex* chi;
	double    maxabs_cand;
	double    maxabs_temp;
	int       i;

	bli_zdabsval2( x, &maxabs_cand );

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bli_zdabsval2( chi, &maxabs_temp );
		
		if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
	}

	*maxabs = maxabs_cand;
}


