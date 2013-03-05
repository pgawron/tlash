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

void bli_sewinvscalv( conj_t conj, int n, float* x, int incx, float* y, int incy )
{
	float*    chi;
	float*    psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bli_sinvscals( chi, psi );
	}
}

void bli_dewinvscalv( conj_t conj, int n, double* x, int incx, double* y, int incy )
{
	double*   chi;
	double*   psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bli_dinvscals( chi, psi );
	}
}

void bli_csewinvscalv( conj_t conj, int n, float* x, int incx, scomplex* y, int incy )
{
	float*    chi;
	scomplex* psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bli_csinvscals( chi, psi );
	}
}

void bli_cewinvscalv( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy )
{
	scomplex* chi;
	scomplex* psi;
	scomplex  conjchi;
	int       i;

	if ( bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bli_ccopyconj( chi, &conjchi );
			bli_cinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bli_cinvscals( chi, psi );
		}
	}
}

void bli_zdewinvscalv( conj_t conj, int n, double* x, int incx, dcomplex* y, int incy )
{
	double*   chi;
	dcomplex* psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bli_zdinvscals( chi, psi );
	}
}

void bli_zewinvscalv( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy )
{
	dcomplex* chi;
	dcomplex* psi;
	dcomplex  conjchi;
	int       i;

	if ( bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bli_zcopyconj( chi, &conjchi );
			bli_zinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bli_zinvscals( chi, psi );
		}
	}
}

