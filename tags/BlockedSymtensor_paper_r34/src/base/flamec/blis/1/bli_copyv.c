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

void bli_icopyv( conj_t conj, int m, int* x, int incx, int* y, int incy )
{
	int*      chi;
	int*      psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

void bli_scopyv( conj_t conj, int m, float* x, int incx, float* y, int incy )
{
	bli_scopy( m,
	           x, incx, 
	           y, incy );
}

void bli_dcopyv( conj_t conj, int m, double* x, int incx, double* y, int incy )
{
	bli_dcopy( m,
	           x, incx, 
	           y, incy );
}

void bli_ccopyv( conj_t conj, int m, scomplex* x, int incx, scomplex* y, int incy )
{
	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	bli_ccopy( m,
	           x, incx, 
	           y, incy );

	if ( bli_is_conj( conj ) )
		bli_cconjv( m,
	                y, incy );
}

void bli_zcopyv( conj_t conj, int m, dcomplex* x, int incx, dcomplex* y, int incy )
{
	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	bli_zcopy( m,
	           x, incx, 
	           y, incy );

	if ( bli_is_conj( conj ) )
		bli_zconjv( m,
	                y, incy );
}

// --- Mixed-datatype and general stride copy routines---------------

// sd ds
void bli_sdcopyv( conj_t conj, int m, float* x, int incx, double* y, int incy )
{
	float*    chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}
void bli_dscopyv( conj_t conj, int m, double* x, int incx, float* y, int incy )
{
	double*   chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

// sc cs
void bli_sccopyv( conj_t conj, int m, float* x, int incx, scomplex* y, int incy )
{
	float*    chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bli_cscopyv( conj_t conj, int m, scomplex* x, int incx, float* y, int incy )
{
	scomplex* chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// sz zs
void bli_szcopyv( conj_t conj, int m, float* x, int incx, dcomplex* y, int incy )
{
	float*    chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bli_zscopyv( conj_t conj, int m, dcomplex* x, int incx, float* y, int incy )
{
	dcomplex* chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dc cd
void bli_dccopyv( conj_t conj, int m, double* x, int incx, scomplex* y, int incy )
{
	double*   chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bli_cdcopyv( conj_t conj, int m, scomplex* x, int incx, double* y, int incy )
{
	scomplex* chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dz zd
void bli_dzcopyv( conj_t conj, int m, double* x, int incx, dcomplex* y, int incy )
{
	double*   chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bli_zdcopyv( conj_t conj, int m, dcomplex* x, int incx, double* y, int incy )
{
	dcomplex* chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// cz zc
void bli_czcopyv( conj_t conj, int m, scomplex* x, int incx, dcomplex* y, int incy )
{
	scomplex* chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bli_is_conj( conj ) )
		bli_zconjv( m,
	                y, incy );
}
void bli_zccopyv( conj_t conj, int m, dcomplex* x, int incx, scomplex* y, int incy )
{
	dcomplex* chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bli_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bli_is_conj( conj ) )
		bli_cconjv( m,
	                y, incy );
}

