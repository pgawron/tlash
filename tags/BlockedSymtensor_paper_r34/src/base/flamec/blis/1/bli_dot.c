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

void bli_sdot( conj_t conj, int n, float* x, int incx, float* y, int incy, float* rho )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*rho = cblas_sdot( n,
	                   x, incx,
	                   y, incy );
#else
	*rho = F77_sdot( &n,
	                 x, &incx,
	                        y, &incy );
#endif
}

void bli_ddot( conj_t conj, int n, double* x, int incx, double* y, int incy, double* rho )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*rho = cblas_ddot( n,
	                   x, incx,
	                   y, incy );
#else
	*rho = F77_ddot( &n,
	                 x, &incx,
	                        y, &incy );
#endif
}

void bli_cdot( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	if ( bli_is_conj( conj ) )
	{
	    cblas_cdotc_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
	else // if ( !bli_is_conj( conj ) )
	{
	    cblas_cdotu_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
#else
	bli_cdot_in( conj,
	             n,
	             x, incx,
	             y, incy,
	             rho );
#endif
}

void bli_zdot( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	if ( bli_is_conj( conj ) )
	{
	    cblas_zdotc_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
	else // if ( !bli_is_conj( conj ) )
	{
	    cblas_zdotu_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
#else
	bli_zdot_in( conj,
	             n,
	             x, incx,
	             y, incy,
	             rho );
#endif
}


// --- Inlined helper implementations ---

void bli_cdot_in( conj_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho )
{
	scomplex* xip;
	scomplex* yip;
	scomplex  xi;
	scomplex  yi;
	scomplex  rho_temp;
	int       i;

	rho_temp.real = 0.0F;
	rho_temp.imag = 0.0F;
		
	xip = x;
	yip = y;
		
	if ( bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - -xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + -xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	else // if ( !bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	
	rho->real = rho_temp.real;
	rho->imag = rho_temp.imag;
}

void bli_zdot_in( conj_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho )
{
	dcomplex* xip;
	dcomplex* yip;
	dcomplex  xi;
	dcomplex  yi;
	dcomplex  rho_temp;
	int       i;

	rho_temp.real = 0.0;
	rho_temp.imag = 0.0;
		
	xip = x;
	yip = y;
		
	if ( bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - -xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + -xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	else // if ( !bli_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	
	rho->real = rho_temp.real;
	rho->imag = rho_temp.imag;
}

