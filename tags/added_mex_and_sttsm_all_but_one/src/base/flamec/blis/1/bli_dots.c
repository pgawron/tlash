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

void bli_sdots( conj_t conj, int n, float* alpha, float* x, int incx, float* y, int incy, float* beta, float* rho )
{
	float dot_prod;

	bli_sdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	*rho = (*beta) * (*rho) + (*alpha) * dot_prod;
}

void bli_ddots( conj_t conj, int n, double* alpha, double* x, int incx, double* y, int incy, double* beta, double* rho )
{
	double dot_prod;

	bli_ddot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	*rho = (*beta) * (*rho) + (*alpha) * dot_prod;
}

void bli_cdots( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho )
{
	scomplex rho_orig = *rho;
	scomplex dot_prod;

	bli_cdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	rho->real = beta->real  * rho_orig.real - beta->imag  * rho_orig.imag +
	            alpha->real * dot_prod.real - alpha->imag * dot_prod.imag;
	rho->imag = beta->real  * rho_orig.imag + beta->imag  * rho_orig.real +
	            alpha->real * dot_prod.imag + alpha->imag * dot_prod.real;
}

void bli_zdots( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho )
{
	dcomplex rho_orig = *rho;
	dcomplex dot_prod;

	bli_zdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	rho->real = beta->real  * rho_orig.real - beta->imag  * rho_orig.imag +
	            alpha->real * dot_prod.real - alpha->imag * dot_prod.imag;
	rho->imag = beta->real  * rho_orig.imag + beta->imag  * rho_orig.real +
	            alpha->real * dot_prod.imag + alpha->imag * dot_prod.real;
}

