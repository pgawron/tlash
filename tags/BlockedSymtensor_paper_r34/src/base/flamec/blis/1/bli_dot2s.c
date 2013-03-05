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

void bli_sdot2s( conj_t conj, int n, float* alpha, float* x, int incx, float* y, int incy, float* beta, float* rho )
{
	float dot;

	bli_sdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0F * (*alpha) * dot;
}

void bli_ddot2s( conj_t conj, int n, double* alpha, double* x, int incx, double* y, int incy, double* beta, double* rho )
{
	double dot;

	bli_ddot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0 * (*alpha) * dot;
}

void bli_cdot2s( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho )
{
	scomplex dotxy;
	scomplex dotyx;
	scomplex alpha_d    = *alpha;
	scomplex alphac_d   = *alpha;
	scomplex beta_d     = *beta;
	scomplex rho_d      = *rho;

	alphac_d.imag *= -1.0F;

	bli_cdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dotxy );

	bli_cdot( conj,
	          n,
	          y, incy,
	          x, incx,
	          &dotyx );

	rho->real = beta_d.real   * rho_d.real - beta_d.imag   * rho_d.imag +
	            alpha_d.real  * dotxy.real - alpha_d.imag  * dotxy.imag +
	            alphac_d.real * dotyx.real - alphac_d.imag * dotyx.imag; 
	rho->imag = beta_d.real   * rho_d.imag + beta_d.imag   * rho_d.real +
	            alpha_d.real  * dotxy.imag + alpha_d.imag  * dotxy.real +
	            alphac_d.real * dotyx.imag + alphac_d.imag * dotyx.real; 
}

void bli_zdot2s( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho )
{
	dcomplex dotxy;
	dcomplex dotyx;
	dcomplex alpha_d    = *alpha;
	dcomplex alphac_d   = *alpha;
	dcomplex beta_d     = *beta;
	dcomplex rho_d      = *rho;

	alphac_d.imag *= -1.0;

	bli_zdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dotxy );

	bli_zdot( conj,
	          n,
	          y, incy,
	          x, incx,
	          &dotyx );

	rho->real = beta_d.real   * rho_d.real - beta_d.imag   * rho_d.imag +
	            alpha_d.real  * dotxy.real - alpha_d.imag  * dotxy.imag +
	            alphac_d.real * dotyx.real - alphac_d.imag * dotyx.imag; 
	rho->imag = beta_d.real   * rho_d.imag + beta_d.imag   * rho_d.real +
	            alpha_d.real  * dotxy.imag + alpha_d.imag  * dotxy.real +
	            alphac_d.real * dotyx.imag + alphac_d.imag * dotyx.real; 
}

