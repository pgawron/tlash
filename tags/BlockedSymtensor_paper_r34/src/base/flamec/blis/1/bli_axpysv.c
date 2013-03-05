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

void bli_saxpysv( int n, float* alpha0, float* alpha1, float* x, int incx, float* beta, float* y, int incy )
{
	float    alpha_prod;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bli_sscal( n,
	           beta,
	           y, incy );

	bli_saxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bli_daxpysv( int n, double* alpha0, double* alpha1, double* x, int incx, double* beta, double* y, int incy )
{
	double   alpha_prod;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bli_dscal( n,
	           beta,
	           y, incy );

	bli_daxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bli_caxpysv( int n, scomplex* alpha0, scomplex* alpha1, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex alpha_prod;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bli_cscal( n,
	           beta,
	           y, incy );

	bli_caxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bli_zaxpysv( int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex alpha_prod;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bli_zscal( n,
	           beta,
	           y, incy );

	bli_zaxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

