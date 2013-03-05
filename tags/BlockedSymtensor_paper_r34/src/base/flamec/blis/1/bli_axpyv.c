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

void bli_saxpyv( conj_t conj, int n, float* alpha, float* x, int incx, float* y, int incy )
{
	bli_saxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bli_daxpyv( conj_t conj, int n, double* alpha, double* x, int incx, double* y, int incy )
{
	bli_daxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bli_caxpyv( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy )
{
	scomplex* x_copy;
	int       incx_copy;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	x_copy    = x;
	incx_copy = incx;
	
	if ( bli_is_conj( conj ) )
	{
		x_copy    = bli_callocv( n );
		incx_copy = 1;
	
		bli_ccopyv( conj,
		            n,
		            x,      incx,
		            x_copy, incx_copy );
	}

	bli_caxpy( n,
	           alpha,
	           x_copy, incx_copy,
	           y,      incy );

	if ( bli_is_conj( conj ) )
		bli_cfree( x_copy );
}

void bli_zaxpyv( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy )
{
	dcomplex* x_copy;
	int       incx_copy;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;

	x_copy    = x;
	incx_copy = incx;
	
	if ( bli_is_conj( conj ) )
	{
		x_copy    = bli_zallocv( n );
		incx_copy = 1;
	
		bli_zcopyv( conj,
		            n,
		            x,      incx,
		            x_copy, incx_copy );
	}

	bli_zaxpy( n,
	           alpha,
	           x_copy, incx_copy,
	           y,      incy );

	if ( bli_is_conj( conj ) )
		bli_zfree( x_copy );
}

