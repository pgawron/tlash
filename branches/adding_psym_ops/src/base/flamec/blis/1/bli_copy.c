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

void bli_scopy( int m, float* x, int incx, float* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	cblas_scopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_scopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bli_dcopy( int m, double* x, int incx, double* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	cblas_dcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_dcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bli_ccopy( int m, scomplex* x, int incx, scomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	cblas_ccopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_ccopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bli_zcopy( int m, dcomplex* x, int incx, dcomplex* y, int incy )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	cblas_zcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_zcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

