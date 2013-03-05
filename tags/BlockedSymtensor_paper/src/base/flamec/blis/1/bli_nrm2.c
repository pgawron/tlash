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

void bli_snrm2( int n, float* x, int incx, float* norm )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*norm = cblas_snrm2( n,
	                     x, incx );
#else
	*norm = F77_snrm2( &n,
	                   x, &incx );
#endif
}

void bli_dnrm2( int n, double* x, int incx, double* norm )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dnrm2( n,
	                     x, incx );
#else
	*norm = F77_dnrm2( &n,
	                   x, &incx );
#endif
}

void bli_cnrm2( int n, scomplex* x, int incx, float* norm )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*norm = cblas_scnrm2( n,
	                      x, incx );
#else
	*norm = F77_scnrm2( &n,
	                    x, &incx );
#endif
}

void bli_znrm2( int n, dcomplex* x, int incx, double* norm )
{
#ifdef BLIS_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dznrm2( n,
	                      x, incx );
#else
	*norm = F77_dznrm2( &n,
	                    x, &incx );
#endif
}

