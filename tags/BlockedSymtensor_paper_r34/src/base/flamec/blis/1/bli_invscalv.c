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

void bli_sinvscalv( conj_t conj, int n, float* alpha, float* x, int incx )
{
	float alpha_inv;

	if ( bli_seq1( alpha ) ) return;

	alpha_inv = 1.0F / *alpha;

	bli_sscal( n,
	           &alpha_inv,
	           x, incx );
}

void bli_dinvscalv( conj_t conj, int n, double* alpha, double* x, int incx )
{
	double alpha_inv;

	if ( bli_deq1( alpha ) ) return;

	alpha_inv = 1.0 / *alpha;

	bli_dscal( n,
	           &alpha_inv,
	           x, incx );
}

void bli_csinvscalv( conj_t conj, int n, float* alpha, scomplex* x, int incx )
{
	float alpha_inv;

	if ( bli_seq1( alpha ) ) return;

	alpha_inv = 1.0F / *alpha;

	bli_csscal( n,
	            &alpha_inv,
	            x, incx );
}

void bli_cinvscalv( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx )
{
	scomplex alpha_inv;

	if ( bli_ceq1( alpha ) ) return;

	bli_cinvert2s( conj, alpha, &alpha_inv );

	bli_cscal( n,
	           &alpha_inv,
	           x, incx );
}

void bli_zdinvscalv( conj_t conj, int n, double* alpha, dcomplex* x, int incx )
{
	double alpha_inv;

	if ( bli_deq1( alpha ) ) return;

	alpha_inv = 1.0 / *alpha;

	bli_zdscal( n,
	            &alpha_inv,
	            x, incx );
}

void bli_zinvscalv( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx )
{
	dcomplex alpha_inv;

	if ( bli_zeq1( alpha ) ) return;

	bli_zinvert2s( conj, alpha, &alpha_inv );

	bli_zscal( n,
	           &alpha_inv,
	           x, incx );
}

