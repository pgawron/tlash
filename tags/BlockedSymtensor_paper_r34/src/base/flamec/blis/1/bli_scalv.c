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

void bli_sscalv( conj_t conj, int n, float* alpha, float* x, int incx )
{
	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_seq1( alpha ) ) return;

	bli_sscal( n,
	           alpha,
	           x, incx );
}

void bli_dscalv( conj_t conj, int n, double* alpha, double* x, int incx )
{
	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_deq1( alpha ) ) return;

	bli_dscal( n,
	           alpha,
	           x, incx );
}

void bli_csscalv( conj_t conj, int n, float* alpha, scomplex* x, int incx )
{
	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_seq1( alpha ) ) return;

	bli_csscal( n,
	            alpha,
	            x, incx );
}

void bli_cscalv( conj_t conj, int n, scomplex* alpha, scomplex* x, int incx )
{
	scomplex alpha_conj;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_ceq1( alpha ) ) return;

	bli_ccopys( conj, alpha, &alpha_conj );

	bli_cscal( n,
	           &alpha_conj,
	           x, incx );
}

void bli_zdscalv( conj_t conj, int n, double* alpha, dcomplex* x, int incx )
{
	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_deq1( alpha ) ) return;

	bli_zdscal( n,
	            alpha,
	            x, incx );
}

void bli_zscalv( conj_t conj, int n, dcomplex* alpha, dcomplex* x, int incx )
{
	dcomplex alpha_conj;

	// Return early if possible.
	if ( bli_zero_dim1( n ) ) return;
	if ( bli_zeq1( alpha ) ) return;

	bli_zcopys( conj, alpha, &alpha_conj );

	bli_zscal( n,
	           &alpha_conj,
	           x, incx );
}

