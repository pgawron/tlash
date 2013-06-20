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

void bli_sinvert2s( conj_t conj, float* alpha, float* beta )
{
	float  one = 1.0F;

	*beta = one / *alpha;
}

void bli_dinvert2s( conj_t conj, double* alpha, double* beta )
{
	double one = 1.0;

	*beta = one / *alpha;
}

void bli_cinvert2s( conj_t conj, scomplex* alpha, scomplex* beta )
{
	float  one = 1.0F;
	float  temp;

    temp = one / ( alpha->real * alpha->real +
                   alpha->imag * alpha->imag );
    beta->real = alpha->real *  temp;
    beta->imag = alpha->imag * -temp;

	if ( bli_is_conj( conj ) )
		bli_cconjs( beta );
}

void bli_zinvert2s( conj_t conj, dcomplex* alpha, dcomplex* beta )
{
	double one = 1.0;
	double temp;

    temp = one / ( alpha->real * alpha->real +
                   alpha->imag * alpha->imag );
    beta->real = alpha->real *  temp;
    beta->imag = alpha->imag * -temp;

	if ( bli_is_conj( conj ) )
		bli_zconjs( beta );
}

