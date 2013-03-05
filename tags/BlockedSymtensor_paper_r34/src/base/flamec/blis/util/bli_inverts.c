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

void bli_sinverts( conj_t conj, float* alpha )
{
	float  one = 1.0F;

	*alpha = one / *alpha;
}

void bli_dinverts( conj_t conj, double* alpha )
{
	double one = 1.0;

	*alpha = one / *alpha;
}

void bli_cinverts( conj_t conj, scomplex* alpha )
{
	float  one = 1.0F;
	float  temp;

    temp = one / ( alpha->real * alpha->real +
                   alpha->imag * alpha->imag );
    alpha->real = alpha->real *  temp;
    alpha->imag = alpha->imag * -temp;

	if ( bli_is_conj( conj ) )
		bli_cconjs( alpha );
}

void bli_zinverts( conj_t conj, dcomplex* alpha )
{
	double one = 1.0;
	double temp;

    temp = one / ( alpha->real * alpha->real +
                   alpha->imag * alpha->imag );
    alpha->real = alpha->real *  temp;
    alpha->imag = alpha->imag * -temp;

	if ( bli_is_conj( conj ) )
		bli_zconjs( alpha );
}

