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

/*
   Effective computation:

     rho_xz = beta * rho_xz + x * z;
     rho_yz = beta * rho_yz + y * z;

   where x and y are optionally conjugated.
*/

dcomplex zdotc_( int*      n,
                 dcomplex* x, int* inc_x,
                 dcomplex* z, int* inc_z )
{
	dcomplex* restrict x1;
	dcomplex* restrict z1;
	int                i;
	v2df_t rho1v;
	v2df_t z11v, z12v;
	v2df_t x1v, x1rv;
	dcomplex rho;
	int    n1 = *n;
	int    incx = *inc_x;
	int    incz = *inc_z;

	x1 = x;
	z1 = z;

	rho1v.v = _mm_setzero_pd();

	{
		v2df_t bcac, adbd;

		for ( i = 0; i < n1; ++i )
		{
			z11v.v = _mm_loaddup_pd( ( double* )&(z1->real) );
			z12v.v = _mm_loaddup_pd( ( double* )&(z1->imag) );

			x1v.v  = _mm_load_pd( ( double* )x1 );
			x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
			bcac.v = x1rv.v * z11v.v;
			adbd.v = x1v.v  * z12v.v;
			rho1v.v = rho1v.v + _mm_addsub_pd( bcac.v, adbd.v );

			x1 += incx;
			z1 += incz;
		}

		rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );

		rho1v.d[1] = -rho1v.d[1];
	}

	rho.real = rho1v.d[0];
	rho.imag = rho1v.d[1];

	return rho;
}

