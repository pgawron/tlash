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

void bli_sident( int m, float* a, int a_rs, int a_cs )
{
	float* alpha;
	int    i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0F;

			if ( i == j )
				*alpha = 1.0F;
		}
	}
}

void bli_dident( int m, double* a, int a_rs, int a_cs )
{
	double* alpha;
	int     i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0;

			if ( i == j )
				*alpha = 1.0;
		}
	}
}

void bli_cident( int m, scomplex* a, int a_rs, int a_cs )
{
	scomplex* alpha;
	int       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0F;
			alpha->imag = 0.0F;

			if ( i == j )
				alpha->real = 1.0F;
		}
	}
}

void bli_zident( int m, dcomplex* a, int a_rs, int a_cs )
{
	dcomplex* alpha;
	int       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0;
			alpha->imag = 0.0;

			if ( i == j )
				alpha->real = 1.0;
		}
	}
}

