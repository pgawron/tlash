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

#ifdef BLIS_ENABLE_USE_OF_FLA_MALLOC
  #include "FLAME.h"
  #define BLIS_MALLOC FLA_malloc
#else
  #define BLIS_MALLOC malloc
#endif

void*     bli_vallocv( unsigned int n_elem, unsigned int elem_size )
{
	return ( void*  ) BLIS_MALLOC( n_elem * elem_size );
}

int*      bli_iallocv( unsigned int n_elem )
{
	return ( int*   ) BLIS_MALLOC( n_elem * sizeof( int ) );
}

float*    bli_sallocv( unsigned int n_elem )
{
	return ( float* ) BLIS_MALLOC( n_elem * sizeof( float ) );
}

double*   bli_dallocv( unsigned int n_elem )
{
	return ( double* ) BLIS_MALLOC( n_elem * sizeof( double ) );
}

scomplex* bli_callocv( unsigned int n_elem )
{
	return ( scomplex* ) BLIS_MALLOC( n_elem * sizeof( scomplex ) );
}

dcomplex* bli_zallocv( unsigned int n_elem )
{
	return ( dcomplex* ) BLIS_MALLOC( n_elem * sizeof( dcomplex ) );
}

