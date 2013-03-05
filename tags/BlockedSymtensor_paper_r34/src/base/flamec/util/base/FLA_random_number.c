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

#include "FLAME.h"

float FLA_random_float()
{
  return ( float )( ( ( double ) rand() / ( double ) RAND_MAX ) * 2.0 - 1.0 );
}

double FLA_random_double()
{
  return ( ( double ) rand() / ( double ) RAND_MAX ) * 2.0 - 1.0;
}

scomplex FLA_random_scomplex()
{
  scomplex z;

  z.real = FLA_random_float();
  z.imag = FLA_random_float();

  return z;
}

dcomplex FLA_random_dcomplex()
{
  dcomplex z;

  z.real = FLA_random_double();
  z.imag = FLA_random_double();

  return z;
}

