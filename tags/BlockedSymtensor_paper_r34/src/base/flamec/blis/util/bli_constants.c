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

// --- two ---

float bli_s2( void )
{
	float x;
	x = 2.0F;
	return x;
}

double bli_d2( void )
{
	double x;
	x = 2.0;
	return x;
}

scomplex bli_c2( void )
{
	scomplex x;
	x.real = bli_s2();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_z2( void )
{
	dcomplex x;
	x.real = bli_d2();
	x.imag = bli_d0();
	return x;
}

// --- one ---

float bli_s1( void )
{
	float x;
	x = 1.0F;
	return x;
}

double bli_d1( void )
{
	double x;
	x = 1.0;
	return x;
}

scomplex bli_c1( void )
{
	scomplex x;
	x.real = bli_s1();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_z1( void )
{
	dcomplex x;
	x.real = bli_d1();
	x.imag = bli_d0();
	return x;
}

// --- one half ---

float bli_s1h( void )
{
	float x;
	x = 0.5F;
	return x;
}

double bli_d1h( void )
{
	double x;
	x = 0.5;
	return x;
}

scomplex bli_c1h( void )
{
	scomplex x;
	x.real = bli_s1h();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_z1h( void )
{
	dcomplex x;
	x.real = bli_d1h();
	x.imag = bli_d0();
	return x;
}

// --- zero ---

float bli_s0( void )
{
	float x;
	x = 0.0F;
	return x;
}

double bli_d0( void )
{
	double x;
	x = 0.0;
	return x;
}

scomplex bli_c0( void )
{
	scomplex x;
	x.real = bli_s0();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_z0( void )
{
	dcomplex x;
	x.real = bli_d0();
	x.imag = bli_d0();
	return x;
}

// --- minus one half ---

float bli_sm1h( void )
{
	float x;
	x = -0.5F;
	return x;
}

double bli_dm1h( void )
{
	double x;
	x = -0.5;
	return x;
}

scomplex bli_cm1h( void )
{
	scomplex x;
	x.real = bli_sm1h();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_zm1h( void )
{
	dcomplex x;
	x.real = bli_dm1h();
	x.imag = bli_d0();
	return x;
}

// --- minus one ---

float bli_sm1( void )
{
	float x;
	x = -1.0F;
	return x;
}

double bli_dm1( void )
{
	double x;
	x = -1.0;
	return x;
}

scomplex bli_cm1( void )
{
	scomplex x;
	x.real = bli_sm1();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_zm1( void )
{
	dcomplex x;
	x.real = bli_dm1();
	x.imag = bli_d0();
	return x;
}

// --- minus two ---

float bli_sm2( void )
{
	float x;
	x = -2.0F;
	return x;
}

double bli_dm2( void )
{
	double x;
	x = -2.0;
	return x;
}

scomplex bli_cm2( void )
{
	scomplex x;
	x.real = bli_sm2();
	x.imag = bli_s0();
	return x;
}

dcomplex bli_zm2( void )
{
	dcomplex x;
	x.real = bli_dm2();
	x.imag = bli_d0();
	return x;
}

