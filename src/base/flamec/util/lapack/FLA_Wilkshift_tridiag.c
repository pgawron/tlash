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

FLA_Error FLA_Wilkshift_tridiag( FLA_Obj delta1, FLA_Obj epsilon, FLA_Obj delta2, FLA_Obj kappa )
/*
  Compute a Wilkinson shift, kappa, which is the eigenvalue of the symmetric
  matrix

    / a - d   c \
    \   c     0 /

  that is closest to d.

  The algorithm upon which this code is based was derived from a similar
  algorithm for general matrices on page 99 of G. W. (Pete) Stewart's
  "Matrix Algorithms Volume II: Eigensystems".

  Stewart notes that the characteristic equation of this matrix is

    lambda^2 - (a - d)*lambda - c^2 = 0

  Let 2p = (a - d) and q = c^2:

    lambda^2 - 2p*lambda - q = 0

  Thus, the roots of the equation are

    lambda = p +/- sqrt( p^2 - q )
           = p +/- r

  Stewart proposes computing the smallest root indirectly to avoid
  cancellation via the relation

    lambda_min * lambda_max = -q

  To compute the larger root, Stewart notes that

    |lambda|^2 = |p +/- r|^2
               = |p|^2 +/- 2*p*r + |r|^2

  and that we should choose the +/- sign so that the middle term above
  is positive, and use this sign in

    p +/- r

  to find the larger root.

  -FGVZ
*/
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( delta1 );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Wilkshift_tridiag_check( delta1, epsilon, delta2, kappa );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* delta1_p  = ( float* ) FLA_FLOAT_PTR( delta1 );
			float* epsilon_p = ( float* ) FLA_FLOAT_PTR( epsilon );
			float* delta2_p  = ( float* ) FLA_FLOAT_PTR( delta2 );
			float* kappa_p   = ( float* ) FLA_FLOAT_PTR( kappa );

			FLA_Wilkshift_tridiag_ops( *delta1_p,
			                           *epsilon_p,
			                           *delta2_p,
			                           kappa_p );

			break;
		}

		case FLA_DOUBLE:
		{
			double* delta1_p  = ( double* ) FLA_DOUBLE_PTR( delta1 );
			double* epsilon_p = ( double* ) FLA_DOUBLE_PTR( epsilon );
			double* delta2_p  = ( double* ) FLA_DOUBLE_PTR( delta2 );
			double* kappa_p   = ( double* ) FLA_DOUBLE_PTR( kappa );

			FLA_Wilkshift_tridiag_opd( *delta1_p,
			                           *epsilon_p,
			                           *delta2_p,
			                           kappa_p );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Wilkshift_tridiag_ops( float   delta1,
                                     float   epsilon,
                                     float   delta2,
                                     float*  kappa )
{
	float  a = delta1;
	float  c = epsilon;
	float  d = delta2;
	float  p, q, r, s, k;

	// Begin with kappa equal to d.
	k = d;

	// Compute a scaling factor to promote numerical stability.
	s = fabs( a ) + 2.0F * fabs( c ) + fabs( d );

	if ( s == 0.0F ) return FLA_SUCCESS;

	// Compute q with scaling applied.
	q = ( c / s ) * ( c / s );

	if ( q != 0.0F )
	{
		// Compute p = 0.5*( a - d ) with scaling applied.
		p = 0.5F * ( ( a / s ) - ( d / s ) );

		// Compute r = sqrt( p*p - q ).
		r = sqrt( p * p + q );

		// If p*r is negative, then we need to negate r to ensure we use the
		// larger of the two eigenvalues.
		if ( p * r < 0.0F ) r = -r;

		// Compute the Wilkinson shift with scaling removed:
		//   k = lambda_min + d
		//     = d + lambda_min
		//     = d + (-q / lambda_max)
		//     = d - q / lambda_max
		//     = d - q / (p + r)
		k = k - s * ( q / ( p + r ) );
	}

	// Save the result.
	*kappa = k;

	return FLA_SUCCESS;
}

double fla_dlapy2( double x, double y );

FLA_Error FLA_Wilkshift_tridiag_opd( double  delta1,
                                     double  epsilon,
                                     double  delta2,
                                     double* kappa )
{
	double a = delta1;
	double c = epsilon;
	double d = delta2;
	double p, q, r, s, k;

	// Begin with kappa equal to d.
	k = d;

	// Compute a scaling factor to promote numerical stability.
	s = fabs( a ) + 2.0 * fabs( c ) + fabs( d );

	if ( s == 0.0 ) return FLA_SUCCESS;

	// Compute q with scaling applied.
	q = ( c / s ) * ( c / s );

	if ( q != 0.0 )
	{

		// Compute p = 0.5*( a - d ) with scaling applied.
		p = 0.5 * ( ( a / s ) - ( d / s ) );

		// Compute r = sqrt( p*p - q ).
		r = sqrt( p * p + q );

		// If p*r is negative, then we need to negate r to ensure we use the
		// larger of the two eigenvalues.
		if ( p * r < 0.0 ) r = -r;

		// Compute the Wilkinson shift with scaling removed:
		//   k = lambda_min + d
		//     = d + lambda_min
		//     = d + (-q / lambda_max)
		//     = d - q / lambda_max
		//     = d - q / (p + r)
		k = k - s * ( q / ( p + r ) );

/*
        // LAPACK method:

		p = 0.5 * ( ( a ) - ( d ) ) / c ;
		//r = sqrt( p * p + 1.0 );
		r = fla_dlapy2( p, 1.0 );
		if ( p < 0.0 ) r = -r;
		k = k - ( c / ( p + r ) );
*/
	}

	// Save the result.
	*kappa = k;

	return FLA_SUCCESS;
}

/*
double fla_dlapy2( double x, double y )
{
	double zero = 0.0;
	double one  = 1.0;
	double xabs, yabs, w, z, r;

	xabs = fabs( x );
	yabs = fabs( y );
	w    = max( xabs, yabs );
	z    = min( xabs, yabs );

	if ( z == zero )
		r = w;
	else
		r = w * sqrt( one + pow( ( z / w ), 2.0 ) );

	return r;
}
*/
