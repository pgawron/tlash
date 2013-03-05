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

FLA_Error FLA_Sv_2x2( FLA_Obj alpha11, FLA_Obj alpha12, FLA_Obj alpha22,
                      FLA_Obj sigma1, FLA_Obj sigma2 )
/*
  Compute the singular value decomposition of a 2x2 triangular matrix A
  such that

    / alpha11 alpha12 \
    \    0    alpha22 /

  is equal to

    / gammaL -sigmaL \ / sigma1    0    \ / gammaR -sigmaR \'
    \ sigmaL  gammaL / \    0    sigma2 / \ sigmaR  gammaR /

  Upon completion, sigma1 and sigma2 are overwritten with the
  singular values of smaller and larger absolute values, respectively.
  gammaL, sigmaL, gammaR, and sigmaR are not computed.

  This routine is a nearly-verbatim translation of slas2() and dlas2()
  from the netlib distribution of LAPACK.

  -FGVZ
*/
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( alpha11 );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  buff_alpha11 = FLA_FLOAT_PTR( alpha11 );
			float*  buff_alpha12 = FLA_FLOAT_PTR( alpha12 );
			float*  buff_alpha22 = FLA_FLOAT_PTR( alpha22 );
			float*  buff_sigma1  = FLA_FLOAT_PTR( sigma1 );
			float*  buff_sigma2  = FLA_FLOAT_PTR( sigma2 );

			FLA_Sv_2x2_ops( buff_alpha11,
			                buff_alpha12,
			                buff_alpha22,
			                buff_sigma1,
			                buff_sigma2 );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_alpha11 = FLA_DOUBLE_PTR( alpha11 );
			double* buff_alpha12 = FLA_DOUBLE_PTR( alpha12 );
			double* buff_alpha22 = FLA_DOUBLE_PTR( alpha22 );
			double* buff_sigma1  = FLA_DOUBLE_PTR( sigma1 );
			double* buff_sigma2  = FLA_DOUBLE_PTR( sigma2 );

			FLA_Sv_2x2_opd( buff_alpha11,
			                buff_alpha12,
			                buff_alpha22,
			                buff_sigma1,
			                buff_sigma2 );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Sv_2x2_ops( float*    alpha11,
                          float*    alpha12,
                          float*    alpha22,
                          float*    sigma1,
                          float*    sigma2 )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Sv_2x2_opd( double*   alpha11,
                          double*   alpha12,
                          double*   alpha22,
                          double*   sigma1,
                          double*   sigma2 )
{
	double zero = 0.0;
	double one  = 1.0;
	double two  = 2.0;

	double f, g, h;
	double as, at, au, c, fa, fhmin, fhmax, ga, ha;
	double ssmin, ssmax;
	double temp, temp2;

	f = *alpha11;
	g = *alpha12;
	h = *alpha22;

	fa = fabs( f );
	ga = fabs( g );
	ha = fabs( h );

	fhmin = min( fa, ha );
	fhmax = max( fa, ha );

	if ( fhmin == zero )
	{
		ssmin = zero;

		if ( fhmax == zero )
			ssmax = ga;
		else
		{
			temp = min( fhmax, ga ) / max( fhmax, ga );
			ssmax = max( fhmax, ga ) * sqrt( one + temp * temp );
		}
	}
	else
	{
		if ( ga < fhmax )
		{
			as = one + fhmin / fhmax;
			at = ( fhmax - fhmin ) / fhmax;
			au = ( ga / fhmax ) * ( ga / fhmax );
			c  = two / ( sqrt( as * as + au ) + sqrt( at * at + au ) );
			ssmin = fhmin * c;
			ssmax = fhmax / c;
		}
		else
		{
			au = fhmax / ga;

			if ( au == zero )
			{
				ssmin = ( fhmin * fhmax ) / ga;
				ssmax = ga;
			}
			else
			{
				as = one + fhmin / fhmax;
				at = ( fhmax - fhmin ) / fhmax;
				temp  = as * au;
				temp2 = at * au;
				c  = one / ( sqrt( one + temp * temp ) +
				             sqrt( one + temp2 * temp2 ) );
				ssmin = ( fhmin * c ) * au;
				ssmin = ssmin + ssmin;
				ssmax = ga / ( c + c );
			}
		}
	}

	// Save the output values.

	*sigma1 = ssmin;
	*sigma2 = ssmax;

	return FLA_SUCCESS;
}

