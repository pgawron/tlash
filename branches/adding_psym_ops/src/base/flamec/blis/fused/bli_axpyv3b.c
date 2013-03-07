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

/*
   Effective computation:

   y = y + alpha1 * x1;
   y = y + alpha2 * x2;
   y = y + alpha3 * x3;
*/

void bli_saxpyv3b( int       n,
                   float*    alpha1,
                   float*    alpha2,
                   float*    alpha3,
                   float*    x1, int inc_x1,
                   float*    x2, int inc_x2,
                   float*    x3, int inc_x3,
                   float*    y,  int inc_y )
{
	bli_abort();
}


void bli_daxpyv3b( int       n,
                   double*   alpha1,
                   double*   alpha2,
                   double*   alpha3,
                   double*   x1, int inc_x1,
                   double*   x2, int inc_x2,
                   double*   x3, int inc_x3,
                   double*   y,  int inc_y )
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict chi2;
	double*   restrict chi3;
	double*   restrict psi1;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	v2df_t    a1v, a2v, a3v;
	v2df_t    x11v, x12v, x13v;
	v2df_t    x21v, x22v, x23v;
	v2df_t    y1v;
	v2df_t    y2v;

	if ( inc_x1 != 1 ||
	     inc_x2 != 1 ||
	     inc_x3 != 1 ||
	     inc_y  != 1 ) bli_abort();

	n_pre = 0;
	if ( ( unsigned long ) y % 16 != 0 )
	{
		if ( ( unsigned long ) x1 % 16 == 0 ||
		     ( unsigned long ) x2 % 16 == 0 ||
		     ( unsigned long ) x3 % 16 == 0 ) bli_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	chi1 = x1;
	chi2 = x2;
	chi3 = x3;
	psi1 = y;

	if ( n_pre == 1 )
	{
		double   alpha1_c = *alpha1;
		double   alpha2_c = *alpha2;
		double   alpha3_c = *alpha3;
		double   chi11_c = *chi1;
		double   chi12_c = *chi2;
		double   chi13_c = *chi3;

		*psi1 += alpha1_c * chi11_c + alpha2_c * chi12_c + alpha3_c * chi13_c;

		chi1 += inc_x1;
		chi2 += inc_x2;
		chi3 += inc_x3;
		psi1 += inc_y;
	}

	a1v.v = _mm_loaddup_pd( ( double* )alpha1 );
	a2v.v = _mm_loaddup_pd( ( double* )alpha2 );
	a3v.v = _mm_loaddup_pd( ( double* )alpha3 );

	for ( i = 0; i < n_run; ++i )
	{
		x11v.v = _mm_load_pd( ( double* )chi1 );
		x12v.v = _mm_load_pd( ( double* )chi2 );
		x13v.v = _mm_load_pd( ( double* )chi3 );
		y1v.v  = _mm_load_pd( ( double* )psi1 );

		y1v.v += a1v.v * x11v.v + a2v.v * x12v.v + a3v.v * x13v.v;

		_mm_store_pd( ( double* )psi1, y1v.v );

		x21v.v = _mm_load_pd( ( double* )(chi1 + 2) );
		x22v.v = _mm_load_pd( ( double* )(chi2 + 2) );
		x23v.v = _mm_load_pd( ( double* )(chi3 + 2) );
		y2v.v  = _mm_load_pd( ( double* )(psi1 + 2) );

		y2v.v += a1v.v * x21v.v + a2v.v * x22v.v + a3v.v * x23v.v;

		_mm_store_pd( ( double* )(psi1 + 2), y2v.v );

		chi1 += 4;
		chi2 += 4;
		chi3 += 4;
		psi1 += 4;
	}

	if ( n_left > 0 )
	{
		double   alpha1_c = *alpha1;
		double   alpha2_c = *alpha2;
		double   alpha3_c = *alpha3;

		for ( i = 0; i < n_left; ++i )
		{
			double   chi11_c = *chi1;
			double   chi12_c = *chi2;
			double   chi13_c = *chi3;

			*psi1 += alpha1_c * chi11_c + alpha2_c * chi12_c + alpha3_c * chi13_c;

			chi1 += inc_x1;
			chi2 += inc_x2;
			chi3 += inc_x3;
			psi1 += inc_y;
		}
	}
}
#elif BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict chi2;
	double*   restrict psi1;
	double    alpha1_c;
	double    alpha2_c;
	double    temp1;
	double    temp2;
	int       i;

	int       n_run  = n / 2;
	int       n_left = n % 2;
	int       twoinc_x1 = 2*inc_x1;
	int       twoinc_x2 = 2*inc_x2;
	int       twoinc_y  = 2*inc_y;

	chi1 = x1;
	chi2 = x2;
	psi1 = y;

	alpha1_c = *alpha1;
	alpha2_c = *alpha2;

	for ( i = 0; i < n_run; ++i )
	{
		double   chi11_c = *chi1;
		double   chi21_c = *(chi1 + inc_x1);
		double   chi12_c = *chi2;
		double   chi22_c = *(chi2 + inc_x2);
		double   psi1_c  = *psi1;
		double   psi2_c  = *(psi1 + inc_y);

		// psi1 = psi1 + alpha1 * chi11 + alpha2 * chi12;
		// psi2 = psi2 + alpha1 * chi21 + alpha2 * chi22;
		temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;
		temp2 = alpha1_c * chi21_c + alpha2_c * chi22_c;

		*psi1           = psi1_c + temp1;
		*(psi1 + inc_y) = psi2_c + temp2;

		chi1 += twoinc_x1;
		chi2 += twoinc_x2;
		psi1 += twoinc_y;
	}

	if ( n_left == 1 )
	{
		double   chi11_c = *chi1;
		double   chi12_c = *chi2;

		// psi1 = psi1 + alpha1 * chi11 + alpha2 * chi12;
		temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;

		*psi1 = *psi1 + temp1;
	}
}
#endif


void bli_caxpyv3b( int       n,
                   scomplex* alpha1,
                   scomplex* alpha2,
                   scomplex* alpha3,
                   scomplex* x1, int inc_x1,
                   scomplex* x2, int inc_x2,
                   scomplex* x3, int inc_x3,
                   scomplex* y,  int inc_y )
{
	bli_abort();
}


void bli_zaxpyv3b( int       n,
                   dcomplex* alpha1,
                   dcomplex* alpha2,
                   dcomplex* alpha3,
                   dcomplex* x1, int inc_x1,
                   dcomplex* x2, int inc_x2,
                   dcomplex* x3, int inc_x3,
                   dcomplex* y,  int inc_y )
{
	bli_abort();
}

