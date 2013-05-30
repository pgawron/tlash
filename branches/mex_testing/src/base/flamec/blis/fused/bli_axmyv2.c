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

     y = y - alpha * x;
     z = z - beta  * x;

   where x is optionally conjugated.
*/

void bli_saxmyv2( conj_t    conjx,
                  int       n,
                  float*    alpha,
                  float*    beta,
                  float*    x, int inc_x,
                  float*    y, int inc_y,
                  float*    z, int inc_z )
{
	bli_abort();
}


void bli_daxmyv2( conj_t    conjx,
                  int       n,
                  double*   alpha,
                  double*   beta,
                  double*   x, int inc_x,
                  double*   y, int inc_y,
                  double*   z, int inc_z )
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict psi1;
	double*   restrict zeta1;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	v2df_t    a1v, b1v;
	v2df_t    x1v, y1v, z1v;
	v2df_t    x2v, y2v, z2v;

	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bli_abort();

	n_pre = 0;
	if ( ( unsigned long ) z % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) y % 16 == 0 ) bli_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	chi1  = x;
	psi1  = y;
	zeta1 = z;

	if ( n_pre == 1 )
	{
		double   alpha_c = *alpha;
		double   beta_c  = *beta;
		double   chi1_c  = *chi1;

		*psi1  -= alpha_c * chi1_c;
		*zeta1 -= beta_c  * chi1_c;

		chi1  += inc_x;
		psi1  += inc_y;
		zeta1 += inc_z;
	}

	a1v.v = _mm_loaddup_pd( ( double* )alpha );
	b1v.v = _mm_loaddup_pd( ( double* )beta );

	for ( i = 0; i < n_run; ++i )
	{
		x1v.v = _mm_load_pd( ( double* )chi1 );
		y1v.v = _mm_load_pd( ( double* )psi1 );
		z1v.v = _mm_load_pd( ( double* )zeta1 );

		x2v.v = _mm_load_pd( ( double* )(chi1 + 2) );
		y2v.v = _mm_load_pd( ( double* )(psi1 + 2) );
		z2v.v = _mm_load_pd( ( double* )(zeta1 + 2) );

		y1v.v = y1v.v - a1v.v * x1v.v;
		z1v.v = z1v.v - b1v.v * x1v.v;

		_mm_store_pd( ( double* )psi1,  y1v.v );
		_mm_store_pd( ( double* )zeta1, z1v.v );

		y2v.v = y2v.v - a1v.v * x2v.v;
		z2v.v = z2v.v - b1v.v * x2v.v;

		_mm_store_pd( ( double* )(psi1 + 2),  y2v.v );
		_mm_store_pd( ( double* )(zeta1 + 2), z2v.v );

		chi1  += 4;
		psi1  += 4;
		zeta1 += 4;
	}

	if ( n_left > 0 )
	{
		double   alpha_c = *alpha;
		double   beta_c  = *beta;

		for( i = 0; i < n_left; ++i )
		{
			double   chi1_c = *chi1;

			*psi1  -= alpha_c * chi1_c;
			*zeta1 -= beta_c  * chi1_c;

			chi1  += inc_x;
			psi1  += inc_y;
			zeta1 += inc_z;
		}
	}
}
#elif BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict psi1;
	double*   restrict zeta1;
	double    alpha_c;
	double    beta_c;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bli_abort();

	n_pre = 0;
	if ( ( unsigned long ) z % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) y % 16 == 0 ) bli_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	chi1  = x;
	psi1  = y;
	zeta1 = z;

	alpha_c = *alpha;
	beta_c  = *beta;

	if ( n_pre == 1 )
	{
		//double   alpha_c = *alpha;
		//double   beta_c  = *beta;
		double   chi1_c  = *chi1;

		*psi1  -= alpha_c * chi1_c;
		*zeta1 -= beta_c  * chi1_c;

		chi1  += inc_x;
		psi1  += inc_y;
		zeta1 += inc_z;
	}

	for ( i = 0; i < n_run; ++i )
	{
		double   chi1_c = *chi1;
		double   chi2_c = *(chi1 + 1);
		double   psi1_c = *psi1;
		double   psi2_c = *(psi1 + 1);
		double   zeta1_c = *zeta1;
		double   zeta2_c = *(zeta1 + 1);

		// psi1  = psi1  - alpha * chi1;
		// psi2  = psi2  - alpha * chi2;
		psi1_c  -= alpha_c * chi1_c;
		psi2_c  -= alpha_c * chi2_c;

		// zeta1 = zeta1 - beta  * chi1;
		// zeta2 = zeta2 - beta  * chi2;
		zeta1_c -= beta_c  * chi1_c;
		zeta2_c -= beta_c  * chi2_c;

		*psi1        = psi1_c;
		*(psi1 + 1)  = psi2_c;
		*zeta1       = zeta1_c;
		*(zeta1 + 1) = zeta2_c;

		chi1  += 2;
		psi1  += 2;
		zeta1 += 2;
	}

	if ( n_left > 0 )
	{
		//double   alpha_c = *alpha;
		//double   beta_c  = *beta;

		for( i = 0; i < n_left; ++i )
		{
			double   chi1_c = *chi1;

			*psi1  -= alpha_c * chi1_c;
			*zeta1 -= beta_c  * chi1_c;

			chi1  += inc_x;
			psi1  += inc_y;
			zeta1 += inc_z;
		}
	}
}
#endif


void bli_caxmyv2( conj_t    conjx,
                  int       n,
                  scomplex* alpha,
                  scomplex* beta,
                  scomplex* x, int inc_x,
                  scomplex* y, int inc_y,
                  scomplex* z, int inc_z )
{
	bli_abort();
}


void bli_zaxmyv2( conj_t    conjx,
                  int       n,
                  dcomplex* alpha,
                  dcomplex* beta,
                  dcomplex* x, int inc_x,
                  dcomplex* y, int inc_y,
                  dcomplex* z, int inc_z )
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict psi1;
	dcomplex* restrict zeta1;
	dcomplex  alpha_c;
	dcomplex  beta_c;
	int       i;
	v2df_t    alphav, alpharv;
	v2df_t    betav,  betarv;
	v2df_t    x11v, x12v, y1v, z1v;
	v2df_t    acbc, bdad;

	chi1  = x;
	psi1  = y;
	zeta1 = z;

	alphav.v  = _mm_load_pd( ( double* )alpha );
	betav.v   = _mm_load_pd( ( double* )beta );
	alpharv.v = _mm_shuffle_pd( alphav.v, alphav.v, _MM_SHUFFLE2 (0,1) );
	betarv.v  = _mm_shuffle_pd( betav.v, betav.v, _MM_SHUFFLE2 (0,1) );

	if ( bli_is_conj( conjx ) )
	{
		alpha_c = *alpha;
		beta_c  = *beta;

		for ( i = 0; i < n; ++i )
		{
			dcomplex chi1_c = *chi1;

			// psi1  = psi1  + alpha * chi1;
			psi1->real  += alpha_c.real *  chi1_c.real - alpha_c.imag * -chi1_c.imag;
			psi1->imag  += alpha_c.real * -chi1_c.imag + alpha_c.imag *  chi1_c.real;

			// zeta1 = zeta1 + beta  * chi1;
			zeta1->real += beta_c.real  *  chi1_c.real - beta_c.imag  * -chi1_c.imag;
			zeta1->imag += beta_c.real  * -chi1_c.imag + beta_c.imag  *  chi1_c.real;

			chi1  += inc_x;
			psi1  += inc_y;
			zeta1 += inc_z;
		}
	}
	else
	{
		if ( inc_x == 1 &&
		     inc_y == 1 &&
		     inc_z == 1 )
		{
			for ( i = 0; i < n; ++i )
			{
				x11v.v = _mm_load_pd( ( double* )chi1 );
				x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
				x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );

				acbc.v = alphav.v * x11v.v;
				bdad.v = alpharv.v * x12v.v;
				y1v.v = _mm_load_pd( ( double* )psi1 );
				y1v.v = y1v.v - _mm_addsub_pd( acbc.v, bdad.v );
				_mm_store_pd( ( double* )psi1, y1v.v );

				acbc.v = betav.v * x11v.v;
				bdad.v = betarv.v * x12v.v;
				z1v.v = _mm_load_pd( ( double* )zeta1 );
				z1v.v = z1v.v - _mm_addsub_pd( acbc.v, bdad.v );
				_mm_store_pd( ( double* )zeta1, z1v.v );

				chi1  += 1;
				psi1  += 1;
				zeta1 += 1;
			}
		}
		else
		{
			for ( i = 0; i < n; ++i )
			{
				x11v.v = _mm_load_pd( ( double* )chi1 );
				x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
				x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );

				acbc.v = alphav.v * x11v.v;
				bdad.v = alpharv.v * x12v.v;
				y1v.v = _mm_load_pd( ( double* )psi1 );
				y1v.v = y1v.v - _mm_addsub_pd( acbc.v, bdad.v );
				_mm_store_pd( ( double* )psi1, y1v.v );

				acbc.v = betav.v * x11v.v;
				bdad.v = betarv.v * x12v.v;
				z1v.v = _mm_load_pd( ( double* )zeta1 );
				z1v.v = z1v.v - _mm_addsub_pd( acbc.v, bdad.v );
				_mm_store_pd( ( double* )zeta1, z1v.v );

				chi1  += inc_x;
				psi1  += inc_y;
				zeta1 += inc_z;
			}
		}
	}
}
#elif BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict psi1;
	dcomplex* restrict zeta1;
	dcomplex  alpha_c;
	dcomplex  beta_c;
	int       i;

	chi1  = x;
	psi1  = y;
	zeta1 = z;

	alpha_c = *alpha;
	beta_c  = *beta;

	if ( bli_is_conj( conjx ) )
	{
		for ( i = 0; i < n; ++i )
		{
			dcomplex chi1_c = *chi1;

			// psi1  = psi1  - alpha * chi1;
			psi1->real  -= alpha_c.real *  chi1_c.real - alpha_c.imag * -chi1_c.imag;
			psi1->imag  -= alpha_c.real * -chi1_c.imag + alpha_c.imag *  chi1_c.real;

			// zeta1 = zeta1 - beta  * chi1;
			zeta1->real -= beta_c.real  *  chi1_c.real - beta_c.imag  * -chi1_c.imag;
			zeta1->imag -= beta_c.real  * -chi1_c.imag + beta_c.imag  *  chi1_c.real;

			chi1  += inc_x;
			psi1  += inc_y;
			zeta1 += inc_z;
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			dcomplex chi1_c  = *chi1;

			// psi1  = psi1  - alpha * chi1;
			psi1->real  -= alpha_c.real * chi1_c.real - alpha_c.imag * chi1_c.imag;
			psi1->imag  -= alpha_c.real * chi1_c.imag + alpha_c.imag * chi1_c.real;

			// zeta1 = zeta1 - beta  * chi1;
			zeta1->real -= beta_c.real  * chi1_c.real - beta_c.imag  * chi1_c.imag;
			zeta1->imag -= beta_c.real  * chi1_c.imag + beta_c.imag  * chi1_c.real;

			chi1  += inc_x;
			psi1  += inc_y;
			zeta1 += inc_z;
		}
	}
}
#endif

