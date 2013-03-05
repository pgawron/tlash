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

     rho = conj(a) * x;
     w   = w + kappa * a;
*/

void bli_sdotaxpy( int       n,
                   float*    a, int inc_a,
                   float*    x, int inc_x,
                   float*    kappa,
                   float*    rho,
                   float*    w, int inc_w )
{
	bli_abort();
}


void bli_ddotaxpy( int       n,
                   double*   a, int inc_a,
                   double*   x, int inc_x,
                   double*   kappa,
                   double*   rho,
                   double*   w, int inc_w )
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
{
	double*   restrict alpha1;
	double*   restrict chi1;
	double*   restrict omega1;
	double             rho_c;
	int                i;

	int                n_pre;
	int                n_run;
	int                n_left;

	v2df_t    k1v, rho1v;
	v2df_t    a1v, x1v, w1v;
	v2df_t    a2v, x2v, w2v;
	
	if ( inc_a != 1 ||
	     inc_x != 1 ||
	     inc_w != 1 ) bli_abort();

	n_pre = 0;
	if ( ( unsigned long ) a % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) w % 16 == 0 ) bli_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	alpha1   = a;
	chi1     = x;
	omega1   = w;

	rho_c = 0.0;

	if ( n_pre == 1 )
	{
		double   kappa_c    = *kappa;
		double   alpha1_c   = *alpha1;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		rho_c += alpha1_c * chi1_c;
		omega1_c += kappa_c * alpha1_c;

		*omega1 = omega1_c;

		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	rho1v.v = _mm_setzero_pd();

	k1v.v = _mm_loaddup_pd( ( double* )kappa );

	for ( i = 0; i < n_run; ++i )
	{
		a1v.v = _mm_load_pd( ( double* )alpha1 );
		x1v.v = _mm_load_pd( ( double* )chi1 );
		w1v.v = _mm_load_pd( ( double* )omega1 );

		a2v.v = _mm_load_pd( ( double* )(alpha1 + 2) );
		x2v.v = _mm_load_pd( ( double* )(chi1 + 2) );
		w2v.v = _mm_load_pd( ( double* )(omega1 + 2) );

		rho1v.v += a1v.v * x1v.v;
		w1v.v += k1v.v * a1v.v;

		_mm_store_pd( ( double* )omega1, w1v.v );

		rho1v.v += a2v.v * x2v.v;
		w2v.v += k1v.v * a2v.v;

		_mm_store_pd( ( double* )(omega1 + 2), w2v.v );

		alpha1   += 4;
		chi1     += 4;
		omega1   += 4;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			double   kappa_c    = *kappa;
			double   alpha1_c   = *alpha1;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			rho_c += alpha1_c * chi1_c;
			omega1_c += kappa_c * alpha1_c;

			*omega1 = omega1_c;

			alpha1   += inc_a;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	rho_c += rho1v.d[0] + rho1v.d[1];

	*rho = rho_c;
}
#elif BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
{
	double*   restrict alpha1;
	double*   restrict chi1;
	double*   restrict omega1;
	double             kappa_c;
	double             rho_c;
	int                i;

	int                n_pre;
	int                n_run;
	int                n_left;
	
	if ( inc_a != 1 ||
	     inc_x != 1 ||
	     inc_w != 1 ) bli_abort();

	n_pre = 0;
	if ( ( unsigned long ) a % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) w % 16 == 0 ) bli_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	alpha1   = a;
	chi1     = x;
	omega1   = w;

	rho_c = 0.0;

	kappa_c = *kappa;

	if ( n_pre == 1 )
	{
		double   alpha1_c   = *alpha1;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		rho_c += alpha1_c * chi1_c;
		omega1_c += kappa_c * alpha1_c;

		*omega1 = omega1_c;

		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	for ( i = 0; i < n_run; ++i )
	{
		double   alpha1_c   = *alpha1;
		double   alpha2_c   = *(alpha1 + 1);
		double   chi1_c     = *chi1;
		double   chi2_c     = *(chi1 + 1);
		double   omega1_c   = *omega1;
		double   omega2_c   = *(omega1 + 1);

		// rho += conj(alpha1) * chi1;
		rho_c += alpha1_c * chi1_c;
		rho_c += alpha2_c * chi2_c;

		// omega1 += kappa * alpha1;
		omega1_c += kappa_c * alpha1_c;
		omega2_c += kappa_c * alpha2_c;

		*omega1       = omega1_c;
		*(omega1 + 1) = omega2_c;

		alpha1   += 2;
		chi1     += 2;
		omega1   += 2;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			double   alpha1_c   = *alpha1;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			rho_c += alpha1_c * chi1_c;
			omega1_c += kappa_c * alpha1_c;

			*omega1 = omega1_c;

			alpha1   += inc_a;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	*rho = rho_c;
}
#endif


void bli_cdotaxpy( int       n,
                   scomplex* a, int inc_a,
                   scomplex* x, int inc_x,
                   scomplex* kappa,
                   scomplex* rho,
                   scomplex* w, int inc_w )
{
	bli_abort();
}


void bli_zdotaxpy( int       n,
                   dcomplex* a, int inc_a,
                   dcomplex* x, int inc_x,
                   dcomplex* kappa,
                   dcomplex* rho,
                   dcomplex* w, int inc_w )
#if BLIS_VECTOR_INTRINSIC_TYPE == BLIS_SSE_INTRINSICS
{
	dcomplex* restrict alpha1;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	int                i;

	v2df_t    kappa1v, kappa1rv;
	v2df_t    rho1v;
	v2df_t    a11v, a12v;
	v2df_t    x1v, x1rv;
	v2df_t    w1v;
	v2df_t    acbc, bdad;
	v2df_t    adac, bcbd;

	alpha1   = a;
	chi1     = x;
	omega1   = w;

	if ( inc_a != 1 ||
	     inc_x != 1 ||
	     inc_w != 1 ) bli_abort();

	kappa1v.v  = _mm_load_pd( ( double* )kappa );
	kappa1rv.v = _mm_shuffle_pd( kappa1v.v, kappa1v.v, _MM_SHUFFLE2 (0,1) );

	rho1v.v = _mm_setzero_pd();

	for ( i = 0; i < n; ++i )
	{
		//alpha_c = *alpha1;
		a11v.v  = _mm_loaddup_pd( ( double* )&(alpha1->real) );
		a12v.v  = _mm_loaddup_pd( ( double* )&(alpha1->imag) );

		//rho_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		//rho_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;
		x1v.v  = _mm_load_pd( ( double* )chi1 );
		x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
		adac.v = a11v.v * x1rv.v;
		bcbd.v = a12v.v * x1v.v;
		rho1v.v = rho1v.v + _mm_addsub_pd( adac.v, bcbd.v );

		//omega_c = *omega1;
		w1v.v  = _mm_load_pd( ( double* )omega1 );

		//omega1_c.real += kappa_c.real * alpha1_c.real - kappa_c.imag * alpha1_c.imag;
		//omega1_c.imag += kappa_c.real * alpha1_c.imag + kappa_c.imag * alpha1_c.real;
		acbc.v = kappa1v.v  * a11v.v;
		bdad.v = kappa1rv.v * a12v.v;
		w1v.v += _mm_addsub_pd( acbc.v, bdad.v );

		//*omega1 = omega1_c;
		_mm_store_pd( ( double* )omega1, w1v.v );

		alpha1   += 1;
		chi1     += 1;
		omega1   += 1;
	}

	rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );

	//rho->real = rho_c.real;
	//rho->imag = rho_c.imag;
	_mm_store_pd( ( double* )rho, rho1v.v );
}
#elif BLIS_VECTOR_INTRINSIC_TYPE == BLIS_NO_INTRINSICS
{
	dcomplex* restrict alpha1;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	dcomplex           kappa_c;
	dcomplex           rho_c;
	int                i;

	alpha1   = a;
	chi1     = x;
	omega1   = w;

	rho_c.real = 0.0;
	rho_c.imag = 0.0;

	kappa_c = *kappa;

	for ( i = 0; i < n; ++i )
	{
		dcomplex alpha1_c   = *alpha1;
		dcomplex chi1_c     = *chi1;
		dcomplex omega1_c   = *omega1;

		// rho += conj(alpha1) * chi1;
		rho_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		rho_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;

		// omega1 += kappa * alpha1;
		omega1_c.real += kappa_c.real * alpha1_c.real - kappa_c.imag * alpha1_c.imag;
		omega1_c.imag += kappa_c.real * alpha1_c.imag + kappa_c.imag * alpha1_c.real;

		*omega1 = omega1_c;

		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	rho->real = rho_c.real;
	rho->imag = rho_c.imag;
}
#endif

