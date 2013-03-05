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


#if FLA_VECTOR_INTRINSIC_TYPE == FLA_NO_INTRINSICS

#define MAC_Apply_G_mx3b_ass MAC_Apply_G_mx3b_ops
#define MAC_Apply_G_mx3b_asd MAC_Apply_G_mx3b_opd
#define MAC_Apply_G_mx3b_asc MAC_Apply_G_mx3b_opc
#define MAC_Apply_G_mx3b_asz MAC_Apply_G_mx3b_opz

#elif FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS

#define MAC_Apply_G_mx3b_ass( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{\
	int              n_iter32 = m_A / ( 4 * 8 ); \
	int              n_left32 = m_A % ( 4 * 8 ); \
	int              n_iter4  = n_left32 / ( 4 * 1 ); \
	int              n_left   = n_left32 % ( 4 * 1 ); \
	int              i; \
\
	const int        step_a1 = inc_a1 * 4; \
	const int        step_a2 = inc_a2 * 4; \
	const int        step_a3 = inc_a3 * 4; \
\
	float*  restrict alpha1 = a1; \
	float*  restrict alpha2 = a2; \
	float*  restrict alpha3 = a3; \
\
	v4sf_t           a1v, a2v, a3v; \
	v4sf_t           g12v, s12v; \
	v4sf_t           g23v, s23v; \
	v4sf_t           t1v, t2v; \
\
	g12v.v = _mm_load1_ps( gamma12 ); \
	s12v.v = _mm_load1_ps( sigma12 ); \
	g23v.v = _mm_load1_ps( gamma23 ); \
	s23v.v = _mm_load1_ps( sigma23 ); \
\
	for ( i = 0; i < n_iter32; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	for ( i = 0; i < n_iter4; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
		float ga12 = *gamma12; \
		float si12 = *sigma12; \
		float ga23 = *gamma23; \
		float si23 = *sigma23; \
		float temp1; \
		float temp2; \
		float temp3; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
\
		alpha1 += 1; \
		alpha2 += 1; \
		alpha3 += 1; \
	} \
}

#define MAC_Apply_G_mx3b_asd( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{\
	int              n_iter16 = m_A / ( 2 * 8 ); \
	int              n_left16 = m_A % ( 2 * 8 ); \
	int              n_iter2  = n_left16 / ( 2 * 1 ); \
	int              n_left   = n_left16 % ( 2 * 1 ); \
	int              i; \
\
	const int        step_a1 = inc_a1 * 2; \
	const int        step_a2 = inc_a2 * 2; \
	const int        step_a3 = inc_a3 * 2; \
\
	double* restrict alpha1 = a1; \
	double* restrict alpha2 = a2; \
	double* restrict alpha3 = a3; \
\
	v2df_t           a1v, a2v, a3v; \
	v2df_t           g12v, s12v; \
	v2df_t           g23v, s23v; \
	v2df_t           t1v, t2v; \
\
	g12v.v = _mm_loaddup_pd( gamma12 ); \
	s12v.v = _mm_loaddup_pd( sigma12 ); \
	g23v.v = _mm_loaddup_pd( gamma23 ); \
	s23v.v = _mm_loaddup_pd( sigma23 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		double ga12 = *gamma12; \
		double si12 = *sigma12; \
		double ga23 = *gamma23; \
		double si23 = *sigma23; \
		double temp1; \
		double temp2; \
		double temp3; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
	} \
}

#define MAC_Apply_G_mx3b_asc( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{\
	int                n_iter16 = m_A / ( 2 * 8 ); \
	int                n_left16 = m_A % ( 2 * 8 ); \
	int                n_iter2  = n_left16 / ( 2 * 1 ); \
	int                n_left   = n_left16 % ( 2 * 1 ); \
	int                i; \
\
	const int          step_a1 = inc_a1 * 2; \
	const int          step_a2 = inc_a2 * 2; \
	const int          step_a3 = inc_a3 * 2; \
\
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex* restrict alpha3 = a3; \
\
	v4sf_t             a1v, a2v, a3v; \
	v4sf_t             g12v, s12v; \
	v4sf_t             g23v, s23v; \
	v4sf_t             t1v, t2v; \
\
	g12v.v = _mm_load1_ps( gamma12 ); \
	s12v.v = _mm_load1_ps( sigma12 ); \
	g23v.v = _mm_load1_ps( gamma23 ); \
	s23v.v = _mm_load1_ps( sigma23 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		float ga12 = *gamma12; \
		float si12 = *sigma12; \
		float ga23 = *gamma23; \
		float si23 = *sigma23; \
		scomplex temp1; \
		scomplex temp2; \
		scomplex temp3; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real = temp1.real * ga12 + temp2.real * si12; \
		alpha2->real = temp2.real * ga12 - temp1.real * si12; \
\
		alpha1->imag = temp1.imag * ga12 + temp2.imag * si12; \
		alpha2->imag = temp2.imag * ga12 - temp1.imag * si12; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23 + temp3.real * si23; \
		alpha3->real = temp3.real * ga23 - temp2.real * si23; \
\
		alpha2->imag = temp2.imag * ga23 + temp3.imag * si23; \
		alpha3->imag = temp3.imag * ga23 - temp2.imag * si23; \
	} \
}

#define MAC_Apply_G_mx3b_asz( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{\
	int                n_iter = m_A / 8; \
	int                n_left = m_A % 8; \
	int                i; \
\
	const int          step_a1 = inc_a1 * 1; \
	const int          step_a2 = inc_a2 * 1; \
	const int          step_a3 = inc_a3 * 1; \
\
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex* restrict alpha3 = a3; \
\
	v2df_t             a1v, a2v, a3v; \
	v2df_t             g12v, s12v; \
	v2df_t             g23v, s23v; \
	v2df_t             t1v, t2v; \
\
	g12v.v = _mm_loaddup_pd( gamma12 ); \
	s12v.v = _mm_loaddup_pd( sigma12 ); \
	g23v.v = _mm_loaddup_pd( gamma23 ); \
	s23v.v = _mm_loaddup_pd( sigma23 ); \
\
	for ( i = 0; i < n_iter; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
	} \
}

#endif
