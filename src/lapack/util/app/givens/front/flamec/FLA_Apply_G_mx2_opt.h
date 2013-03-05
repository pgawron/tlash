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

#define MAC_Apply_G_mx2_ops( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	float             ga     = *gamma12; \
	float             si     = *sigma12; \
	float*  restrict  alpha1 = a1; \
	float*  restrict  alpha2 = a2; \
	float             temp1; \
	float             temp2; \
	int               i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 =  ga * temp1 + si * temp2; \
		*alpha2 = -si * temp1 + ga * temp2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}

#define MAC_Apply_G_mx2_opc( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	float              ga12   = *gamma12; \
	float              si12   = *sigma12; \
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex           temp1; \
	scomplex           temp2; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real =  ga12 * temp1.real + si12 * temp2.real; \
		alpha1->imag =  ga12 * temp1.imag + si12 * temp2.imag; \
\
		alpha2->real = -si12 * temp1.real + ga12 * temp2.real; \
		alpha2->imag = -si12 * temp1.imag + ga12 * temp2.imag; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}

#define MAC_Apply_G_mx2_opd( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	double            ga     = *gamma12; \
	double            si     = *sigma12; \
	double* restrict  alpha1 = a1; \
	double* restrict  alpha2 = a2; \
	double            temp1; \
	double            temp2; \
	int               i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 =  ga * temp1 + si * temp2; \
		*alpha2 = -si * temp1 + ga * temp2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}

#define MAC_Apply_G_mx2_opz( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{\
	double             ga12   = *gamma12; \
	double             si12   = *sigma12; \
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex           temp1; \
	dcomplex           temp2; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real =  ga12 * temp1.real + si12 * temp2.real; \
		alpha1->imag =  ga12 * temp1.imag + si12 * temp2.imag; \
\
		alpha2->real = -si12 * temp1.real + ga12 * temp2.real; \
		alpha2->imag = -si12 * temp1.imag + ga12 * temp2.imag; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}
