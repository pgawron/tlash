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

#define MAC_Apply_G_2x2_ops( gamma, sigma, delta1, beta, epsilon1, delta2 ) \
{ \
	float g, s; \
	float d1, e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(delta1)    =  g * d1 + s * e1; \
	*(epsilon1)  = -s * d1 + g * e1; \
\
	*(beta)      = s * d2; \
	*(delta2)    = g * d2; \
}

#define MAC_Apply_G_2x2_opd( gamma, sigma, delta1, beta, epsilon1, delta2 ) \
{ \
	double g, s; \
	double d1, e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(delta1)    =  g * d1 + s * e1; \
	*(epsilon1)  = -s * d1 + g * e1; \
\
	*(beta)      = s * d2; \
	*(delta2)    = g * d2; \
}

