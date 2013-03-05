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

FLA_Error FLA_Givens2( FLA_Obj chi_1, FLA_Obj chi_2, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj chi_1_new );
FLA_Error FLA_Givens2_ops( float*  chi_1,
                           float*  chi_2,
                           float*  gamma,
                           float*  sigma,
                           float*  chi_1_new );
FLA_Error FLA_Givens2_opd( double* chi_1,
                           double* chi_2,
                           double* gamma,
                           double* sigma,
                           double* chi_1_new );
#define MAC_Givens2_ops( chi_1, chi_2, gamma, sigma, chi_1_new ) \
{ \
	float  chi_1_orig = *(chi_1); \
	float  chi_2_orig = *(chi_2); \
	float  g, s; \
	float  norm_x; \
\
	norm_x = ( float  ) sqrt( ( float ) ( chi_1_orig * chi_1_orig + \
	                                      chi_2_orig * chi_2_orig ) ); \
\
	g = chi_1_orig / norm_x; \
	s = chi_2_orig / norm_x; \
\
	if ( fabs( chi_1_orig ) > fabs( chi_2_orig ) && g < 0.0F ) \
	{ \
		g      = -g; \
		s      = -s; \
		norm_x = -norm_x; \
	} \
\
	*(gamma)     = g; \
	*(sigma)     = s; \
	*(chi_1_new) = norm_x; \
\
}

#define MAC_Givens2_opd( chi_1, chi_2, gamma, sigma, chi_1_new ) \
{ \
	double chi_1_orig = *(chi_1); \
	double chi_2_orig = *(chi_2); \
	double g, s; \
	double norm_x; \
\
	norm_x = ( double ) sqrt( chi_1_orig * chi_1_orig + \
	                          chi_2_orig * chi_2_orig ); \
\
	g = chi_1_orig / norm_x; \
	s = chi_2_orig / norm_x; \
\
	if ( fabs( chi_1_orig ) > fabs( chi_2_orig ) && g < 0.0 ) \
	{ \
		g      = -g; \
		s      = -s; \
		norm_x = -norm_x; \
	} \
\
	*(gamma)     = g; \
	*(sigma)     = s; \
	*(chi_1_new) = norm_x; \
\
}

