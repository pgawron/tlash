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

#ifdef FLA_ENABLE_LAPACK2FLAME

void F77_sgebd2( int*      m,
                 int*      n,
                 float*    buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 float*    buff_tu,
                 float*    buff_tv,
                 float*    buff_w,
                 int*      info )
{
	int max_m_n = max( *m, *n );

	F77_sgebrd( m,
	            n,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_tu,
	            buff_tv,
	            buff_w, &max_m_n,
	            info );
}

void F77_dgebd2( int*      m,
                 int*      n,
                 double*   buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 double*   buff_tu,
                 double*   buff_tv,
                 double*   buff_w,
                 int*      info )
{
	int max_m_n = max( *m, *n );

	F77_dgebrd( m,
	            n,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_tu,
	            buff_tv,
	            buff_w, &max_m_n,
	            info );
}

void F77_cgebd2( int*      m,
                 int*      n,
                 scomplex* buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 scomplex* buff_tu,
                 scomplex* buff_tv,
                 scomplex* buff_w,
                 int*      info )
{
	int max_m_n = max( *m, *n );

	F77_cgebrd( m,
	            n,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_tu,
	            buff_tv,
	            buff_w, &max_m_n,
	            info );
}

void F77_zgebd2( int*      m,
                 int*      n,
                 dcomplex* buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 dcomplex* buff_tu,
                 dcomplex* buff_tv,
                 dcomplex* buff_w,
                 int*      info )
{
	int max_m_n = max( *m, *n );

	F77_zgebrd( m,
	            n,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_tu,
	            buff_tv,
	            buff_w, &max_m_n,
	            info );
}

#endif
