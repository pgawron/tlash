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

void F77_ssytd2( char*     uplo,
                 int*      m,
                 float*    buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 float*    buff_t,
                 int*      info )
{
	F77_ssytrd( uplo,
	            m,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_t,
	            NULL, NULL,
	            info );
}

void F77_dsytd2( char*     uplo,
                 int*      m,
                 double*   buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 double*   buff_t,
                 int*      info )
{
	F77_dsytrd( uplo,
	            m,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_t,
	            NULL, NULL,
	            info );
}

void F77_chetd2( char*     uplo,
                 int*      m,
                 scomplex* buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 scomplex* buff_t,
                 int*      info )
{
	F77_chetrd( uplo,
	            m,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_t,
	            NULL, NULL,
	            info );
}

void F77_zhetd2( char*     uplo,
                 int*      m,
                 dcomplex* buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 dcomplex* buff_t,
                 int*      info )
{
	F77_zhetrd( uplo,
	            m,
	            buff_A, ldim_A,
	            buff_d,
	            buff_e,
	            buff_t,
	            NULL, NULL,
	            info );
}

#endif
