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

void F77_sgehrd( int*      m,
                 int*      ilo,
                 int*      ihi,
                 float*    buff_A, int* ldim_A,
                 float*    buff_t,
                 float*    buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_FLOAT;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Hess_UT_create_T( A, &T );

	FLA_Hess_UT( A, T );
	FLA_Hess_UT_recover_tau( T, t );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_dgehrd( int*      m,
                 int*      ilo,
                 int*      ihi,
                 double*   buff_A, int* ldim_A,
                 double*   buff_t,
                 double*   buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Hess_UT_create_T( A, &T );

	FLA_Hess_UT( A, T );
	FLA_Hess_UT_recover_tau( T, t );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_cgehrd( int*      m,
                 int*      ilo,
                 int*      ihi,
                 scomplex* buff_A, int* ldim_A,
                 scomplex* buff_t,
                 scomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_COMPLEX;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Hess_UT_create_T( A, &T );

	FLA_Hess_UT( A, T );
	FLA_Hess_UT_recover_tau( T, t );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_zgehrd( int*      m,
                 int*      ilo,
                 int*      ihi,
                 dcomplex* buff_A, int* ldim_A,
                 dcomplex* buff_t,
                 dcomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE_COMPLEX;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Hess_UT_create_T( A, &T );

	FLA_Hess_UT( A, T );
	FLA_Hess_UT_recover_tau( T, t );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

#endif