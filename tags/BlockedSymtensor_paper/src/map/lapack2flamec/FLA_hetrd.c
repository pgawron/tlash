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

void F77_ssytrd( char*     uplo,
                 int*      m,
                 float*    buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 float*    buff_t,
                 float*    buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_FLOAT;
	FLA_Datatype dtype_re = FLA_FLOAT;
	FLA_Uplo     uplo_fla;
	dim_t        m_d      = *m;
	dim_t        m_e      = *m - 1;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, d, e, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Tridiag_UT_create_T( A, &T );

	FLA_Tridiag_UT( uplo_fla, A, T );
	FLA_Tridiag_UT_recover_tau( T, t );
	FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_dsytrd( char*     uplo,
                 int*      m,
                 double*   buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 double*   buff_t,
                 double*   buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE;
	FLA_Datatype dtype_re = FLA_DOUBLE;
	FLA_Uplo     uplo_fla;
	dim_t        m_d      = *m;
	dim_t        m_e      = *m - 1;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, d, e, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Tridiag_UT_create_T( A, &T );

	FLA_Tridiag_UT( uplo_fla, A, T );
	FLA_Tridiag_UT_recover_tau( T, t );
	FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_chetrd( char*     uplo,
                 int*      m,
                 scomplex* buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 scomplex* buff_t,
                 scomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_COMPLEX;
	FLA_Datatype dtype_re = FLA_FLOAT;
	FLA_Uplo     uplo_fla;
	dim_t        m_d      = *m;
	dim_t        m_e      = *m - 1;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, d, e, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Tridiag_UT_create_T( A, &T );

	FLA_Tridiag_UT( uplo_fla, A, T );
	FLA_Tridiag_UT_recover_tau( T, t );
	FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_zhetrd( char*     uplo,
                 int*      m,
                 dcomplex* buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 dcomplex* buff_t,
                 dcomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE_COMPLEX;
	FLA_Datatype dtype_re = FLA_DOUBLE;
	FLA_Uplo     uplo_fla;
	dim_t        m_d      = *m;
	dim_t        m_e      = *m - 1;
	dim_t        m_t      = *m - 1;
	FLA_Obj      A, d, e, t, T;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );

	FLA_Tridiag_UT_create_T( A, &T );

	FLA_Tridiag_UT( uplo_fla, A, T );
	FLA_Tridiag_UT_recover_tau( T, t );
	FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

#endif
