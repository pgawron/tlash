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

void F77_sgebrd( int*      m,
                 int*      n,
                 float*    buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 float*    buff_tu,
                 float*    buff_tv,
                 float*    buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_FLOAT;
	FLA_Datatype dtype_re = FLA_FLOAT;
    dim_t        min_m_n  = min( *m, *n );
	dim_t        m_d      = min_m_n;
	dim_t        m_e      = min_m_n - 1;
	dim_t        m_tu     = min_m_n;
	dim_t        m_tv     = min_m_n;
	FLA_Obj      A, d, e, tu, tv, TU, TV;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_tu, 1, &tu );
	FLA_Obj_attach_buffer( buff_tu, 1, m_tu, &tu );

	FLA_Obj_create_without_buffer( datatype, m_tv, 1, &tv );
	FLA_Obj_attach_buffer( buff_tv, 1, m_tv, &tv );

	FLA_Bidiag_UT_create_T( A, &TU, &TV );

	FLA_Bidiag_UT( A, TU, TV );
	FLA_Bidiag_UT_recover_tau( TU, TV, tu, tv );
	FLA_Bidiag_UT_extract_diagonals( A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &tu );
	FLA_Obj_free_without_buffer( &tv );
	FLA_Obj_free( &TU );
	FLA_Obj_free( &TV );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_dgebrd( int*      m,
                 int*      n,
                 double*   buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 double*   buff_tu,
                 double*   buff_tv,
                 double*   buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE;
	FLA_Datatype dtype_re = FLA_DOUBLE;
    dim_t        min_m_n  = min( *m, *n );
	dim_t        m_d      = min_m_n;
	dim_t        m_e      = min_m_n - 1;
	dim_t        m_tu     = min_m_n;
	dim_t        m_tv     = min_m_n;
	FLA_Obj      A, d, e, tu, tv, TU, TV;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_tu, 1, &tu );
	FLA_Obj_attach_buffer( buff_tu, 1, m_tu, &tu );

	FLA_Obj_create_without_buffer( datatype, m_tv, 1, &tv );
	FLA_Obj_attach_buffer( buff_tv, 1, m_tv, &tv );

	FLA_Bidiag_UT_create_T( A, &TU, &TV );

	FLA_Bidiag_UT( A, TU, TV );
	FLA_Bidiag_UT_recover_tau( TU, TV, tu, tv );
	FLA_Bidiag_UT_extract_diagonals( A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &tu );
	FLA_Obj_free_without_buffer( &tv );
	FLA_Obj_free( &TU );
	FLA_Obj_free( &TV );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_cgebrd( int*      m,
                 int*      n,
                 scomplex* buff_A, int* ldim_A,
                 float*    buff_d,
                 float*    buff_e,
                 scomplex* buff_tu,
                 scomplex* buff_tv,
                 scomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_COMPLEX;
	FLA_Datatype dtype_re = FLA_FLOAT;
    dim_t        min_m_n  = min( *m, *n );
	dim_t        m_d      = min_m_n;
	dim_t        m_e      = min_m_n - 1;
	dim_t        m_tu     = min_m_n;
	dim_t        m_tv     = min_m_n;
	FLA_Obj      A, d, e, tu, tv, TU, TV;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_tu, 1, &tu );
	FLA_Obj_attach_buffer( buff_tu, 1, m_tu, &tu );

	FLA_Obj_create_without_buffer( datatype, m_tv, 1, &tv );
	FLA_Obj_attach_buffer( buff_tv, 1, m_tv, &tv );

	FLA_Bidiag_UT_create_T( A, &TU, &TV );

	FLA_Bidiag_UT( A, TU, TV );
	FLA_Bidiag_UT_recover_tau( TU, TV, tu, tv );
	FLA_Bidiag_UT_extract_diagonals( A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &tu );
	FLA_Obj_free_without_buffer( &tv );
	FLA_Obj_free( &TU );
	FLA_Obj_free( &TV );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_zgebrd( int*      m,
                 int*      n,
                 dcomplex* buff_A, int* ldim_A,
                 double*   buff_d,
                 double*   buff_e,
                 dcomplex* buff_tu,
                 dcomplex* buff_tv,
                 dcomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype = FLA_DOUBLE_COMPLEX;
	FLA_Datatype dtype_re = FLA_DOUBLE;
    dim_t        min_m_n  = min( *m, *n );
	dim_t        m_d      = min_m_n;
	dim_t        m_e      = min_m_n - 1;
	dim_t        m_tu     = min_m_n;
	dim_t        m_tv     = min_m_n;
	FLA_Obj      A, d, e, tu, tv, TU, TV;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );
	FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );

	FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );
	FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );

	FLA_Obj_create_without_buffer( datatype, m_tu, 1, &tu );
	FLA_Obj_attach_buffer( buff_tu, 1, m_tu, &tu );

	FLA_Obj_create_without_buffer( datatype, m_tv, 1, &tv );
	FLA_Obj_attach_buffer( buff_tv, 1, m_tv, &tv );

	FLA_Bidiag_UT_create_T( A, &TU, &TV );

	FLA_Bidiag_UT( A, TU, TV );
	FLA_Bidiag_UT_recover_tau( TU, TV, tu, tv );
	FLA_Bidiag_UT_extract_diagonals( A, d, e );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &d );
	FLA_Obj_free_without_buffer( &e );
	FLA_Obj_free_without_buffer( &tu );
	FLA_Obj_free_without_buffer( &tv );
	FLA_Obj_free( &TU );
	FLA_Obj_free( &TV );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

#endif
