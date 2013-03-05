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

void F77_sgetrf( int*      m,
                 int*      n,
                 float*    buff_A, int* ldim_A,
                 int*      buff_p,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_FLOAT;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	FLA_Error    e_val;
	int          min_m_n    = min( *m, *n );
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, min_m_n, 1, &p );
	FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );

	e_val = FLA_LU_piv( A, p );
	FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = e_val + 1;
	else                        *info = 0;
}

void F77_dgetrf( int*      m,
                 int*      n,
                 double*   buff_A, int* ldim_A,
                 int*      buff_p,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_DOUBLE;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	FLA_Error    e_val;
	int          min_m_n    = min( *m, *n );
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, min_m_n, 1, &p );
	FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );

	e_val = FLA_LU_piv( A, p );
	FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = e_val + 1;
	else                        *info = 0;
}

void F77_cgetrf( int*      m,
                 int*      n,
                 scomplex* buff_A, int* ldim_A,
                 int*      buff_p,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_COMPLEX;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	FLA_Error    e_val;
	int          min_m_n    = min( *m, *n );
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, min_m_n, 1, &p );
	FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );

	e_val = FLA_LU_piv( A, p );
	FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = e_val + 1;
	else                        *info = 0;
}

void F77_zgetrf( int*      m,
                 int*      n,
                 dcomplex* buff_A, int* ldim_A,
                 int*      buff_p,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_DOUBLE_COMPLEX;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	FLA_Error    e_val;
	int          min_m_n    = min( *m, *n );
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, min_m_n, 1, &p );
	FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );

	e_val = FLA_LU_piv( A, p );
	FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = e_val + 1;
	else                        *info = 0;
}

#endif
