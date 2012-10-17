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

void F77_slaswp( int*      n,
                 float*    buff_A, int* ldim_A,
                 int*      k1,
                 int*      k2,
                 int*      buff_p, int* inc_p )
{
	FLA_Datatype datatype   = FLA_FLOAT;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	int          n_elem_p   = *k2 - *k1 + 1;
	int*         p_begin    = buff_p + *k1 - 1;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *ldim_A, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, n_elem_p, 1, &p );
	FLA_Obj_attach_buffer( p_begin, 1, n_elem_p, &p );

	FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );
	FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );
}

void F77_dlaswp( int*      n,
                 double*   buff_A, int* ldim_A,
                 int*      k1,
                 int*      k2,
                 int*      buff_p, int* inc_p )
{
	FLA_Datatype datatype   = FLA_DOUBLE;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	int          n_elem_p   = *k2 - *k1 + 1;
	int*         p_begin    = buff_p + *k1 - 1;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *ldim_A, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, n_elem_p, 1, &p );
	FLA_Obj_attach_buffer( p_begin, 1, n_elem_p, &p );

	FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );
	FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );
}

void F77_claswp( int*      n,
                 scomplex* buff_A, int* ldim_A,
                 int*      k1,
                 int*      k2,
                 int*      buff_p, int* inc_p )
{
	FLA_Datatype datatype   = FLA_COMPLEX;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	int          n_elem_p   = *k2 - *k1 + 1;
	int*         p_begin    = buff_p + *k1 - 1;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *ldim_A, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, n_elem_p, 1, &p );
	FLA_Obj_attach_buffer( p_begin, 1, n_elem_p, &p );

	FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );
	FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );
}

void F77_zlaswp( int*      n,
                 dcomplex* buff_A, int* ldim_A,
                 int*      k1,
                 int*      k2,
                 int*      buff_p, int* inc_p )
{
	FLA_Datatype datatype   = FLA_DOUBLE_COMPLEX;
	FLA_Datatype datatype_p = FLA_INT;
	FLA_Obj      A, p;
	int          n_elem_p   = *k2 - *k1 + 1;
	int*         p_begin    = buff_p + *k1 - 1;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Obj_create_without_buffer( datatype, *ldim_A, *n, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype_p, n_elem_p, 1, &p );
	FLA_Obj_attach_buffer( p_begin, 1, n_elem_p, &p );

	FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );
	FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &p );

	FLA_Finalize_safe( init_result );
}

#endif
