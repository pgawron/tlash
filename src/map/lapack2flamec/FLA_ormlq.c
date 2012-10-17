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

void F77_sormlq( char*     side,
                 char*     trans,
                 int*      m,
                 int*      n,
                 int*      k,
                 float*    buff_A, int* ldim_A,
                 float*    buff_t,
                 float*    buff_B, int* ldim_B,
                 float*    buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_FLOAT;
	FLA_Direct   direct_fla = FLA_FORWARD;
	FLA_Store    storev_fla = FLA_ROWWISE;
	FLA_Side     side_fla;
	FLA_Trans    trans_fla;
	FLA_Obj      A, t, B;
	FLA_Obj      T, W;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_side( side, &side_fla );
	FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );

	if      ( trans_fla == FLA_CONJ_TRANSPOSE ) trans_fla = FLA_NO_TRANSPOSE;
	else if ( trans_fla == FLA_NO_TRANSPOSE   ) trans_fla = FLA_CONJ_TRANSPOSE;

	FLA_Obj_create_without_buffer( datatype, *m, *k, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *k, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, *k, &t );

	FLA_LQ_UT_create_T( A, &T );
	FLA_Obj_create( datatype, FLA_Obj_length( T ), FLA_Obj_width( B ), 0, 0, &W );

	FLA_Accum_T_UT( direct_fla, storev_fla, A, t, T );

	FLA_Apply_Q_UT( side_fla, trans_fla, direct_fla, storev_fla, A, T, W, B );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_dormlq( char*     side,
                 char*     trans,
                 int*      m,
                 int*      n,
                 int*      k,
                 double*   buff_A, int* ldim_A,
                 double*   buff_t,
                 double*   buff_B, int* ldim_B,
                 double*   buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_DOUBLE;
	FLA_Direct   direct_fla = FLA_FORWARD;
	FLA_Store    storev_fla = FLA_ROWWISE;
	FLA_Side     side_fla;
	FLA_Trans    trans_fla;
	FLA_Obj      A, t, B;
	FLA_Obj      T, W;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_side( side, &side_fla );
	FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );

	if      ( trans_fla == FLA_CONJ_TRANSPOSE ) trans_fla = FLA_NO_TRANSPOSE;
	else if ( trans_fla == FLA_NO_TRANSPOSE   ) trans_fla = FLA_CONJ_TRANSPOSE;

	FLA_Obj_create_without_buffer( datatype, *m, *k, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *k, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, *k, &t );

	FLA_LQ_UT_create_T( A, &T );
	FLA_Obj_create( datatype, FLA_Obj_length( T ), FLA_Obj_width( B ), 0, 0, &W );

	FLA_Accum_T_UT( direct_fla, storev_fla, A, t, T );

	FLA_Apply_Q_UT( side_fla, trans_fla, direct_fla, storev_fla, A, T, W, B );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_cunmlq( char*     side,
                 char*     trans,
                 int*      m,
                 int*      n,
                 int*      k,
                 scomplex* buff_A, int* ldim_A,
                 scomplex* buff_t,
                 scomplex* buff_B, int* ldim_B,
                 scomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_COMPLEX;
	FLA_Direct   direct_fla = FLA_FORWARD;
	FLA_Store    storev_fla = FLA_ROWWISE;
	FLA_Side     side_fla;
	FLA_Trans    trans_fla;
	FLA_Obj      A, t, B;
	FLA_Obj      T, W;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_side( side, &side_fla );
	FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );

	if      ( trans_fla == FLA_CONJ_TRANSPOSE ) trans_fla = FLA_NO_TRANSPOSE;
	else if ( trans_fla == FLA_NO_TRANSPOSE   ) trans_fla = FLA_CONJ_TRANSPOSE;

	FLA_Obj_create_without_buffer( datatype, *m, *k, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *k, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, *k, &t );

	FLA_LQ_UT_create_T( A, &T );
	FLA_Obj_create( datatype, FLA_Obj_length( T ), FLA_Obj_width( B ), 0, 0, &W );

	FLA_Accum_T_UT( direct_fla, storev_fla, A, t, T );

	FLA_Apply_Q_UT( side_fla, trans_fla, direct_fla, storev_fla, A, T, W, B );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

void F77_zunmlq( char*     side,
                 char*     trans,
                 int*      m,
                 int*      n,
                 int*      k,
                 dcomplex* buff_A, int* ldim_A,
                 dcomplex* buff_t,
                 dcomplex* buff_B, int* ldim_B,
                 dcomplex* buff_w, int* lwork,
                 int*      info )
{
	FLA_Datatype datatype   = FLA_DOUBLE_COMPLEX;
	FLA_Direct   direct_fla = FLA_FORWARD;
	FLA_Store    storev_fla = FLA_ROWWISE;
	FLA_Side     side_fla;
	FLA_Trans    trans_fla;
	FLA_Obj      A, t, B;
	FLA_Obj      T, W;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_side( side, &side_fla );
	FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );

	if      ( trans_fla == FLA_CONJ_TRANSPOSE ) trans_fla = FLA_NO_TRANSPOSE;
	else if ( trans_fla == FLA_NO_TRANSPOSE   ) trans_fla = FLA_CONJ_TRANSPOSE;

	FLA_Obj_create_without_buffer( datatype, *m, *k, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *k, 1, &t );
	FLA_Obj_attach_buffer( buff_t, 1, *k, &t );

	FLA_LQ_UT_create_T( A, &T );
	FLA_Obj_create( datatype, FLA_Obj_length( T ), FLA_Obj_width( B ), 0, 0, &W );

	FLA_Accum_T_UT( direct_fla, storev_fla, A, t, T );

	FLA_Apply_Q_UT( side_fla, trans_fla, direct_fla, storev_fla, A, T, W, B );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &t );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );

	FLA_Finalize_safe( init_result );

	*info = 0;
}

#endif
