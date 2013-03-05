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

void F77_strsyl( char*     transa,
                 char*     transb,
                 int*      sgn,
                 int*      m,
                 int*      n,
                 float*    buff_A, int* ldim_A,
                 float*    buff_B, int* ldim_B,
                 float*    buff_C, int* ldim_C,
                 float*    scale,
                 int*      info )
{
	FLA_Datatype datatype       = FLA_FLOAT;
	FLA_Datatype datatype_scale = FLA_FLOAT;
	FLA_Trans    transa_fla;
	FLA_Trans    transb_fla;
	FLA_Obj      sgn_fla;
	FLA_Obj      A, B, C;
	FLA_Obj      scale_fla;
	FLA_Error    e_val;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_trans( transa, &transa_fla );
	FLA_Param_map_netlib_to_flame_trans( transb, &transb_fla );

	if ( *sgn == 1 ) sgn_fla = FLA_ONE;
	else             sgn_fla = FLA_MINUS_ONE;

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *n, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &C );
	FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );

	FLA_Obj_create_without_buffer( datatype_scale, 1, 1, &scale_fla );
	FLA_Obj_attach_buffer( &scale, 1, 1, &scale_fla );

	e_val = FLA_Sylv( transa_fla, transb_fla, sgn_fla, A, B, C, scale_fla );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &C );
	FLA_Obj_free_without_buffer( &scale_fla );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = 1;
	else                        *info = 0;
}

void F77_dtrsyl( char*     transa,
                 char*     transb,
                 int*      sgn,
                 int*      m,
                 int*      n,
                 double*   buff_A, int* ldim_A,
                 double*   buff_B, int* ldim_B,
                 double*   buff_C, int* ldim_C,
                 double*   scale,
                 int*      info )
{
	FLA_Datatype datatype       = FLA_DOUBLE;
	FLA_Datatype datatype_scale = FLA_DOUBLE;
	FLA_Trans    transa_fla;
	FLA_Trans    transb_fla;
	FLA_Obj      sgn_fla;
	FLA_Obj      A, B, C;
	FLA_Obj      scale_fla;
	FLA_Error    e_val;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_trans( transa, &transa_fla );
	FLA_Param_map_netlib_to_flame_trans( transb, &transb_fla );

	if ( *sgn == 1 ) sgn_fla = FLA_ONE;
	else             sgn_fla = FLA_MINUS_ONE;

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A,  1,*ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *n, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &C );
	FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );

	FLA_Obj_create_without_buffer( datatype_scale, 1, 1, &scale_fla );
	FLA_Obj_attach_buffer( &scale, 1, 1, &scale_fla );

	e_val = FLA_Sylv( transa_fla, transb_fla, sgn_fla, A, B, C, scale_fla );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &C );
	FLA_Obj_free_without_buffer( &scale_fla );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = 1;
	else                        *info = 0;
}

void F77_ctrsyl( char*     transa,
                 char*     transb,
                 int*      sgn,
                 int*      m,
                 int*      n,
                 scomplex* buff_A, int* ldim_A,
                 scomplex* buff_B, int* ldim_B,
                 scomplex* buff_C, int* ldim_C,
                 float*    scale,
                 int*      info )
{
	FLA_Datatype datatype       = FLA_COMPLEX;
	FLA_Datatype datatype_scale = FLA_FLOAT;
	FLA_Trans    transa_fla;
	FLA_Trans    transb_fla;
	FLA_Obj      sgn_fla;
	FLA_Obj      A, B, C;
	FLA_Obj      scale_fla;
	FLA_Error    e_val;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_trans( transa, &transa_fla );
	FLA_Param_map_netlib_to_flame_trans( transb, &transb_fla );

	if ( *sgn == 1 ) sgn_fla = FLA_ONE;
	else             sgn_fla = FLA_MINUS_ONE;

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *n, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &C );
	FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );

	FLA_Obj_create_without_buffer( datatype_scale, 1, 1, &scale_fla );
	FLA_Obj_attach_buffer( &scale, 1, 1, &scale_fla );

	e_val = FLA_Sylv( transa_fla, transb_fla, sgn_fla, A, B, C, scale_fla );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &C );
	FLA_Obj_free_without_buffer( &scale_fla );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = 1;
	else                        *info = 0;
}

void F77_ztrsyl( char*     transa,
                 char*     transb,
                 int*      sgn,
                 int*      m,
                 int*      n,
                 dcomplex* buff_A, int* ldim_A,
                 dcomplex* buff_B, int* ldim_B,
                 dcomplex* buff_C, int* ldim_C,
                 double*   scale,
                 int*      info )
{
	FLA_Datatype datatype       = FLA_DOUBLE_COMPLEX;
	FLA_Datatype datatype_scale = FLA_DOUBLE;
	FLA_Trans    transa_fla;
	FLA_Trans    transb_fla;
	FLA_Obj      sgn_fla;
	FLA_Obj      A, B, C;
	FLA_Obj      scale_fla;
	FLA_Error    e_val;
	FLA_Error    init_result;

	FLA_Init_safe( &init_result );

	FLA_Param_map_netlib_to_flame_trans( transa, &transa_fla );
	FLA_Param_map_netlib_to_flame_trans( transb, &transb_fla );

	if ( *sgn == 1 ) sgn_fla = FLA_ONE;
	else             sgn_fla = FLA_MINUS_ONE;

	FLA_Obj_create_without_buffer( datatype, *m, *m, &A );
	FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );

	FLA_Obj_create_without_buffer( datatype, *n, *n, &B );
	FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );

	FLA_Obj_create_without_buffer( datatype, *m, *n, &C );
	FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );

	FLA_Obj_create_without_buffer( datatype_scale, 1, 1, &scale_fla );
	FLA_Obj_attach_buffer( &scale, 1, 1, &scale_fla );

	e_val = FLA_Sylv( transa_fla, transb_fla, sgn_fla, A, B, C, scale_fla );

	FLA_Obj_free_without_buffer( &A );
	FLA_Obj_free_without_buffer( &B );
	FLA_Obj_free_without_buffer( &C );
	FLA_Obj_free_without_buffer( &scale_fla );

	FLA_Finalize_safe( init_result );

	if ( e_val != FLA_SUCCESS ) *info = 1;
	else                        *info = 0;
}

#endif
