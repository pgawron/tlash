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

FLA_Error FLA_Sylv_hn_blk_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_blk_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );

FLA_Error FLA_Sylv_hn_opt_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_opt_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

FLA_Error FLA_Sylv_hn_ops_var1( float sgn,
                                int m_C,
                                int n_C,
                                float* buff_A, int rs_A, int cs_A,
                                float* buff_B, int rs_B, int cs_B,
                                float* buff_C, int rs_C, int cs_C,
                                float* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hn_opd_var1( double sgn,
                                int m_C,
                                int n_C,
                                double* buff_A, int rs_A, int cs_A,
                                double* buff_B, int rs_B, int cs_B,
                                double* buff_C, int rs_C, int cs_C,
                                double* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hn_opc_var1( float sgn,
                                int m_C,
                                int n_C,
                                scomplex* buff_A, int rs_A, int cs_A,
                                scomplex* buff_B, int rs_B, int cs_B,
                                scomplex* buff_C, int rs_C, int cs_C,
                                scomplex* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hn_opz_var1( double sgn,
                                int m_C,
                                int n_C,
                                dcomplex* buff_A, int rs_A, int cs_A,
                                dcomplex* buff_B, int rs_B, int cs_B,
                                dcomplex* buff_C, int rs_C, int cs_C,
                                dcomplex* buff_scale,
                                int* info );
