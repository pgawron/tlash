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

FLA_Error FLA_Eig_gest_nl_blk_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl_blk_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl_blk_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl_blk_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl_blk_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );

FLA_Error FLA_Eig_gest_nl_unb_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_unb_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_unb_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B );

FLA_Error FLA_Eig_gest_nl_opt_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_ops_var1( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y,
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opd_var1( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y,
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opc_var1( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y,
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opz_var1( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y,
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nl_opt_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_ops_var2( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y,
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opd_var2( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y,
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opc_var2( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y,
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opz_var2( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y,
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nl_opt_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_ops_var3( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y,
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opd_var3( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y,
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opc_var3( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y,
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opz_var3( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y,
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nl_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_ops_var4( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y, 
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opd_var4( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y, 
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opc_var4( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opz_var4( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nl_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_ops_var5( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y, 
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opd_var5( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y, 
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opc_var5( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nl_opz_var5( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B );

