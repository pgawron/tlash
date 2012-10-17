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

FLA_Error FLA_Tridiag_UT_l_blk_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var1( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_blk_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_blf_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var2( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_blk_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_blf_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_unb_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_unb_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_l_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var1( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var1( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var1( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var1( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_T, int rs_T, int cs_T );

FLA_Error FLA_Tridiag_UT_l_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var2( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var2( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var2( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var2( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_T, int rs_T, int cs_T );

FLA_Error FLA_Tridiag_UT_l_opt_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_opt_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ops_var3( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_Z, int rs_Z, int cs_Z, 
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opd_var3( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_Z, int rs_Z, int cs_Z, 
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opc_var3( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_Z, int rs_Z, int cs_Z, 
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_opz_var3( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                          dcomplex* buff_T, int rs_T, int cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var1( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var1( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var1( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var1( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_T, int rs_T, int cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var2( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var2( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var2( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var2( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_T, int rs_T, int cs_T );

FLA_Error FLA_Tridiag_UT_l_ofu_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofu_var3( FLA_Obj A, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_l_step_ofs_var3( int m_A,
                                          int m_T,
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_Z, int rs_Z, int cs_Z,
                                          float* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofd_var3( int m_A,
                                          int m_T,
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_Z, int rs_Z, int cs_Z,
                                          double* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofc_var3( int m_A,
                                          int m_T,
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_Z, int rs_Z, int cs_Z,
                                          scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_l_step_ofz_var3( int m_A,
                                          int m_T,
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_Z, int rs_Z, int cs_Z,
                                          dcomplex* buff_T, int rs_T, int cs_T );

// --- Fused operations ---

FLA_Error FLA_Fused_Her2_Ax_l_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj z, FLA_Obj A, FLA_Obj x, FLA_Obj w );
FLA_Error FLA_Fused_Her2_Ax_l_ops_var1( int m_A,
                                        float* buff_alpha, 
                                        float* buff_u, int inc_u, 
                                        float* buff_z, int inc_z, 
                                        float* buff_A, int rs_A, int cs_A, 
                                        float* buff_x, int inc_x, 
                                        float* buff_w, int inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opd_var1( int m_A,
                                        double* buff_alpha, 
                                        double* buff_u, int inc_u, 
                                        double* buff_z, int inc_z, 
                                        double* buff_A, int rs_A, int cs_A, 
                                        double* buff_x, int inc_x, 
                                        double* buff_w, int inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opc_var1( int m_A,
                                        scomplex* buff_alpha, 
                                        scomplex* buff_u, int inc_u, 
                                        scomplex* buff_z, int inc_z, 
                                        scomplex* buff_A, int rs_A, int cs_A, 
                                        scomplex* buff_x, int inc_x, 
                                        scomplex* buff_w, int inc_w );
FLA_Error FLA_Fused_Her2_Ax_l_opz_var1( int m_A,
                                        dcomplex* buff_alpha, 
                                        dcomplex* buff_u, int inc_u, 
                                        dcomplex* buff_z, int inc_z, 
                                        dcomplex* buff_A, int rs_A, int cs_A, 
                                        dcomplex* buff_x, int inc_x, 
                                        dcomplex* buff_w, int inc_w );

FLA_Error FLA_Fused_UZhu_ZUhu_opt_var1( FLA_Obj delta, FLA_Obj U, FLA_Obj Z, FLA_Obj t, FLA_Obj u, FLA_Obj w );
FLA_Error FLA_Fused_UZhu_ZUhu_ops_var1( int m_U,
                                        int n_U,
                                        float* buff_delta, 
                                        float* buff_U, int rs_U, int cs_U, 
                                        float* buff_Z, int rs_Z, int cs_Z, 
                                        float* buff_t, int inc_t, 
                                        float* buff_u, int inc_u, 
                                        float* buff_w, int inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opd_var1( int m_U,
                                        int n_U,
                                        double* buff_delta, 
                                        double* buff_U, int rs_U, int cs_U, 
                                        double* buff_Z, int rs_Z, int cs_Z, 
                                        double* buff_t, int inc_t, 
                                        double* buff_u, int inc_u, 
                                        double* buff_w, int inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opc_var1( int m_U,
                                        int n_U,
                                        scomplex* buff_delta, 
                                        scomplex* buff_U, int rs_U, int cs_U, 
                                        scomplex* buff_Z, int rs_Z, int cs_Z, 
                                        scomplex* buff_t, int inc_t, 
                                        scomplex* buff_u, int inc_u, 
                                        scomplex* buff_w, int inc_w );
FLA_Error FLA_Fused_UZhu_ZUhu_opz_var1( int m_U,
                                        int n_U,
                                        dcomplex* buff_delta, 
                                        dcomplex* buff_U, int rs_U, int cs_U, 
                                        dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                        dcomplex* buff_t, int inc_t, 
                                        dcomplex* buff_u, int inc_u, 
                                        dcomplex* buff_w, int inc_w );
