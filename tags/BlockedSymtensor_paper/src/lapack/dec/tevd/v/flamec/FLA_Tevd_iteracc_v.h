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

// --- FLA_Tevd_iteracc_v_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_iteracc_v_ops_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );
FLA_Error FLA_Tevd_iteracc_v_opd_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );

FLA_Error FLA_Tevd_iteracc_v_ops_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       float*    buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       float*    buff_pu, int inc_pu,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );
FLA_Error FLA_Tevd_iteracc_v_opd_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       double*   buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       double*   buff_pu, int inc_pu,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );

