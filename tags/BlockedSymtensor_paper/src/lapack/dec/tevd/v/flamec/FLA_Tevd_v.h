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

#include "FLA_Tevd_iteracc_v.h"
#include "FLA_Tevd_eigval_v.h"
#include "FLA_Tevd_francis_v.h"

// --- FLA_Tevd_compute_scaling() ----------------------------------------------

FLA_Error FLA_Tevd_compute_scaling_ops( int       m_A,
                                        float*    buff_d, int inc_d, 
                                        float*    buff_e, int inc_e,
                                        float*    sigma );
FLA_Error FLA_Tevd_compute_scaling_opd( int       m_A,
                                        double*   buff_d, int inc_d, 
                                        double*   buff_e, int inc_e,
                                        double*   sigma );

// --- FLA_Tevd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Tevd_find_submatrix_ops( int       m_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( int       m_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );

// --- FLA_Tevd_find_perfshift() -----------------------------------------------

FLA_Error FLA_Tevd_find_perfshift_ops( int       m_d,
                                       int       m_l,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e, 
                                       float*    buff_l, int inc_l, 
                                       int*      buff_lstat, int inc_lstat, 
                                       float*    buff_pu, int inc_pu, 
                                       int*      ij_shift );
FLA_Error FLA_Tevd_find_perfshift_opd( int       m_d,
                                       int       m_l,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e, 
                                       double*   buff_l, int inc_l, 
                                       int*      buff_lstat, int inc_lstat, 
                                       double*   buff_pu, int inc_pu, 
                                       int*      ij_shift );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( int       m_A,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( int       m_A,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 double*   norm );

// --- FLA_Tevd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

// --- FLA_Tevd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               float*    buff_W, int rs_W, int cs_W,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               double*   buff_W, int rs_W, int cs_W,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               scomplex* buff_W, int rs_W, int cs_W,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               dcomplex* buff_W, int rs_W, int cs_W,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

// --- FLA_Tevd_v_opt_var3() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var3( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj ls, FLA_Obj pu, FLA_Obj G, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

// --- FLA_Tevd_v_opt_var4() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var4( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj ls, FLA_Obj pu, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var4( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               float*    buff_W, int rs_W, int cs_W,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var4( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               double*   buff_W, int rs_W, int cs_W,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var4( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               scomplex* buff_W, int rs_W, int cs_W,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var4( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
			                   int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               dcomplex* buff_W, int rs_W, int cs_W,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

