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

// --- MAC_Bsvd_sinval_is_converged() ------------------------------------------

#define MAC_Bsvd_sinval_is_converged_ops( tol, d1, e1 ) \
	fabsf( (e1) ) <= fabsf( (tol) * (d1) )

#define MAC_Bsvd_sinval_is_converged_opd( tol, d1, e1 ) \
	fabs(  (e1) ) <= fabs(  (tol) * (d1) )

// --- FLA_Bsvd_sinval_v_opt_var1() --------------------------------------------

FLA_Error FLA_Bsvd_sinval_v_opt_var1( FLA_Obj tol, FLA_Obj thresh, FLA_Obj G, FLA_Obj H, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Bsvd_sinval_v_ops_var1( int       m_A,
                                      int       n_GH,
                                      int       n_iter_allowed,
                                      float     tol, 
                                      float     thresh, 
                                      scomplex* buff_G, int rs_G, int cs_G,
                                      scomplex* buff_H, int rs_H, int cs_H,
                                      float*    buff_d, int inc_d, 
                                      float*    buff_e, int inc_e,
                                      int*      n_iter );
FLA_Error FLA_Bsvd_sinval_v_opd_var1( int       m_A,
                                      int       n_GH,
                                      int       n_iter_allowed,
                                      double    tol, 
                                      double    thresh, 
                                      dcomplex* buff_G, int rs_G, int cs_G,
                                      dcomplex* buff_H, int rs_H, int cs_H,
                                      double*   buff_d, int inc_d, 
                                      double*   buff_e, int inc_e,
                                      int*      n_iter );

