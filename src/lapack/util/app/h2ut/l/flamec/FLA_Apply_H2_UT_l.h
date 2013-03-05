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

FLA_Error FLA_Apply_H2_UT_l_unb_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_l_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_l_ops_var1( int m_u2_A2,
                                      int n_a1t,
                                      float* tau,
                                      float* u2, int inc_u2,
                                      float* a1t, int inc_a1t,
                                      float* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opd_var1( int m_u2_A2,
                                      int n_a1t,
                                      double* tau,
                                      double* u2, int inc_u2,
                                      double* a1t, int inc_a1t,
                                      double* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opc_var1( int m_u2_A2,
                                      int n_a1t,
                                      scomplex* tau,
                                      scomplex* u2, int inc_u2,
                                      scomplex* a1t, int inc_a1t,
                                      scomplex* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opz_var1( int m_u2_A2,
                                      int n_a1t,
                                      dcomplex* tau,
                                      dcomplex* u2, int inc_u2,
                                      dcomplex* a1t, int inc_a1t,
                                      dcomplex* A2, int rs_A2, int cs_A2 );

