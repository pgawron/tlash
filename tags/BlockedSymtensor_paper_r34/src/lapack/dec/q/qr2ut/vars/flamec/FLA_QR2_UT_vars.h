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

FLA_Error FLA_QR2_UT_blk_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );
FLA_Error FLA_QR2_UT_blk_var2( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );

FLA_Error FLA_QR2_UT_unb_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T );

FLA_Error FLA_QR2_UT_opt_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T );

FLA_Error FLA_QR2_UT_ops_var1( int m_UT,
                               int m_D,
                               float* U, int rs_U, int cs_U,
                               float* D, int rs_D, int cs_D,
                               float* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opd_var1( int m_UT,
                               int m_D,
                               double* U, int rs_U, int cs_U,
                               double* D, int rs_D, int cs_D,
                               double* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opc_var1( int m_UT,
                               int m_D,
                               scomplex* U, int rs_U, int cs_U,
                               scomplex* D, int rs_D, int cs_D,
                               scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opz_var1( int m_UT,
                               int m_D,
                               dcomplex* U, int rs_U, int cs_U,
                               dcomplex* D, int rs_D, int cs_D,
                               dcomplex* T, int rs_T, int cs_T );
