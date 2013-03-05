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

FLA_Error FLA_Apply_G_rb_opt_var1( FLA_Obj c, FLA_Obj s, FLA_Obj A );
FLA_Error FLA_Apply_G_rb_ops_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opd_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opc_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opz_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   dcomplex* buff_A, int rs_A, int cs_A );

