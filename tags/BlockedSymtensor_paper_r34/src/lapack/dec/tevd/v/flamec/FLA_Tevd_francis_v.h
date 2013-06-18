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

// --- FLA_Tevd_francis_v_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_francis_v_opt_var1( FLA_Obj shift, FLA_Obj g, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tevd_francis_v_ops_var1( int       m_A,
                                       float*    buff_shift,
                                       scomplex* buff_g, int inc_g, 
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e ); 
FLA_Error FLA_Tevd_francis_v_opd_var1( int       m_A,
                                       double*   buff_shift,
                                       dcomplex* buff_g, int inc_g, 
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e ); 
