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

#include "FLA_Tridiag_UT_l.h"
//#include "FLA_Tridiag_UT_u.h"

FLA_Error FLA_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );

FLA_Error FLA_Tridiag_UT_l( FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );
FLA_Error FLA_Tridiag_UT_u( FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );

FLA_Error FLA_Tridiag_UT_create_T( FLA_Obj A, FLA_Obj* T );

FLA_Error FLA_Tridiag_UT_recover_tau( FLA_Obj T, FLA_Obj t );

FLA_Error FLA_Tridiag_UT_extract_diagonals( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tridiag_UT_l_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tridiag_UT_u_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Tridiag_UT_realify( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_l_realify_unb( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_l_realify_opt( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_u_realify_unb( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_u_realify_opt( FLA_Obj A, FLA_Obj d );


FLA_Error FLA_Tridiag_UT_shift_U( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Tridiag_UT_shift_U_l_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Tridiag_UT_form_Q( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_form_Q_l_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_Tridiag_UT_form_Q_u_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_Tridiag_UT_form_Q_l_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_form_Q_l_ops_var1( int       m_A,
                                            int       n_AT,
                                            float*    buff_A, int rs_A, int cs_A,
                                            float*    buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opd_var1( int       m_A,
                                            int       n_AT,
                                            double*   buff_A, int rs_A, int cs_A,
                                            double*   buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opc_var1( int       m_A,
                                            int       n_AT,
                                            scomplex* buff_A, int rs_A, int cs_A,
                                            scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opz_var1( int       m_A,
                                            int       n_AT,
                                            dcomplex* buff_A, int rs_A, int cs_A,
                                            dcomplex* buff_T, int rs_T, int cs_T );
