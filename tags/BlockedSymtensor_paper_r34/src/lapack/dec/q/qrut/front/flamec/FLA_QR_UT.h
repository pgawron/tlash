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

#include "FLA_QR_UT_vars.h"

FLA_Error FLA_QR_UT( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_QR_UT_internal( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_copy_internal( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl );

FLA_Error FLA_QR_UT_create_T( FLA_Obj A, FLA_Obj* T );

FLA_Error FLA_QR_UT_recover_tau( FLA_Obj T, FLA_Obj tau );

FLA_Error FLA_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );

FLA_Error FLASH_QR_UT( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, FLA_Obj* A, FLA_Obj* TW );
FLA_Error FLASH_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );


FLA_Error FLA_QR_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_QR_UT_form_Q_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_QR_UT_form_Q_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_form_Q_ops_var1( int       m_A,
                                     int       n_AT,
                                     float*    buff_A, int rs_A, int cs_A,
                                     float*    buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opd_var1( int       m_A,
                                     int       n_AT,
                                     double*   buff_A, int rs_A, int cs_A,
                                     double*   buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opc_var1( int       m_A,
                                     int       n_AT,
                                     scomplex* buff_A, int rs_A, int cs_A,
                                     scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opz_var1( int       m_A,
                                     int       n_AT,
                                     dcomplex* buff_A, int rs_A, int cs_A,
                                     dcomplex* buff_T, int rs_T, int cs_T );
