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

FLA_Error FLA_Trmm_task( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( side, uplo, trans, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_llt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_luc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_luh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lun_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_lut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rln_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rlt_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_ruc_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_ruh_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_run_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, diag, alpha, A, B );
}

FLA_Error FLA_Trmm_rut_task( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trmm_t* cntl )
{
  return FLA_Trmm_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, FLA_TRANSPOSE, diag, alpha, A, B );
}

