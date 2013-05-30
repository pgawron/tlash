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

FLA_Error FLA_Gemm_task( FLA_Trans transa, FLA_Trans transb, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( transa, transb, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_cc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ch_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_cn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ct_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_hn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_ht_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nh_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_nt_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tc_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_CONJ_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_th_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tn_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );
}

FLA_Error FLA_Gemm_tt_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
  return FLA_Gemm_external( FLA_TRANSPOSE, FLA_TRANSPOSE, alpha, A, B, beta, C );
}

