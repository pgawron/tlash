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

#include "FLA_Trsm_llc.h"
#include "FLA_Trsm_llh.h"
#include "FLA_Trsm_lln.h"
#include "FLA_Trsm_llt.h"
#include "FLA_Trsm_luc.h"
#include "FLA_Trsm_luh.h"
#include "FLA_Trsm_lun.h"
#include "FLA_Trsm_lut.h"
#include "FLA_Trsm_rlc.h"
#include "FLA_Trsm_rlh.h"
#include "FLA_Trsm_rln.h"
#include "FLA_Trsm_rlt.h"
#include "FLA_Trsm_ruc.h"
#include "FLA_Trsm_ruh.h"
#include "FLA_Trsm_run.h"
#include "FLA_Trsm_rut.h"

FLA_Error FLA_Trsm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );

FLA_Error FLA_Trsm_llc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_llh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_llt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_luc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_luh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lun( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_lut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rln( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rlt( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_ruc( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_ruh( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_run( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );
FLA_Error FLA_Trsm_rut( FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl );

