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

#include "FLA_Apply_Q_UT_lnfc.h"
#include "FLA_Apply_Q_UT_lnfr.h"
#include "FLA_Apply_Q_UT_lnbc.h"
#include "FLA_Apply_Q_UT_lnbr.h"
#include "FLA_Apply_Q_UT_lhfc.h"
#include "FLA_Apply_Q_UT_lhfr.h"
#include "FLA_Apply_Q_UT_lhbc.h"
#include "FLA_Apply_Q_UT_lhbr.h"

#include "FLA_Apply_Q_UT_rhbc.h"
#include "FLA_Apply_Q_UT_rhbr.h"
#include "FLA_Apply_Q_UT_rhfc.h"
#include "FLA_Apply_Q_UT_rhfr.h"
#include "FLA_Apply_Q_UT_rnbc.h"
#include "FLA_Apply_Q_UT_rnbr.h"
#include "FLA_Apply_Q_UT_rnfc.h"
#include "FLA_Apply_Q_UT_rnfr.h"

FLA_Error FLA_Apply_Q_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_lnfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_rhbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfc( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );

FLA_Error FLA_Apply_Q_UT_create_workspace( FLA_Obj T, FLA_Obj B, FLA_Obj* W );

FLA_Error FLASH_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLASH_Apply_Q_UT_create_workspace( FLA_Obj TW, FLA_Obj B, FLA_Obj* W );

