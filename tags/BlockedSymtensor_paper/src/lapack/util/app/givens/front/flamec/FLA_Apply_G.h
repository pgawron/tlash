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

#include "FLA_Apply_G_lf.h"
#include "FLA_Apply_G_lb.h"
#include "FLA_Apply_G_rf.h"
#include "FLA_Apply_G_rb.h"

FLA_Error FLA_Apply_G( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_internal( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A );

#include "FLA_Givens2.h"

#include "FLA_Apply_GTG.h"

#include "FLA_Apply_GT_2x2.h"
#include "FLA_Apply_G_2x2.h"
#include "FLA_Apply_G_1x2.h"

#include "FLA_Apply_G_mx2_opt.h"
#include "FLA_Apply_G_mx2_asm.h"

#include "FLA_Apply_G_mx3_opt.h"
#include "FLA_Apply_G_mx3b_opt.h"
#include "FLA_Apply_G_mx3_asm.h"
#include "FLA_Apply_G_mx3b_asm.h"

#include "FLA_Apply_G_mx4s_opt.h"
#include "FLA_Apply_G_mx4s_asm.h"

