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

// Factorizations
#include "FLA_Chol.h"
#include "FLA_LU_nopiv.h"
#include "FLA_LU_piv.h"
#include "FLA_LU_incpiv.h"
#include "FLA_QR_UT.h"
#include "FLA_QR2_UT.h"
#include "FLA_QR_UT_inc.h"
#include "FLA_LQ_UT.h"
#include "FLA_CAQR2_UT.h"
#include "FLA_CAQR_UT_inc.h"

// Other Decompositions
#include "FLA_Hevd.h"
#include "FLA_Tevd.h"
#include "FLA_Svd.h"
#include "FLA_Bsvd.h"

// Inversions
#include "FLA_Trinv.h"
#include "FLA_SPDinv.h"

// Reductions
#include "FLA_Hess_UT.h"
#include "FLA_Tridiag_UT.h"
#include "FLA_Bidiag_UT.h"

// Solves
#include "FLA_Lyap.h"
#include "FLA_Sylv.h"

// Miscellaneous
#include "FLA_Ttmm.h"
#include "FLA_UDdate_UT.h"
#include "FLA_UDdate_UT_inc.h"

// Utility
#include "FLA_Accum_T_UT.h"
#include "FLA_Apply_G.h"
#include "FLA_Apply_H2_UT.h"
#include "FLA_Apply_HUD_UT.h"
#include "FLA_Apply_Q_UT.h"
#include "FLA_Apply_Q2_UT.h"
#include "FLA_Apply_CAQ2_UT.h"
#include "FLA_Apply_QUD_UT.h"
#include "FLA_Apply_Q_UT_inc.h"
#include "FLA_Apply_CAQ_UT_inc.h"
#include "FLA_Apply_QUD_UT_inc.h"
#include "FLA_Apply_pivots.h"

// Eigensolvers
#include "FLA_Eig_gest.h"
