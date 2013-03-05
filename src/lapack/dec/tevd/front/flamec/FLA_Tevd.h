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

#include "FLA_Tevd_n.h"
#include "FLA_Tevd_v.h"

// --- MAC_Tevd_eigval_converged() ---------------------------------------------

#define MAC_Tevd_eigval_converged_ops( eps, safmin, d1, e1, d2 ) \
	fabsf( e1 ) <= (eps) * sqrt( fabsf( d1 ) ) * sqrt( fabsf( d2 ) ) + (safmin)

#define MAC_Tevd_eigval_converged_opd( eps, safmin, d1, e1, d2 ) \
	fabs( e1 )  <= (eps) * sqrt( fabs( d1 ) )  * sqrt( fabs( d2 ) )  + (safmin)

// --- MAC_Tevd_eigval_converged2() ---------------------------------------------

#define MAC_Tevd_eigval_converged2_ops( eps2, safmin, d1, e1, d2 ) \
	(e1) * (e1) <=        (eps2) * fabsf( (d1) * (d2) ) + (safmin)

#define MAC_Tevd_eigval_converged2_opd( eps2, safmin, d1, e1, d2 ) \
	(e1) * (e1) <=        (eps2) * fabs( (d1) * (d2) ) + (safmin)

FLA_Error FLA_Tevd( FLA_Evd_type jobz, FLA_Obj U, FLA_Obj d, FLA_Obj e, FLA_Obj l );

