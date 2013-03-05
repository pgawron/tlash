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

extern fla_appiv_t* fla_appiv_cntl_leaf;

FLA_Error FLA_Apply_pivots_task( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   //return FLA_Apply_pivots_unb_external( side, trans, p, A );
   return FLA_Apply_pivots_internal( side, trans,
                                     p, A,
                                     fla_appiv_cntl_leaf );
}

FLA_Error FLA_Apply_pivots_ln_task( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   //return FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );
   return FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE,
                                     p, A,
                                     fla_appiv_cntl_leaf );
}

