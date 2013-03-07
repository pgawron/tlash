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

extern fla_ttmm_t* fla_ttmm_cntl_leaf;

FLA_Error FLA_Ttmm_task( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl )
{
  return FLA_Ttmm_internal( uplo, A,
                            fla_ttmm_cntl_leaf );
}

FLA_Error FLA_Ttmm_l_task( FLA_Obj A, fla_ttmm_t* cntl )
{
  //return FLA_Ttmm_unb_external( FLA_LOWER_TRIANGULAR, A );
  return FLA_Ttmm_internal( FLA_LOWER_TRIANGULAR, A,
                            fla_ttmm_cntl_leaf );
}

FLA_Error FLA_Ttmm_u_task( FLA_Obj A, fla_ttmm_t* cntl )
{
  //return FLA_Ttmm_unb_external( FLA_UPPER_TRIANGULAR, A );
  return FLA_Ttmm_internal( FLA_UPPER_TRIANGULAR, A,
                            fla_ttmm_cntl_leaf );
}
