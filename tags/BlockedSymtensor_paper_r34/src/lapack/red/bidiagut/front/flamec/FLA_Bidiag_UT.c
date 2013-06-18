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

extern fla_bidiagut_t* fla_bidiagut_cntl_fused;
extern fla_bidiagut_t* fla_bidiagut_cntl_nofus;

FLA_Error FLA_Bidiag_UT( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_check( A, TU, TV );

  if ( FLA_Obj_row_stride( A ) == 1 &&
       FLA_Obj_row_stride( TU ) == 1 &&
       FLA_Obj_row_stride( TV ) == 1 &&
       FLA_Obj_is_double_precision( A ) )
    r_val = FLA_Bidiag_UT_internal( A, TU, TV, fla_bidiagut_cntl_fused );
  else
    r_val = FLA_Bidiag_UT_internal( A, TU, TV, fla_bidiagut_cntl_nofus );

  return r_val;
}
