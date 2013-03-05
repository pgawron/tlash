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

extern fla_qr2ut_t* flash_qr2ut_cntl;
extern fla_qr2ut_t* flash_qr2ut_cntl_leaf;
extern fla_qr2ut_t* fla_qr2ut_cntl_leaf;

FLA_Error FLASH_QR2_UT( FLA_Obj B, FLA_Obj D, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR2_UT_check( B, D, T );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Invoke FLA_QR2_UT_internal() with the standard control tree.
  r_val = FLA_QR2_UT_internal( B, D, T, flash_qr2ut_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

