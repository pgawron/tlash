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

extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_copy_task( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl )
{
  FLA_Error r_val;
  FLA_Obj   AT,
            AB;

  // Perform a QR factorization as we normally would.
  r_val = FLA_QR_UT_internal( A, T,
                              fla_qrut_cntl_leaf );

  // Partition away the bottom part of the matrix, if there is any, so that
  // the dimensions match that of U. This step is only necessary so that
  // the copyr operation below works normally for the last iteration of
  // incremental QR. The whole point of making a copy of the lower triangle
  // of A is to allow stage 2 to proceed without a dependency leading into
  // stage 3. But the last iteration does not perform stage 2, and thus U
  // is never read and so the copyr does not need to happen.
  FLA_Part_2x1( A,   &AT,
                     &AB,    FLA_Obj_min_dim( A ), FLA_TOP );

  // Copy the Householder vectors into U.
  FLA_Copyr_external( FLA_LOWER_TRIANGULAR, AT, U );

  return r_val;
}

