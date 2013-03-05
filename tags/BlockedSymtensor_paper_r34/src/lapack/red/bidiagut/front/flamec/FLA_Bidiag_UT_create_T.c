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

FLA_Error FLA_Bidiag_UT_create_T( FLA_Obj A, FLA_Obj* TU, FLA_Obj* TV )
{
  FLA_Datatype datatype;
  dim_t        b_alg, k;
  dim_t        rs_T, cs_T;

  // Query the datatype of A.
  datatype = FLA_Obj_datatype( A );

  // Query the blocksize from the library.
  b_alg = FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  // Scale the blocksize by a pre-set global constant.
  b_alg = ( dim_t )( ( ( double ) b_alg ) * FLA_BIDIAG_INNER_TO_OUTER_B_RATIO );

  // Query the minimum dimension of A.
  k = FLA_Obj_min_dim( A );

  // Figure out whether TU and TV should be row-major or column-major.
  if ( FLA_Obj_row_stride( A ) == 1 )
  {
    rs_T = 1;
    cs_T = b_alg;
  }
  else // if ( FLA_Obj_col_stride( A ) == 1 )
  {
    rs_T = k;
    cs_T = 1;
  }

  // Create two b_alg x k matrices to hold the block Householder transforms
  // that will be accumulated within the bidiagonal reduction algorithm.
  FLA_Obj_create( datatype, b_alg, k, rs_T, cs_T, TU );
  FLA_Obj_create( datatype, b_alg, k, rs_T, cs_T, TV );

  return FLA_SUCCESS;
}

