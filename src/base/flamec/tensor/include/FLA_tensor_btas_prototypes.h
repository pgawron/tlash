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

// --- top-level wrapper prototypes --------------------------------------------

FLA_Error FLA_Permute_single( FLA_Obj A, dim_t permtuation[], FLA_Obj* B);
FLA_Error FLA_Permute_hier( FLA_Obj A, dim_t permtuation[], FLA_Obj* B);
FLA_Error FLA_Permute_single_inplace( FLA_Obj* A, dim_t permutation[]);
FLA_Error FLA_Permute( dim_t permutation[], dim_t order, dim_t size[], FLA_Obj A[], dim_t* B_order, dim_t (*B_size)[], FLA_Obj (*B)[] );
FLA_Error FLA_Ttm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm_single_no_permC( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm_single_mode( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm_single_mode_no_permC( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C );
FLA_Error FLA_Ttm( FLA_Obj alpha, FLA_Obj A, dim_t nModes, dim_t mode[nModes], FLA_Obj beta, FLA_Obj B[nModes], FLA_Obj C );
FLA_Error FLA_Ttm_hierC_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C);
FLA_Error FLA_Ttm_hierA_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C);
FLA_Error FLA_Sttsm_single( FLA_Obj alpha, FLA_Obj A, dim_t mode, FLA_Obj beta, FLA_Obj B, FLA_Obj C, dim_t maxIndex );
FLA_Error FLA_Sttsm( FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj B, FLA_Obj C );

// --- check routine prototypes ------------------------------------------------

