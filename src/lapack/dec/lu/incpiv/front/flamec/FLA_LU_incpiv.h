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

#include "FLA_LU_incpiv_aux.h"

FLA_Error FLASH_LU_incpiv_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* p, FLA_Obj* L );
dim_t     FLASH_LU_incpiv_determine_alg_blocksize( FLA_Obj A );

FLA_Error FLASH_LU_incpiv_noopt( FLA_Obj A, FLA_Obj p, FLA_Obj L );
FLA_Error FLASH_LU_incpiv_opt1( FLA_Obj A, FLA_Obj p, FLA_Obj L );

FLA_Error FLASH_LU_incpiv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj B, FLA_Obj X );