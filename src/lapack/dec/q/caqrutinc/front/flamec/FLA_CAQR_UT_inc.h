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

FLA_Error FLASH_CAQR_UT_inc( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW );

FLA_Error FLASH_CAQR_UT_inc_noopt( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW );

FLA_Error FLASH_CAQR_UT_inc_create_hier_matrices( dim_t p, FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* ATW, FLA_Obj* R, FLA_Obj* RTW );
dim_t     FLASH_CAQR_UT_inc_determine_alg_blocksize( FLA_Obj A );
FLA_Error FLASH_CAQR_UT_inc_adjust_views( FLA_Obj A, FLA_Obj TW );

void      FLA_CAQR_UT_inc_init_structure( dim_t p, dim_t nb_part, FLA_Obj R );

dim_t     FLA_CAQR_UT_inc_compute_blocks_per_part( dim_t p, FLA_Obj A );

FLA_Error FLA_CAQR_UT_inc_factorize_panels( dim_t nb_part, FLA_Obj A, FLA_Obj ATW );

FLA_Error FLA_CAQR_UT_inc_copy_triangles( dim_t nb_part, FLA_Obj A, FLA_Obj R );

FLA_Error FLA_CAQR_UT_inc_blk_var1( FLA_Obj R, FLA_Obj TW, fla_caqrutinc_t* cntl );

FLA_Error FLASH_CAQR_UT_inc_solve( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj B, FLA_Obj X );

