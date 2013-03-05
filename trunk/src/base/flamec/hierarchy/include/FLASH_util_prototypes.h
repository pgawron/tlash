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

// --- FLASH utility routine prototypes ----------------------------------------

double    FLASH_Max_elemwise_diff( FLA_Obj A, FLA_Obj B );

FLA_Error FLASH_Random_matrix( FLA_Obj H );
FLA_Error FLASH_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj H );

FLA_Error FLASH_Norm1( FLA_Obj H, FLA_Obj norm );
FLA_Error FLASH_Obj_shift_diagonal( FLA_Conj conj, FLA_Obj sigma, FLA_Obj H );

FLA_Error FLASH_Set( FLA_Obj alpha, FLA_Obj H );

FLA_Error FLASH_Obj_create_diag_panel( FLA_Obj A, FLA_Obj* U );

FLA_Error FLASH_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLASH_Hermitianize( FLA_Uplo uplo, FLA_Obj A );


