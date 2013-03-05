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

// --- Query routine prototypes ------------------------------------------------

// --- trans ---

int bli_does_trans( trans_t trans );
int bli_does_notrans( trans_t trans );
int bli_does_conj( trans_t trans );

int bli_is_notrans( trans_t trans );
int bli_is_trans( trans_t trans );
int bli_is_conjnotrans( trans_t trans );
int bli_is_conjtrans( trans_t trans );

// --- conj ---

int bli_is_noconj( conj_t conj );
int bli_is_conj( conj_t conj );

// --- uplo ---

int bli_is_lower( uplo_t uplo );
int bli_is_upper( uplo_t uplo );

// --- side ---

int bli_is_left( side_t side );
int bli_is_right( side_t side );

// --- diag ---

int bli_is_nonunit_diag( diag_t diag );
int bli_is_unit_diag( diag_t diag );
int bli_is_zero_diag( diag_t diag );

// --- mapping-related ---

conj_t bli_proj_trans_to_conj( trans_t trans );

// --- storage-related ---

void bli_check_storage_3m( int a_rs, int a_cs, int b_rs, int b_cs, int c_rs, int c_cs );
void bli_check_storage_2m( int a_rs, int a_cs, int b_rs, int b_cs );
int bli_is_row_or_col_storage( int rs, int cs );
int bli_is_row_storage( int rs, int cs );
int bli_is_col_storage( int rs, int cs );
int bli_is_gen_storage( int rs, int cs );
int bli_is_vector( int m, int n );

// --- vector-related ---

int bli_vector_dim( int m, int n );
int bli_vector_inc( trans_t trans, int m, int n, int rs, int cs );

// --- dimension-related ---

int bli_zero_dim1( int m );
int bli_zero_dim2( int m, int n );
int bli_zero_dim3( int m, int k, int n );

