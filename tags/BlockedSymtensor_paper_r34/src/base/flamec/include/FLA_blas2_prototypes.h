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

FLA_Error FLA_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );


// --- task wrapper prototypes -------------------------------------------------

FLA_Error FLA_Gemv_task( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Trsv_task( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

FLA_Error FLA_Gemv_h_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_n_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Gemv_t_task( FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );

FLA_Error FLA_Trsv_lc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ln_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_lt_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_uc_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_un_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );
FLA_Error FLA_Trsv_ut_task( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );


// --- external wrapper prototypes ---------------------------------------------

FLA_Error FLA_Gemv_external( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc_external( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger_external( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc_external( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc_external( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc_external( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c_external( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );


// --- gpu wrapper prototypes --------------------------------------------------

FLA_Error FLA_Gemv_external_gpu( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu, FLA_Obj beta, FLA_Obj y, void* y_gpu );
FLA_Error FLA_Trsv_external_gpu( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, void* A_gpu, FLA_Obj x, void* x_gpu );


// --- check routine prototypes ------------------------------------------------

// front-ends
FLA_Error FLA_Gemv_check( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Gemvc_check( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Ger_check( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Gerc_check( FLA_Conj conjx, FLA_Conj conjy, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Hemv_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Hemvc_check( FLA_Uplo uplo, FLA_Conj conja, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Her_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Herc_check( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Her2_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Her2c_check( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Symv_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Syr_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Syr2_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A );
FLA_Error FLA_Trmv_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trmvsx_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );
FLA_Error FLA_Trsv_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x );
FLA_Error FLA_Trsvsx_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y );

// internal back-ends
FLA_Error FLA_Gemv_internal_check( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl );
FLA_Error FLA_Trsv_internal_check( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl );

