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

//--------------------------------------------------------------------------
// Required supermatrix definitions
//--------------------------------------------------------------------------

void FLASH_Queue_begin( void );
void FLASH_Queue_end( void );
FLA_Error FLASH_Queue_enable( void );
FLA_Error FLASH_Queue_disable( void );
FLA_Bool FLASH_Queue_get_enabled( void );

//--------------------------------------------------------------------------
// Required supermatrix macros
//--------------------------------------------------------------------------

// Level-3 BLAS

#define ENQUEUE_FLASH_Gemm( transA, transB, alpha, A, B, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Hemm( side, uplo, alpha, A, B, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Herk( uplo, transA, alpha, A, beta, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Her2k( uplo, transA, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Symm( side, uplo, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Syrk( uplo, transA, alpha, A, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Syr2k( uplo, transA, alpha, A, B, beta, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trmm( side, uplo, trans, diag, alpha, A, C, cntl  ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trsm( side, uplo, trans, diag, alpha, A, C, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Level-2 BLAS

#define ENQUEUE_FLASH_Gemv( transA, alpha, A, x, beta, y, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Trsv( uplo, trans, diag, A, x, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Level-1 BLAS

#define ENQUEUE_FLASH_Axpy( alpha, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Axpyt( trans, alpha, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copy( A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copyt( trans, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Copyr( uplo, A, B, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Scal( alpha, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Scalr( uplo, alpha, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

// Base

#define ENQUEUE_FLASH_Obj_create_buffer( rs, cs, A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )

#define ENQUEUE_FLASH_Obj_free_buffer( A, cntl ) \
        FLA_Check_error_code( FLA_SUPERMATRIX_NOT_ENABLED )
