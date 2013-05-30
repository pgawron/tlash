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

// --- FLAME to BLAS/LAPACK mappings -------------------------------------------

void FLA_Param_map_flame_to_netlib_trans( FLA_Trans trans, void* blas_trans )
{
	if ( trans == FLA_NO_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasNoTrans;
#else
		*( ( char*                 ) blas_trans ) = 'N';
#endif
	}
	else if ( trans == FLA_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasTrans;
#else
		*( ( char*                 ) blas_trans ) = 'T';
#endif
	}
	else if ( trans == FLA_CONJ_TRANSPOSE )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasConjTrans;
#else
		*( ( char*                 ) blas_trans ) = 'C';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_TRANS );
	}
}

void FLA_Param_map_flame_to_netlib_uplo( FLA_Uplo uplo, void* blas_uplo )
{
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasLower;
#else
		*( ( char*            ) blas_uplo ) = 'L';
#endif
	}
	else if ( uplo == FLA_UPPER_TRIANGULAR )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasUpper;
#else
		*( ( char*            ) blas_uplo ) = 'U';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_UPLO );
	}
}

void FLA_Param_map_flame_to_netlib_side( FLA_Side side, void* blas_side )
{
	if ( side == FLA_LEFT )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasLeft;
#else
		*( ( char*            ) blas_side ) = 'L';
#endif
	}
	else if ( side == FLA_RIGHT )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasRight;
#else
		*( ( char*            ) blas_side ) = 'R';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SIDE );
	}
}

void FLA_Param_map_flame_to_netlib_diag( FLA_Diag diag, void* blas_diag )
{
	if ( diag == FLA_NONUNIT_DIAG )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasNonUnit;
#else
		*( ( char*            ) blas_diag ) = 'N';
#endif
	}
	else if ( diag == FLA_UNIT_DIAG )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasUnit;
#else
		*( ( char*            ) blas_diag ) = 'U';
#endif
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIAG );
	}
}

void FLA_Param_map_flame_to_netlib_direct( FLA_Direct direct, void* lapack_direct )
{
	if ( direct == FLA_FORWARD )
	{
		*( ( char* ) lapack_direct ) = 'F';
	}
	else if ( direct == FLA_BACKWARD )
	{
		*( ( char* ) lapack_direct ) = 'B';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIRECT );
	}
}

void FLA_Param_map_flame_to_netlib_storev( FLA_Store storev, void* lapack_storev )
{
	if ( storev == FLA_COLUMNWISE )
	{
		*( ( char* ) lapack_storev ) = 'C';
	}
	else if ( storev == FLA_ROWWISE )
	{
		*( ( char* ) lapack_storev ) = 'R';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_STOREV );
	}
}

void FLA_Param_map_flame_to_netlib_evd_type( FLA_Evd_type evd_type, void* lapack_evd_type )
{
	if ( evd_type == FLA_EVD_WITHOUT_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'N';
	}
	else if ( evd_type == FLA_EVD_WITH_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'V';
	}
	else if ( evd_type == FLA_EVD_OF_TRIDIAG_WITH_VECTORS )
	{
		*( ( char* ) lapack_evd_type ) = 'I';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_EVD_TYPE );
	}
}

void FLA_Param_map_flame_to_netlib_svd_type( FLA_Svd_type svd_type, void* lapack_svd_type )
{
	if      ( svd_type == FLA_SVD_VECTORS_ALL )
	{
		*( ( char* ) lapack_svd_type ) = 'A';
	}
	else if ( svd_type == FLA_SVD_VECTORS_MIN_COPY )
	{
		*( ( char* ) lapack_svd_type ) = 'S';
	}
	else if ( svd_type == FLA_SVD_VECTORS_MIN_OVERWRITE )
	{
		*( ( char* ) lapack_svd_type ) = 'O';
	}
	else if ( svd_type == FLA_SVD_VECTORS_NONE )
	{
		*( ( char* ) lapack_svd_type ) = 'N';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SVD_TYPE );
	}
}

void FLA_Param_map_flame_to_netlib_machval( FLA_Machval machval, void* blas_machval )
{
	if      ( machval == FLA_MACH_EPS )
	{
		*( ( char* ) blas_machval ) = 'E';
	}
	else if ( machval == FLA_MACH_SFMIN )
	{
		*( ( char* ) blas_machval ) = 'S';
	}
	else if ( machval == FLA_MACH_BASE )
	{
		*( ( char* ) blas_machval ) = 'B';
	}
	else if ( machval == FLA_MACH_PREC )
	{
		*( ( char* ) blas_machval ) = 'P';
	}
	else if ( machval == FLA_MACH_NDIGMANT )
	{
		*( ( char* ) blas_machval ) = 'N';
	}
	else if ( machval == FLA_MACH_RND )
	{
		*( ( char* ) blas_machval ) = 'R';
	}
	else if ( machval == FLA_MACH_EMIN )
	{
		*( ( char* ) blas_machval ) = 'M';
	}
	else if ( machval == FLA_MACH_RMIN )
	{
		*( ( char* ) blas_machval ) = 'U';
	}
	else if ( machval == FLA_MACH_EMAX )
	{
		*( ( char* ) blas_machval ) = 'L';
	}
	else if ( machval == FLA_MACH_RMAX )
	{
		*( ( char* ) blas_machval ) = 'O';
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_MACHVAL );
	}
}

// --- FLAME to BLIS mappings --------------------------------------------------

void FLA_Param_map_flame_to_blis_trans( FLA_Trans trans, trans_t* blis_trans )
{
	if ( trans == FLA_NO_TRANSPOSE )
	{
		*blis_trans = BLIS_NO_TRANSPOSE;
	}
	else if ( trans == FLA_TRANSPOSE )
	{
		*blis_trans = BLIS_TRANSPOSE;
	}
	else if ( trans == FLA_CONJ_NO_TRANSPOSE )
	{
		*blis_trans = BLIS_CONJ_NO_TRANSPOSE;
	}
	else if ( trans == FLA_CONJ_TRANSPOSE )
	{
		*blis_trans = BLIS_CONJ_TRANSPOSE;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_TRANS );
	}
}

void FLA_Param_map_flame_to_blis_conj( FLA_Conj conj, conj_t* blis_conj )
{
	if ( conj == FLA_NO_CONJUGATE )
	{
		*blis_conj = BLIS_NO_CONJUGATE;
	}
	else if ( conj == FLA_CONJUGATE )
	{
		*blis_conj = BLIS_CONJUGATE;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_CONJ );
	}
}

void FLA_Param_map_flame_to_blis_uplo( FLA_Uplo uplo, uplo_t* blis_uplo )
{
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		*blis_uplo = BLIS_LOWER_TRIANGULAR;
	}
	else if ( uplo == FLA_UPPER_TRIANGULAR )
	{
		*blis_uplo = BLIS_UPPER_TRIANGULAR;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_UPLO );
	}
}

void FLA_Param_map_flame_to_blis_side( FLA_Side side, side_t* blis_side )
{
	if ( side == FLA_LEFT )
	{
		*blis_side = BLIS_LEFT;
	}
	else if ( side == FLA_RIGHT )
	{
		*blis_side = BLIS_RIGHT;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_SIDE );
	}
}

void FLA_Param_map_flame_to_blis_diag( FLA_Diag diag, diag_t* blis_diag )
{
	if ( diag == FLA_NONUNIT_DIAG )
	{
		*blis_diag = BLIS_NONUNIT_DIAG;
	}
	else if ( diag == FLA_UNIT_DIAG )
	{
		*blis_diag = BLIS_UNIT_DIAG;
	}
	else
	{
		FLA_Check_error_code( FLA_INVALID_DIAG );
	}
}

// --- BLAS/LAPACK to FLAME mappings -------------------------------------------

void FLA_Param_map_netlib_to_flame_trans( char* trans, FLA_Trans* flame_trans )
{
	if      ( *trans == 'n' || *trans == 'N' )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( *trans == 't' || *trans == 'T' )
		*flame_trans = FLA_TRANSPOSE;
	else if ( *trans == 'c' || *trans == 'C' )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_netlib_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo )
{
	if      ( *uplo == 'l' || *uplo == 'L' )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( *uplo == 'u' || *uplo == 'U' )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_netlib_to_flame_side( char* side, FLA_Side* flame_side )
{
	if      ( *side == 'l' || *side == 'L' )
		*flame_side = FLA_LEFT;
	else if ( *side == 'r' || *side == 'R' )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_netlib_to_flame_diag( char* diag, FLA_Diag* flame_diag )
{
	if      ( *diag == 'n' || *diag == 'N' )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( *diag == 'u' || *diag == 'U' )
		*flame_diag = FLA_UNIT_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

void FLA_Param_map_netlib_to_flame_inv( int* itype, FLA_Inv* flame_inv )
{
	if      ( *itype == 1 )
		*flame_inv = FLA_INVERSE;
	else if ( *itype == 2 || *itype == 3 )
		*flame_inv = FLA_NO_INVERSE;
	else
		FLA_Check_error_code( FLA_INVALID_INVERSE );
}

// --- BLIS to FLAME mappings --------------------------------------------------

void FLA_Param_map_blis_to_flame_trans( trans_t trans, FLA_Trans* flame_trans )
{
	if      ( bli_is_notrans( trans ) )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( bli_is_trans( trans ) )
		*flame_trans = FLA_TRANSPOSE;
	else if ( bli_is_conjnotrans( trans ) )
		*flame_trans = FLA_CONJ_NO_TRANSPOSE;
	else if ( bli_is_conjtrans( trans ) )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_blis_to_flame_uplo( uplo_t uplo, FLA_Uplo* flame_uplo )
{
	if      ( bli_is_lower( uplo ) )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( bli_is_upper( uplo ) )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_blis_to_flame_side( side_t side, FLA_Side* flame_side )
{
	if      ( bli_is_left( side ) )
		*flame_side = FLA_LEFT;
	else if ( bli_is_right( side ) )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_blis_to_flame_diag( diag_t diag, FLA_Diag* flame_diag )
{
	if      ( bli_is_nonunit_diag( diag ) )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( bli_is_unit_diag( diag ) )
		*flame_diag = FLA_UNIT_DIAG;
	else if ( bli_is_zero_diag( diag ) )
		*flame_diag = FLA_ZERO_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

// --- FLAME char to FLAME mappings --------------------------------------------

void FLA_Param_map_char_to_flame_trans( char* trans, FLA_Trans* flame_trans )
{
	if      ( *trans == 'n' || *trans == 'N' )
		*flame_trans = FLA_NO_TRANSPOSE;
	else if ( *trans == 't' || *trans == 'T' )
		*flame_trans = FLA_TRANSPOSE;
	else if ( *trans == 'c' || *trans == 'C' )
		*flame_trans = FLA_CONJ_NO_TRANSPOSE;
	else if ( *trans == 'h' || *trans == 'H' )
		*flame_trans = FLA_CONJ_TRANSPOSE;
	else
		FLA_Check_error_code( FLA_INVALID_TRANS );
}

void FLA_Param_map_char_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo )
{
	if      ( *uplo == 'l' || *uplo == 'L' )
		*flame_uplo = FLA_LOWER_TRIANGULAR;
	else if ( *uplo == 'u' || *uplo == 'U' )
		*flame_uplo = FLA_UPPER_TRIANGULAR;
	else
		FLA_Check_error_code( FLA_INVALID_UPLO );
}

void FLA_Param_map_char_to_flame_side( char* side, FLA_Side* flame_side )
{
	if      ( *side == 'l' || *side == 'L' )
		*flame_side = FLA_LEFT;
	else if ( *side == 'r' || *side == 'R' )
		*flame_side = FLA_RIGHT;
	else
		FLA_Check_error_code( FLA_INVALID_SIDE );
}

void FLA_Param_map_char_to_flame_diag( char* diag, FLA_Diag* flame_diag )
{
	if      ( *diag == 'n' || *diag == 'N' )
		*flame_diag = FLA_NONUNIT_DIAG;
	else if ( *diag == 'u' || *diag == 'U' )
		*flame_diag = FLA_UNIT_DIAG;
	else
		FLA_Check_error_code( FLA_INVALID_DIAG );
}

void FLA_Param_map_char_to_flame_direct( char* direct, FLA_Direct* flame_direct )
{
	if      ( *direct == 'b' || *direct == 'B' )
		*flame_direct = FLA_BACKWARD;
	else if ( *direct == 'f' || *direct == 'F' )
		*flame_direct = FLA_FORWARD;
	else
		FLA_Check_error_code( FLA_INVALID_DIRECT );
}

void FLA_Param_map_char_to_flame_storev( char* storev, FLA_Direct* flame_storev )
{
	if      ( *storev == 'c' || *storev == 'C' )
		*flame_storev = FLA_COLUMNWISE;
	else if ( *storev == 'r' || *storev == 'R' )
		*flame_storev = FLA_ROWWISE;
	else
		FLA_Check_error_code( FLA_INVALID_STOREV );
}

void FLA_Param_map_char_to_flame_inv( char* inv, FLA_Inv* flame_inv )
{
	if      ( *inv == 'i' || *inv == 'I' )
		*flame_inv = FLA_INVERSE;
	else if ( *inv == 'n' || *inv == 'N' )
		*flame_inv = FLA_NO_INVERSE;
	else
		FLA_Check_error_code( FLA_INVALID_INVERSE );
}

