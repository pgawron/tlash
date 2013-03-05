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

#include "blis.h"

// --- BLIS to BLAS/LAPACK mappings --------------------------------------------

void bli_param_map_to_netlib_trans( trans_t blis_trans, void* blas_trans )
{
	if ( bli_is_notrans( blis_trans ) || bli_is_conjnotrans( blis_trans ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasNoTrans;
#else
		*( ( char*                 ) blas_trans ) = 'N';
#endif
	}
	else if ( bli_is_trans( blis_trans ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasTrans;
#else
		*( ( char*                 ) blas_trans ) = 'T';
#endif
	}
	else if ( bli_is_conjtrans( blis_trans ))
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_TRANSPOSE* ) blas_trans ) = CblasConjTrans;
#else
		*( ( char*                 ) blas_trans ) = 'C';
#endif
	}
	else
	{
		bli_abort_msg( "Invalid BLIS trans value to map." );
	}
}

void bli_param_map_to_netlib_uplo( uplo_t blis_uplo, void* blas_uplo )
{
	if ( bli_is_lower( blis_uplo ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasLower;
#else
		*( ( char*            ) blas_uplo ) = 'L';
#endif
	}
	else if ( bli_is_upper( blis_uplo ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_UPLO* ) blas_uplo ) = CblasUpper;
#else
		*( ( char*            ) blas_uplo ) = 'U';
#endif
	}
	else
	{
		bli_abort_msg( "Invalid BLIS uplo value to map." );
	}
}

void bli_param_map_to_netlib_side( side_t blis_side, void* blas_side )
{
	if ( bli_is_left( blis_side ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasLeft;
#else
		*( ( char*            ) blas_side ) = 'L';
#endif
	}
	else if ( bli_is_right( blis_side ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_SIDE* ) blas_side ) = CblasRight;
#else
		*( ( char*            ) blas_side ) = 'R';
#endif
	}
	else
	{
		bli_abort_msg( "Invalid BLIS side value to map." );
	}
}

void bli_param_map_to_netlib_diag( diag_t blis_diag, void* blas_diag )
{
	if ( bli_is_nonunit_diag( blis_diag ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasNonUnit;
#else
		*( ( char*            ) blas_diag ) = 'N';
#endif
	}
	else if ( bli_is_unit_diag( blis_diag ) )
	{
#ifdef FLA_ENABLE_CBLAS_INTERFACES
		*( ( enum CBLAS_DIAG* ) blas_diag ) = CblasUnit;
#else
		*( ( char*            ) blas_diag ) = 'U';
#endif
	}
	else
	{
		bli_abort_msg( "Invalid BLIS diag value to map." );
	}
}

