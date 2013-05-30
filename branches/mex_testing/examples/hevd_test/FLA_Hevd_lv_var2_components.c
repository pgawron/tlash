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

FLA_Error FLA_Hevd_lv_var2_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Uplo     uplo = FLA_LOWER_TRIANGULAR;
	FLA_Datatype dt;
	FLA_Datatype dt_real;
	FLA_Datatype dt_comp;
	FLA_Obj      T, r, d, e, G, R, W;
	dim_t        mn_A;
	dim_t        n_G = k_accum;
	double       dtime_temp;

	mn_A    = FLA_Obj_length( A );
	dt      = FLA_Obj_datatype( A );
	dt_real = FLA_Obj_datatype_proj_to_real( A );
	dt_comp = FLA_Obj_datatype_proj_to_complex( A );

	// If the matrix is a scalar, then the EVD is easy.
	if ( mn_A == 1 )
	{
		FLA_Copy( A, l );
		FLA_Set( FLA_ONE, A );

		return FLA_SUCCESS;
	}

	// Create a matrix to hold block Householder transformations.
	FLA_Tridiag_UT_create_T( A, &T );

	// Create a vector to hold the realifying scalars.
	FLA_Obj_create( dt, mn_A, 1, 0, 0, &r );

	// Create vectors to hold the diagonal and sub-diagonal.
	FLA_Obj_create( dt_real, mn_A,     1, 0, 0, &d );
	FLA_Obj_create( dt_real, mn_A-1,   1, 0, 0, &e );
	FLA_Obj_create( dt_comp, mn_A-1, n_G, 0, 0, &G );
	FLA_Obj_create( dt_real, mn_A,  mn_A, 0, 0, &R );
	FLA_Obj_create( dt,      mn_A,  mn_A, 0, 0, &W );

  dtime_temp = FLA_Clock();
  {
	// Reduce the matrix to tridiagonal form.
	FLA_Tridiag_UT( uplo, A, T );
  }
  *dtime_tred = FLA_Clock() - dtime_temp;

	// Apply scalars to rotate elements on the sub-diagonal to the real domain.
	FLA_Tridiag_UT_realify( uplo, A, r );

	// Extract the diagonal and sub-diagonal from A.
	FLA_Tridiag_UT_extract_diagonals( uplo, A, d, e );

  dtime_temp = FLA_Clock();
  {
	// Form Q, overwriting A.
	FLA_Tridiag_UT_form_Q( uplo, A, T );
  }
  *dtime_appq = FLA_Clock() - dtime_temp;

	// Apply the scalars in r to Q.
	FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, r, A );

  dtime_temp = FLA_Clock();
  {
	// Perform an eigenvalue decomposition on the tridiagonal matrix.
	r_val = FLA_Tevd_v_opt_var2( n_iter_max, d, e, G, R, W, A, b_alg );
  }
  *dtime_tevd = FLA_Clock() - dtime_temp;

	// Copy the converged eigenvalues to the output vector.
	FLA_Copy( d, l );

	// Sort the eigenvalues and eigenvectors in ascending order.
	FLA_Sort_evd( FLA_FORWARD, l, A );

	FLA_Obj_free( &T );
	FLA_Obj_free( &r );
	FLA_Obj_free( &d );
	FLA_Obj_free( &e );
	FLA_Obj_free( &G );
	FLA_Obj_free( &R );
	FLA_Obj_free( &W );

	return r_val;
}

