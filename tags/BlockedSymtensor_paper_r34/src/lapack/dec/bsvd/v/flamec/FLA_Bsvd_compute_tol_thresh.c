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

FLA_Error FLA_Bsvd_compute_tol_thresh( FLA_Obj tolmul, FLA_Obj maxitr, FLA_Obj d, FLA_Obj e, FLA_Obj tol, FLA_Obj thresh )
{
	FLA_Datatype datatype;
	int          n_A;
	int          inc_d;
	int          inc_e;

	datatype = FLA_Obj_datatype( d );

	n_A      = FLA_Obj_vector_dim( d );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_tolmul = FLA_FLOAT_PTR( tolmul );
			float*    buff_maxitr = FLA_FLOAT_PTR( maxitr );
			float*    buff_d      = FLA_FLOAT_PTR( d );
			float*    buff_e      = FLA_FLOAT_PTR( e );
			float*    buff_tol    = FLA_FLOAT_PTR( tol );
			float*    buff_thresh = FLA_FLOAT_PTR( thresh );

			FLA_Bsvd_compute_tol_thresh_ops( n_A,
			                                 *buff_tolmul,
			                                 *buff_maxitr,
			                                 buff_d, inc_d,
			                                 buff_e, inc_e,
			                                 buff_tol,
			                                 buff_thresh );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_tolmul = FLA_DOUBLE_PTR( tolmul );
			double*   buff_maxitr = FLA_DOUBLE_PTR( maxitr );
			double*   buff_d      = FLA_DOUBLE_PTR( d );
			double*   buff_e      = FLA_DOUBLE_PTR( e );
			double*   buff_tol    = FLA_DOUBLE_PTR( tol );
			double*   buff_thresh = FLA_DOUBLE_PTR( thresh );

			FLA_Bsvd_compute_tol_thresh_opd( n_A,
			                                 *buff_tolmul,
			                                 *buff_maxitr,
			                                 buff_d, inc_d,
			                                 buff_e, inc_e,
			                                 buff_tol,
			                                 buff_thresh );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_tol_thresh_ops( int       n_A,
                                           float     tolmul,
                                           float     maxitr,
                                           float*    buff_d, int inc_d, 
                                           float*    buff_e, int inc_e, 
                                           float*    tol,
                                           float*    thresh )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_tol_thresh_opd( int       n_A,
                                           double    tolmul,
                                           double    maxitr,
                                           double*   buff_d, int inc_d, 
                                           double*   buff_e, int inc_e, 
                                           double*   tol,
                                           double*   thresh )
{
	double zero = bli_d0();
	double smin;
	double eps, unfl;
	double mu;
	int    i;

	// Query machine epsilon and the safe minimum.
	eps  = FLA_Mach_params_opd( FLA_MACH_EPS );
	unfl = FLA_Mach_params_opd( FLA_MACH_SFMIN );

	// Compute tol as the product of machine epsilon and tolmul.
	*tol = tolmul * eps;

	// Compute the approximate maximum singular value.
	// NOT NEEDED unless we're supporting absolute accuracy.
	//FLA_Bsvd_sinval_find_max( n_A,
	//                          buff_d, inc_d,
	//                          buff_e, inc_e,
	//                          &smax );

	// Compute the approximate minimum singular value.
	smin = fabs( *buff_d );

	// Skip the accumulation of smin if the first element is zero.
	if ( smin != zero )
	{
		mu = smin;
		for ( i = 1; i < n_A; ++i )
		{
			double* epsilon1 = buff_e + (i-1)*inc_e;
			double* delta2   = buff_d + (i  )*inc_d;

			mu   = fabs( *delta2 ) * ( mu / ( mu + fabs( *epsilon1 ) ) );
			smin = min( smin, mu );

			// Stop early if we encountered a zero.
			if ( smin == zero ) break;
		}
	}

	// Compute thresh either in terms of tol or as a function of the
	// maximum total number of iterations, the problem size, and the
	// safe minimum.
	smin = smin / sqrt( ( double ) n_A );
	*thresh = max( *tol * smin, maxitr * n_A * n_A * unfl );

	return FLA_SUCCESS;
}

