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

FLA_Error FLA_Bsvd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj RG, FLA_Obj RH, FLA_Obj W, FLA_Obj U, FLA_Obj V, dim_t b_alg )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Datatype datatype;
	int          m_U, m_V, n_GH;
	int          inc_d;
	int          inc_e;
	int          rs_G, cs_G;
	int          rs_H, cs_H;
	int          rs_RG, cs_RG;
	int          rs_RH, cs_RH;
	int          rs_W, cs_W;
	int          rs_U, cs_U;
	int          rs_V, cs_V;

	datatype = FLA_Obj_datatype( U );

	m_U        = FLA_Obj_length( U );
	m_V        = FLA_Obj_length( V );
	n_GH       = FLA_Obj_width( G );

	inc_d      = FLA_Obj_vector_inc( d );
	inc_e      = FLA_Obj_vector_inc( e );
	
	rs_G       = FLA_Obj_row_stride( G );
	cs_G       = FLA_Obj_col_stride( G );

	rs_H       = FLA_Obj_row_stride( H );
	cs_H       = FLA_Obj_col_stride( H );

	rs_RG      = FLA_Obj_row_stride( RG );
	cs_RG      = FLA_Obj_col_stride( RG );

	rs_RH      = FLA_Obj_row_stride( RH );
	cs_RH      = FLA_Obj_col_stride( RH );

	rs_W       = FLA_Obj_row_stride( W );
	cs_W       = FLA_Obj_col_stride( W );

	rs_U       = FLA_Obj_row_stride( U );
	cs_U       = FLA_Obj_col_stride( U );

	rs_V       = FLA_Obj_row_stride( V );
	cs_V       = FLA_Obj_col_stride( V );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_H = FLA_COMPLEX_PTR( H );
			float*    buff_RG = FLA_FLOAT_PTR( RG );
			float*    buff_RH = FLA_FLOAT_PTR( RH );
			float*    buff_W = FLA_FLOAT_PTR( W );
			float*    buff_U = FLA_FLOAT_PTR( U );
			float*    buff_V = FLA_FLOAT_PTR( V );

			r_val = FLA_Bsvd_v_ops_var2( min( m_U, m_V ),
                                                     m_U,
			                             m_V,
			                             n_GH,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_H, rs_H, cs_H,
			                             buff_RG, rs_RG, cs_RG,
			                             buff_RH, rs_RH, cs_RH,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             buff_V, rs_V, cs_V,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_H = FLA_DOUBLE_COMPLEX_PTR( H );
			double*   buff_RG = FLA_DOUBLE_PTR( RG );
			double*   buff_RH = FLA_DOUBLE_PTR( RH );
			double*   buff_W = FLA_DOUBLE_PTR( W );
			double*   buff_U = FLA_DOUBLE_PTR( U );
			double*   buff_V = FLA_DOUBLE_PTR( V );

			r_val = FLA_Bsvd_v_opd_var2( min( m_U, m_V ),
                                                     m_U,
			                             m_V,
			                             n_GH,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_H, rs_H, cs_H,
			                             buff_RG, rs_RG, cs_RG,
			                             buff_RH, rs_RH, cs_RH,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             buff_V, rs_V, cs_V,
			                             b_alg );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_H = FLA_COMPLEX_PTR( H );
			float*    buff_RG = FLA_FLOAT_PTR( RG );
			float*    buff_RH = FLA_FLOAT_PTR( RH );
			scomplex* buff_W = FLA_COMPLEX_PTR( W );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );
			scomplex* buff_V = FLA_COMPLEX_PTR( V );

			r_val = FLA_Bsvd_v_opc_var2( min( m_U, m_V ),
                                                     m_U,
			                             m_V,
			                             n_GH,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_H, rs_H, cs_H,
			                             buff_RG, rs_RG, cs_RG,
			                             buff_RH, rs_RH, cs_RH,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             buff_V, rs_V, cs_V,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_H = FLA_DOUBLE_COMPLEX_PTR( H );
			double*   buff_RG = FLA_DOUBLE_PTR( RG );
			double*   buff_RH = FLA_DOUBLE_PTR( RH );
			dcomplex* buff_W = FLA_DOUBLE_COMPLEX_PTR( W );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );
			dcomplex* buff_V = FLA_DOUBLE_COMPLEX_PTR( V );

			r_val = FLA_Bsvd_v_opz_var2( min( m_U, m_V ),
                                                     m_U,
			                             m_V,
			                             n_GH,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_G, rs_G, cs_G,
			                             buff_H, rs_H, cs_H,
			                             buff_RG, rs_RG, cs_RG,
			                             buff_RH, rs_RH, cs_RH,
			                             buff_W, rs_W, cs_W,
			                             buff_U, rs_U, cs_U,
			                             buff_V, rs_V, cs_V,
			                             b_alg );

			break;
		}
	}

	return r_val;
}



FLA_Error FLA_Bsvd_v_ops_var2( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_H, int rs_H, int cs_H,
                               float*    buff_RG, int rs_RG, int cs_RG,
                               float*    buff_RH, int rs_RH, int cs_RH,
                               float*    buff_W, int rs_W, int cs_W,
                               float*    buff_U, int rs_U, int cs_U,
                               float*    buff_V, int rs_V, int cs_V,
                               int       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_v_opd_var2( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_H, int rs_H, int cs_H,
                               double*   buff_RG, int rs_RG, int cs_RG,
                               double*   buff_RH, int rs_RH, int cs_RH,
                               double*   buff_W, int rs_W, int cs_W,
                               double*   buff_U, int rs_U, int cs_U,
                               double*   buff_V, int rs_V, int cs_V,
                               int       b_alg )
{
	dcomplex  one        = bli_z1();
	double    rone       = bli_d1();
	double    rzero      = bli_d0();

	int       maxitr     = 6;

	double    eps;
	double    tolmul;
	double    tol;
	double    thresh;

	dcomplex* G;
	dcomplex* H;
	double*   d1;
	double*   e1;
	int       r_val;
	int       done;
	int       m_GH_sweep_max;
	int       ij_begin;
	int       ijTL, ijBR;
	int       m_A11;
	int       n_iter_perf;
	int       n_UV_apply;
	int       total_deflations;
	int       n_deflations;
	int       n_iter_prev;
	int       n_iter_perf_sweep_max;

	// Compute some convergence constants.
	eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	tolmul = max( 10.0, min( 100.0, pow( eps, -0.125 ) ) );
	FLA_Bsvd_compute_tol_thresh_opd( min_m_n,
	                                 tolmul,
	                                 maxitr,
	                                 buff_d, inc_d,
	                                 buff_e, inc_e,
	                                 &tol,
	                                 &thresh );

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// and H that we would need to initialize for the next sweep.
	m_GH_sweep_max = min_m_n - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Initialize RG and RH to identity.
	bli_dident( min_m_n,
	            buff_RG, rs_RG, cs_RG );
	bli_dident( min_m_n,
	            buff_RH, rs_RH, cs_RH );

	// Iterate until the matrix has completely deflated.
	for ( total_deflations = 0; done != TRUE; )
	{

		// Initialize G and H to contain only identity rotations.
		bli_zsetm( m_GH_sweep_max,
		           n_GH,
		           &one,
		           buff_G, rs_G, cs_G );
		bli_zsetm( m_GH_sweep_max,
		           n_GH,
		           &one,
		           buff_H, rs_H, cs_H );

		// Keep track of the maximum number of iterations performed in the
		// current sweep. This is used when applying the sweep's Givens
		// rotations.
		n_iter_perf_sweep_max = 0;

		// Perform a sweep: Move through the matrix and perform a bidiagonal
		// SVD on each non-zero submatrix that is encountered. During the
		// first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
		for ( ij_begin = 0; ij_begin < min_m_n;  )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Bsvd_v_opd_var2: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first submatrix along the diagonal that is
			// bounded by zeroes (or endpoints of the matrix). If no
			// submatrix is found (ie: if the entire superdiagonal is zero
			// then FLA_FAILURE is returned. This function also inspects
			// superdiagonal elements for proximity to zero. If a given
			// element is close enough to zero, then it is deemed
			// converged and manually set to zero.
			r_val = FLA_Bsvd_find_submatrix_opd( min_m_n,
			                                     ij_begin,
			                                     buff_d, inc_d,
			                                     buff_e, inc_e,
			                                     &ijTL,
			                                     &ijBR );

			// Verify that a submatrix was found. If one was not found,
			// then we are done with the current sweep. Furthermore, if
			// a submatrix was not found AND we began our search at the
			// beginning of the matrix (ie: ij_begin == 0), then the
			// matrix has completely deflated and so we are done with
			// Francis step iteration.
			if ( r_val == FLA_FAILURE )
			{
				if ( ij_begin == 0 )
				{
#ifdef PRINTF
printf( "FLA_Bsvd_v_opd_var2: superdiagonal is completely zero.\n" );
printf( "FLA_Bsvd_v_opd_var2: Francis iteration is done!\n" );
#endif
					done = TRUE;
				}

				// Break out of the current sweep so we can apply the last
				// remaining Givens rotations.
				break;
			}

			// If we got this far, then:
			//   (a) ijTL refers to the index of the first non-zero
			//       superdiagonal along the diagonal, and
			//   (b) ijBR refers to either:
			//       - the first zero element that occurs after ijTL, or
			//       - the the last diagonal element.
			// Note that ijTL and ijBR also correspond to the first and
			// last diagonal elements of the submatrix of interest. Thus,
			// we may compute the dimension of this submatrix as:
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Bsvd_v_opd_var2: ij_begin = %d\n", ij_begin );
printf( "FLA_Bsvd_v_opd_var2: ijTL     = %d\n", ijTL );
printf( "FLA_Bsvd_v_opd_var2: ijBR     = %d\n", ijBR );
printf( "FLA_Bsvd_v_opd_var2: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next submatrix
			// search in the current sweep.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;
			H  = buff_H + ijTL * rs_H;

			// Search for a batch of singular values, recursing on deflated
			// subproblems whenever possible. A new singular value search is
			// performed as long as
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Bsvd_iteracc_v_opd_var1( m_A11,
			                                            n_GH,
			                                            ijTL,
			                                            tol,
			                                            thresh,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            G,  rs_G, cs_G,
			                                            H,  rs_H, cs_H,
			                                            &n_iter_perf );

			// Record the number of deflations that occurred.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Bsvd_v_opd_var2: deflations observed       = %d\n", n_deflations );
printf( "FLA_Bsvd_v_opd_var2: total deflations observed = %d\n", total_deflations );
printf( "FLA_Bsvd_v_opd_var2: num iterations performed  = %d\n", n_iter_perf );
#endif

			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G and H we can apply and safely include all
			// non-identity rotations that were computed during the
			// singular value searches.
			m_GH_sweep_max = ijBR;

		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.

		// Recall that the number of columns of U and V to which we apply
		// rotations is one more than the number of rotations.
		n_UV_apply = m_GH_sweep_max + 1;


#ifdef PRINTF
printf( "FLA_Bsvd_v_opd_var2: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
printf( "FLA_Bsvd_v_opd_var2: m_U = %d\n", m_U );
printf( "FLA_Bsvd_v_opd_var2: m_V = %d\n", m_V );
printf( "FLA_Bsvd_v_opd_var2: napp= %d\n", n_UV_apply );
#endif

		// Apply the Givens rotations. Note that we only apply k sets of
		// rotations, where k = n_iter_perf_sweep_max. Also note that we only
		// apply to n_UV_apply columns of U and V since this is the most we
		// need to touch given the most recent value stored to m_GH_sweep_max.
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          min_m_n,
		                          n_UV_apply,
		                          n_iter_prev,
		                          buff_G,  rs_G,  cs_G,
		                          buff_RG, rs_RG, cs_RG,
		                          b_alg );
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          min_m_n,
		                          n_UV_apply,
		                          n_iter_prev,
		                          buff_H,  rs_H,  cs_H,
		                          buff_RH, rs_RH, cs_RH,
		                          b_alg );

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
printf( "FLA_Bsvd_v_opd_var2: total number of iterations performed: %d\n", n_iter_prev );
#endif
	}

	// Copy the contents of Q to temporary storage.
	bli_dcopymt( BLIS_NO_TRANSPOSE,
	             m_U,
	             m_V,
	             buff_U, rs_U, cs_U,
	             buff_W, rs_W, cs_W );
// W needs to be max_m_n-by-min_m_n!!!!!!!!!!!!!!!

	// Multiply U by R, overwriting U.
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE,
	           m_U,
	           m_V,
	           m_V,
	           &rone,
	           ( double* )buff_W,  rs_W,  cs_W,
	                      buff_RG, rs_RG, cs_RG,
	           &rzero,
	           ( double* )buff_U,  rs_U,  cs_U );

	bli_dcopymt( BLIS_NO_TRANSPOSE,
	             m_V,
	             m_V,
	             buff_V, rs_V, cs_V,
	             buff_W, rs_W, cs_W );

	// Multiply V by R, overwriting V.
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE,
	           m_V,
	           m_V,
	           m_V,
	           &rone,
	           ( double* )buff_W,  rs_W,  cs_W,
	                      buff_RH, rs_RH, cs_RH,
	           &rzero,
	           ( double* )buff_V,  rs_V,  cs_V );

	// Make all the singular values positive.
	{
		int    i;
		double minus_one = bli_dm1();

		for ( i = 0; i < min_m_n; ++i )
		{
			if ( buff_d[ (i  )*inc_d ] < rzero )
			{
				buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];
				
				// Scale the right singular vectors.
				bli_dscalv( BLIS_NO_CONJUGATE,
				            m_V,
				            &minus_one,
				            buff_V + (i  )*cs_V, rs_V );
			}
		}
	}

	return n_iter_prev;
}

FLA_Error FLA_Bsvd_v_opc_var2( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_H, int rs_H, int cs_H,
                               float*    buff_RG, int rs_RG, int cs_RG,
                               float*    buff_RH, int rs_RH, int cs_RH,
                               scomplex* buff_W, int rs_W, int cs_W,
                               scomplex* buff_U, int rs_U, int cs_U,
                               scomplex* buff_V, int rs_V, int cs_V,
                               int       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

FLA_Error FLA_Bsvd_v_opz_var2( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_H, int rs_H, int cs_H,
                               double*   buff_RG, int rs_RG, int cs_RG,
                               double*   buff_RH, int rs_RH, int cs_RH,
                               dcomplex* buff_W, int rs_W, int cs_W,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               dcomplex* buff_V, int rs_V, int cs_V,
                               int       b_alg )
{
	dcomplex  one        = bli_z1();
	double    rone       = bli_d1();
	double    rzero      = bli_d0();

	int       maxitr     = 6;

	double    eps;
	double    tolmul;
	double    tol;
	double    thresh;

	dcomplex* G;
	dcomplex* H;
	double*   d1;
	double*   e1;
	int       r_val;
	int       done;
	int       m_GH_sweep_max;
	int       ij_begin;
	int       ijTL, ijBR;
	int       m_A11;
	int       n_iter_perf;
	int       n_UV_apply;
	int       total_deflations;
	int       n_deflations;
	int       n_iter_prev;
	int       n_iter_perf_sweep_max;

	// Compute some convergence constants.
	eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	tolmul = max( 10.0, min( 100.0, pow( eps, -0.125 ) ) );
	FLA_Bsvd_compute_tol_thresh_opd( min_m_n,
	                                 tolmul,
	                                 maxitr,
	                                 buff_d, inc_d,
	                                 buff_e, inc_e,
	                                 &tol,
	                                 &thresh );

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// and H that we would need to initialize for the next sweep.
	m_GH_sweep_max = min_m_n - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Initialize RG and RH to identity.
	bli_dident( min_m_n,
	            buff_RG, rs_RG, cs_RG );
	bli_dident( min_m_n,
	            buff_RH, rs_RH, cs_RH );

	// Iterate until the matrix has completely deflated.
	for ( total_deflations = 0; done != TRUE; )
	{

		// Initialize G and H to contain only identity rotations.
		bli_zsetm( m_GH_sweep_max,
		           n_GH,
		           &one,
		           buff_G, rs_G, cs_G );
		bli_zsetm( m_GH_sweep_max,
		           n_GH,
		           &one,
		           buff_H, rs_H, cs_H );

		// Keep track of the maximum number of iterations performed in the
		// current sweep. This is used when applying the sweep's Givens
		// rotations.
		n_iter_perf_sweep_max = 0;

		// Perform a sweep: Move through the matrix and perform a bidiagonal
		// SVD on each non-zero submatrix that is encountered. During the
		// first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
		for ( ij_begin = 0; ij_begin < min_m_n;  )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Bsvd_v_opz_var2: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first submatrix along the diagonal that is
			// bounded by zeroes (or endpoints of the matrix). If no
			// submatrix is found (ie: if the entire superdiagonal is zero
			// then FLA_FAILURE is returned. This function also inspects
			// superdiagonal elements for proximity to zero. If a given
			// element is close enough to zero, then it is deemed
			// converged and manually set to zero.
			r_val = FLA_Bsvd_find_submatrix_opd( min_m_n,
			                                     ij_begin,
			                                     buff_d, inc_d,
			                                     buff_e, inc_e,
			                                     &ijTL,
			                                     &ijBR );

			// Verify that a submatrix was found. If one was not found,
			// then we are done with the current sweep. Furthermore, if
			// a submatrix was not found AND we began our search at the
			// beginning of the matrix (ie: ij_begin == 0), then the
			// matrix has completely deflated and so we are done with
			// Francis step iteration.
			if ( r_val == FLA_FAILURE )
			{
				if ( ij_begin == 0 )
				{
#ifdef PRINTF
printf( "FLA_Bsvd_v_opz_var2: superdiagonal is completely zero.\n" );
printf( "FLA_Bsvd_v_opz_var2: Francis iteration is done!\n" );
#endif
					done = TRUE;
				}

				// Break out of the current sweep so we can apply the last
				// remaining Givens rotations.
				break;
			}

			// If we got this far, then:
			//   (a) ijTL refers to the index of the first non-zero
			//       superdiagonal along the diagonal, and
			//   (b) ijBR refers to either:
			//       - the first zero element that occurs after ijTL, or
			//       - the the last diagonal element.
			// Note that ijTL and ijBR also correspond to the first and
			// last diagonal elements of the submatrix of interest. Thus,
			// we may compute the dimension of this submatrix as:
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Bsvd_v_opz_var2: ij_begin = %d\n", ij_begin );
printf( "FLA_Bsvd_v_opz_var2: ijTL     = %d\n", ijTL );
printf( "FLA_Bsvd_v_opz_var2: ijBR     = %d\n", ijBR );
printf( "FLA_Bsvd_v_opz_var2: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next submatrix
			// search in the current sweep.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			G  = buff_G + ijTL * rs_G;
			H  = buff_H + ijTL * rs_H;

			// Search for a batch of singular values, recursing on deflated
			// subproblems whenever possible. A new singular value search is
			// performed as long as
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Bsvd_iteracc_v_opd_var1( m_A11,
			                                            n_GH,
			                                            ijTL,
			                                            tol,
			                                            thresh,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            G,  rs_G, cs_G,
			                                            H,  rs_H, cs_H,
			                                            &n_iter_perf );

			// Record the number of deflations that occurred.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Bsvd_v_opz_var2: deflations observed       = %d\n", n_deflations );
printf( "FLA_Bsvd_v_opz_var2: total deflations observed = %d\n", total_deflations );
printf( "FLA_Bsvd_v_opz_var2: num iterations performed  = %d\n", n_iter_perf );
#endif

			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G and H we can apply and safely include all
			// non-identity rotations that were computed during the
			// singular value searches.
			m_GH_sweep_max = ijBR;

		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.

		// Recall that the number of columns of U and V to which we apply
		// rotations is one more than the number of rotations.
		n_UV_apply = m_GH_sweep_max + 1;


#ifdef PRINTF
printf( "FLA_Bsvd_v_opz_var2: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
printf( "FLA_Bsvd_v_opz_var2: m_U = %d\n", m_U );
printf( "FLA_Bsvd_v_opz_var2: m_V = %d\n", m_V );
printf( "FLA_Bsvd_v_opz_var2: napp= %d\n", n_UV_apply );
#endif

		// Apply the Givens rotations. Note that we only apply k sets of
		// rotations, where k = n_iter_perf_sweep_max. Also note that we only
		// apply to n_UV_apply columns of U and V since this is the most we
		// need to touch given the most recent value stored to m_GH_sweep_max.
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          min_m_n,
		                          n_UV_apply,
		                          n_iter_prev,
		                          buff_G,  rs_G,  cs_G,
		                          buff_RG, rs_RG, cs_RG,
		                          b_alg );
		//FLA_Apply_G_rf_bld_var5b( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9b( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6b( n_iter_perf_sweep_max,
		                          min_m_n,
		                          n_UV_apply,
		                          n_iter_prev,
		                          buff_H,  rs_H,  cs_H,
		                          buff_RH, rs_RH, cs_RH,
		                          b_alg );

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
printf( "FLA_Bsvd_v_opz_var2: total number of iterations performed: %d\n", n_iter_prev );
#endif
	}

	// Copy the contents of Q to temporary storage.
	bli_zcopymt( BLIS_NO_TRANSPOSE,
	             m_U,
	             m_V,
	             buff_U, rs_U, cs_U,
	             buff_W, rs_W, cs_W );
// W needs to be max_m_n-by-min_m_n!!!!!!!!!!!!!!!

	// Multiply U by R, overwriting U.
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE,
	           2*m_U,
	           m_V,
	           m_V,
	           &rone,
	           ( double* )buff_W,  rs_W,  2*cs_W,
	                      buff_RG, rs_RG,   cs_RG,
	           &rzero,
	           ( double* )buff_U,  rs_U,  2*cs_U );

	bli_zcopymt( BLIS_NO_TRANSPOSE,
	             m_V,
	             m_V,
	             buff_V, rs_V, cs_V,
	             buff_W, rs_W, cs_W );

	// Multiply V by R, overwriting V.
	bli_dgemm( BLIS_NO_TRANSPOSE,
	           BLIS_NO_TRANSPOSE,
	           2*m_V,
	           m_V,
	           m_V,
	           &rone,
	           ( double* )buff_W,  rs_W,  2*cs_W,
	                      buff_RH, rs_RH,   cs_RH,
	           &rzero,
	           ( double* )buff_V,  rs_V,  2*cs_V );

	// Make all the singular values positive.
	{
		int    i;
		double minus_one = bli_dm1();

		for ( i = 0; i < min_m_n; ++i )
		{
			if ( buff_d[ (i  )*inc_d ] < rzero )
			{
				buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];
				
				// Scale the right singular vectors.
				bli_zdscalv( BLIS_NO_CONJUGATE,
				             m_V,
				             &minus_one,
				             buff_V + (i  )*cs_V, rs_V );
			}
		}
	}

	return n_iter_prev;
}

