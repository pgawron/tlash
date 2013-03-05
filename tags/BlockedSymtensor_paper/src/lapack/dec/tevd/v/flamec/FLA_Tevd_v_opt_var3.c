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

FLA_Error FLA_Tevd_v_opt_var3( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj ls, FLA_Obj pu, FLA_Obj G, FLA_Obj U, dim_t b_alg )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Datatype datatype;
	int          m_A, m_U, n_G;
	int          inc_d;
	int          inc_e;
	int          inc_l;
	int          inc_ls;
	int          inc_pu;
	int          rs_G, cs_G;
	int          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A       = FLA_Obj_vector_dim( d );
	m_U       = FLA_Obj_length( U );
	n_G       = FLA_Obj_width( G );

	inc_d     = FLA_Obj_vector_inc( d );
	inc_e     = FLA_Obj_vector_inc( e );
	inc_l     = FLA_Obj_vector_inc( l );
	inc_ls    = FLA_Obj_vector_inc( ls );
	inc_pu    = FLA_Obj_vector_inc( pu );
	
	rs_G      = FLA_Obj_row_stride( G );
	cs_G      = FLA_Obj_col_stride( G );

	rs_U      = FLA_Obj_row_stride( U );
	cs_U      = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			float*    buff_l = FLA_FLOAT_PTR( l );
			int*      buff_ls= FLA_INT_PTR( ls );
			float*    buff_pu= FLA_FLOAT_PTR( pu );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			float*    buff_U = FLA_FLOAT_PTR( U );

			r_val = FLA_Tevd_v_ops_var3( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_l, inc_l,
			                             buff_ls, inc_ls,
			                             buff_pu, inc_pu,
			                             buff_G, rs_G, cs_G,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			double*   buff_l = FLA_DOUBLE_PTR( l );
			int*      buff_ls= FLA_INT_PTR( ls );
			double*   buff_pu= FLA_DOUBLE_PTR( pu );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			r_val = FLA_Tevd_v_opd_var3( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_l, inc_l,
			                             buff_ls, inc_ls,
			                             buff_pu, inc_pu,
			                             buff_G, rs_G, cs_G,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_COMPLEX:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			float*    buff_l = FLA_FLOAT_PTR( l );
			int*      buff_ls= FLA_INT_PTR( ls );
			float*    buff_pu= FLA_FLOAT_PTR( pu );
			scomplex* buff_G = FLA_COMPLEX_PTR( G );
			scomplex* buff_U = FLA_COMPLEX_PTR( U );

			r_val = FLA_Tevd_v_opc_var3( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_l, inc_l,
			                             buff_ls, inc_ls,
			                             buff_pu, inc_pu,
			                             buff_G, rs_G, cs_G,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			double*   buff_l = FLA_DOUBLE_PTR( l );
			int*      buff_ls= FLA_INT_PTR( ls );
			double*   buff_pu= FLA_DOUBLE_PTR( pu );
			dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
			dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );

			r_val = FLA_Tevd_v_opz_var3( m_A,
			                             m_U,
			                             n_G,
			                             n_iter_max,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_l, inc_l,
			                             buff_ls, inc_ls,
			                             buff_pu, inc_pu,
			                             buff_G, rs_G, cs_G,
			                             buff_U, rs_U, cs_U,
			                             b_alg );

			break;
		}
	}

	return r_val;
}



FLA_Error FLA_Tevd_v_ops_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
                               int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opd_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
                               int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg )
{
	dcomplex  one  = bli_z1();

	dcomplex* G;
	double*   d1;
	double*   e1;
	double*   pu1;
	int       r_val;
	int       done;
	int       m_G_sweep_max;
	int       ij_begin;
	int       ijTL, ijBR;
	int       m_A11;
	int       n_iter_perf;
	int       n_U_apply;
	int       total_deflations;
	int       n_deflations;
	int       n_iter_prev;
	int       n_iter_perf_sweep_max;

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Iterate until the matrix has completely deflated.
	for ( total_deflations = 0; done != TRUE; )
	{

		// Initialize G to contain only identity rotations.
		bli_zsetm( m_G_sweep_max,
		           n_G,
		           &one,
		           buff_G, rs_G, cs_G );

		// Keep track of the maximum number of iterations performed in the
		// current sweep. This is used when applying the sweep's Givens
		// rotations.
		n_iter_perf_sweep_max = 0;

		// Perform a sweep: Move through the matrix and perform a tridiagonal
		// EVD on each non-zero submatrix that is encountered. During the
		// first time through, ijTL will be 0 and ijBR will be m_A - 1.
		for ( ij_begin = 0; ij_begin < m_A; )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_v_opd_var3: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first submatrix along the diagonal that is
			// bounded by zeroes (or endpoints of the matrix). If no
			// submatrix is found (ie: if the entire subdiagonal is zero
			// then FLA_FAILURE is returned. This function also inspects
			// subdiagonal elements for proximity to zero. If a given
			// element is close enough to zero, then it is deemed
			// converged and manually set to zero.
			r_val = FLA_Tevd_find_submatrix_opd( m_A,
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
printf( "FLA_Tevd_v_opd_var3: subdiagonal is completely zero.\n" );
printf( "FLA_Tevd_v_opd_var3: Francis iteration is done!\n" );
#endif
					done = TRUE;
				}

				// Break out of the current sweep so we can apply the last
				// remaining Givens rotations.
				break;
			}

			// If we got this far, then:
			//   (a) ijTL refers to the index of the first non-zero
			//       subdiagonal along the diagonal, and
			//   (b) ijBR refers to either:
			//       - the first zero element that occurs after ijTL, or
			//       - the the last diagonal element.
			// Note that ijTL and ijBR also correspond to the first and
			// last diagonal elements of the submatrix of interest. Thus,
			// we may compute the dimension of this submatrix as:
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var3: ij_begin = %d\n", ij_begin );
printf( "FLA_Tevd_v_opd_var3: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_v_opd_var3: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_v_opd_var3: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next submatrix
			// search in the current sweep.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			pu1 = buff_pu + ijTL * inc_pu;
			G  = buff_G + ijTL * rs_G;

			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever a split occurs. Iteration continues
			// as long as:
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Tevd_iteracc_v_opd_var3( m_A11,
			                                            m_A,
			                                            n_G,
			                                            ijTL,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            buff_l,  inc_l,
			                                            buff_ls,  inc_ls,
			                                            pu1, inc_pu,
			                                            G,  rs_G, cs_G,
			                                            &n_iter_perf );

			// Record the number of deflations that were observed.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var3: deflations observed       = %d\n", n_deflations );
printf( "FLA_Tevd_v_opd_var3: total deflations observed = %d\n", total_deflations );
printf( "FLA_Tevd_v_opd_var3: num iterations performed  = %d\n", n_iter_perf );
#endif

			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G we can apply and safely include all
			// non-identity rotations that were computed during the
			// eigenvalue searches.
			m_G_sweep_max = ijBR;

			// Make sure we haven't exceeded our maximum iteration count.
			if ( n_iter_prev >= m_A * n_iter_max )
			{
#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var3: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
				FLA_Abort();
				//return FLA_FAILURE;
			}
		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.

		// Recall that the number of columns of U to which we apply
		// rotations is one more than the number of rotations.
		n_U_apply = m_G_sweep_max + 1;

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var3: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
#endif

		// Apply the Givens rotations. Note that we optimize the scope
		// of the operation in two ways:
		//   1. We only apply k sets of Givens rotations, where
		//      k = n_iter_perf_sweep_max. We could simply always apply
		//      n_G sets of rotations since G is initialized to contain
		//      identity rotations in every element, but we do this to
		//      save a little bit of time.
		//   2. We only apply to the first n_U_apply columns of A since
		//      this is the most we need to touch given the ijBR index
		//      bound of the last submatrix found in the previous sweep.
		//      Similar to above, we could simply always perform the
		//      application on all m_A columns of A, but instead we apply
		//      only to the first n_U_apply columns to save time.
		//FLA_Apply_G_rf_bld_var5( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_bld_var3( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var9( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_bld_var6( n_iter_perf_sweep_max,
		                         m_U,
		                         n_U_apply,
		                         buff_G, rs_G, cs_G,
		                         buff_U, rs_U, cs_U,
		                         b_alg );

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
printf( "FLA_Tevd_v_opd_var3: total number of iterations performed: %d\n", n_iter_prev );
#endif
	}

	return n_iter_prev;
}

FLA_Error FLA_Tevd_v_opc_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               float*    buff_l, int inc_l,
                               int*      buff_ls, int inc_ls,
                               float*    buff_pu, int inc_pu,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}

FLA_Error FLA_Tevd_v_opz_var3( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               double*   buff_l, int inc_l,
                               int*      buff_ls, int inc_ls,
                               double*   buff_pu, int inc_pu,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg )
{
	dcomplex  one  = bli_z1();

	dcomplex* G;
	double*   d1;
	double*   e1;
	double*   pu1;
	int       r_val;
	int       done;
	int       m_G_sweep_max;
	int       ij_begin;
	int       ijTL, ijBR;
	int       m_A11;
	int       n_iter_perf;
	int       n_U_apply;
	int       total_deflations;
	int       n_deflations;
	int       n_iter_prev;
	int       n_iter_perf_sweep_max;

	// Initialize our completion flag.
	done = FALSE;

	// Initialize a counter that holds the maximum number of rows of G
	// that we would need to initialize for the next sweep.
	m_G_sweep_max = m_A - 1;

	// Initialize a counter for the total number of iterations performed.
	n_iter_prev = 0;

	// Iterate until the matrix has completely deflated.
	for ( total_deflations = 0; done != TRUE; )
	{

		// Initialize G to contain only identity rotations.
		bli_zsetm( m_G_sweep_max,
		           n_G,
		           &one,
		           buff_G, rs_G, cs_G );

		// Keep track of the maximum number of iterations performed in the
		// current sweep. This is used when applying the sweep's Givens
		// rotations.
		n_iter_perf_sweep_max = 0;

		// Perform a sweep: Move through the matrix and perform a tridiagonal
		// EVD on each non-zero submatrix that is encountered. During the
		// first time through, ijTL will be 0 and ijBR will be m_A - 1.
		for ( ij_begin = 0; ij_begin < m_A; )
		{

#ifdef PRINTF
if ( ij_begin == 0 )
printf( "FLA_Tevd_v_opz_var3: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

			// Search for the first submatrix along the diagonal that is
			// bounded by zeroes (or endpoints of the matrix). If no
			// submatrix is found (ie: if the entire subdiagonal is zero
			// then FLA_FAILURE is returned. This function also inspects
			// subdiagonal elements for proximity to zero. If a given
			// element is close enough to zero, then it is deemed
			// converged and manually set to zero.
			r_val = FLA_Tevd_find_submatrix_opd( m_A,
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
printf( "FLA_Tevd_v_opz_var3: subdiagonal is completely zero.\n" );
printf( "FLA_Tevd_v_opz_var3: Francis iteration is done!\n" );
#endif
					done = TRUE;
				}

				// Break out of the current sweep so we can apply the last
				// remaining Givens rotations.
				break;
			}

			// If we got this far, then:
			//   (a) ijTL refers to the index of the first non-zero
			//       subdiagonal along the diagonal, and
			//   (b) ijBR refers to either:
			//       - the first zero element that occurs after ijTL, or
			//       - the the last diagonal element.
			// Note that ijTL and ijBR also correspond to the first and
			// last diagonal elements of the submatrix of interest. Thus,
			// we may compute the dimension of this submatrix as:
			m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: ij_begin = %d\n", ij_begin );
printf( "FLA_Tevd_v_opz_var3: ijTL     = %d\n", ijTL );
printf( "FLA_Tevd_v_opz_var3: ijBR     = %d\n", ijBR );
printf( "FLA_Tevd_v_opz_var3: m_A11    = %d\n", m_A11 );
#endif

			// Adjust ij_begin, which gets us ready for the next submatrix
			// search in the current sweep.
			ij_begin = ijBR + 1;

			// Index to the submatrices upon which we will operate.
			d1 = buff_d + ijTL * inc_d;
			e1 = buff_e + ijTL * inc_e;
			pu1 = buff_pu + ijTL * inc_pu;
			G  = buff_G + ijTL * rs_G;

			// Search for a batch of eigenvalues, recursing on deflated
			// subproblems whenever a split occurs. Iteration continues
			// as long as:
			//   (a) there is still matrix left to operate on, and
			//   (b) the number of iterations performed in this batch is
			//       less than n_G.
			// If/when either of the two above conditions fails to hold,
			// the function returns.
			n_deflations = FLA_Tevd_iteracc_v_opd_var3( m_A11,
			                                            m_A,
			                                            n_G,
			                                            ijTL,
			                                            d1, inc_d,
			                                            e1, inc_e,
			                                            buff_l,  inc_l,
			                                            buff_ls,  inc_ls,
			                                            pu1, inc_pu,
			                                            G,  rs_G, cs_G,
			                                            &n_iter_perf );

			// Record the number of deflations that were observed.
			total_deflations += n_deflations;

			// Update the maximum number of iterations performed in the
			// current sweep.
			n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: deflations observed       = %d\n", n_deflations );
printf( "FLA_Tevd_v_opz_var3: total deflations observed = %d\n", total_deflations );
printf( "FLA_Tevd_v_opz_var3: num iterations performed  = %d\n", n_iter_perf );
#endif

			// Store the most recent value of ijBR in m_G_sweep_max.
			// When the sweep is done, this value will contain the minimum
			// number of rows of G we can apply and safely include all
			// non-identity rotations that were computed during the
			// eigenvalue searches.
			m_G_sweep_max = ijBR;

			// Make sure we haven't exceeded our maximum iteration count.
			if ( n_iter_prev >= m_A * n_iter_max )
			{
#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
				FLA_Abort();
				//return FLA_FAILURE;
			}
		}

		// The sweep is complete. Now we must apply the Givens rotations
		// that were accumulated during the sweep.

		// Recall that the number of columns of U to which we apply
		// rotations is one more than the number of rotations.
		n_U_apply = m_G_sweep_max + 1;

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
#endif

		// Apply the Givens rotations. Note that we optimize the scope
		// of the operation in two ways:
		//   1. We only apply k sets of Givens rotations, where
		//      k = n_iter_perf_sweep_max. We could simply always apply
		//      n_G sets of rotations since G is initialized to contain
		//      identity rotations in every element, but we do this to
		//      save a little bit of time.
		//   2. We only apply to the first n_U_apply columns of A since
		//      this is the most we need to touch given the ijBR index
		//      bound of the last submatrix found in the previous sweep.
		//      Similar to above, we could simply always perform the
		//      application on all m_A columns of A, but instead we apply
		//      only to the first n_U_apply columns to save time.
		//FLA_Apply_G_rf_blz_var5( n_iter_perf_sweep_max,
		FLA_Apply_G_rf_blz_var3( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_blz_var9( n_iter_perf_sweep_max,
		//FLA_Apply_G_rf_blz_var6( n_iter_perf_sweep_max,
		                         m_U,
		                         n_U_apply,
		                         buff_G, rs_G, cs_G,
		                         buff_U, rs_U, cs_U,
		                         b_alg );

		// Increment the total number of iterations previously performed.
		n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
printf( "FLA_Tevd_v_opz_var3: total number of iterations performed: %d\n", n_iter_prev );
#endif
	}

	return n_iter_prev;
}

