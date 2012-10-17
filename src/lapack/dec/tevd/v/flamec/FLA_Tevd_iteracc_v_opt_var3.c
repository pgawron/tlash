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

FLA_Error FLA_Tevd_iteracc_v_ops_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       float*    buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       float*    buff_pu, int inc_pu,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_iteracc_v_opd_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       double*   buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       double*   buff_pu, int inc_pu,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf )
{
	FLA_Error r_val;
	int       i, k;
	int       k_iter       = 0;
	int       n_deflations = 0;
	//double    pshift;
	//double    eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	//double    safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );

	// Iterate from back to front until all that is left is a 2x2.
	for ( i = m_A - 1; i > 1; --i )
	{
		dcomplex* G1     = buff_G + (k_iter)*cs_G;
		int       m_ATL  = i + 1;
		int       k_left = n_G - k_iter;

		/*------------------------------------------------------------*/

		// Search for an eigenvalue of ATL submatrix until
		//   (a) deflation occurs, or
		//   (b) we perform the maximum number of additional iterations
		//       that are allowed within the current sweep
		//       (ie: n_G - k_iter).

		r_val = FLA_Tevd_eigval_v_opd_var3( m_ATL,
		                                    m_U,
		                                    k_left,
		                                    //1,
		                                    G1, rs_G, cs_G,
		                                    buff_d, inc_d,
		                                    buff_e, inc_e,
		                                    buff_l, inc_l,
		                                    buff_ls, inc_ls,
		                                    buff_pu, inc_pu,
		                                    &k );
/*
			FLA_Tevd_find_perfshift_opd( m_ATL,
			                             m_U,
			                             buff_d, inc_d,
			                             buff_e, inc_e,
			                             buff_l, inc_l,
			                             buff_ls, inc_ls,
			                             buff_pu, inc_pu,
			                             &pshift );

		for ( k = 0; k < k_left; ++k )
		{
			// Mark the shift as used.
			//pshift = *(buff_l + ij_shift * inc_l);
			//*(buff_ls + ij_shift * inc_ls) = 0;
			//printf( "using pshift %22.15e\n", pshift );
			//printf( "using pshift %f\n", pshift );

			r_val = FLA_Tevd_francis_v_opd_var1( m_ATL,
			                                     &pshift,
			                                     g1,     rs_G,
			                                     buff_d, inc_d,
			                                     buff_e, inc_e );
			g1 += cs_G;

			// Check for internal deflation.
			if ( r_val != FLA_SUCCESS )
			{
//#ifdef PRINTF
//				printf( "FLA_Tevd_eigval_v_opt_var1: Internal deflation in col %d, eig %d\n", r_val, m_A - 1 );
//				printf( "FLA_Tevd_eigval_v_opt_var1: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
//				printf( "FLA_Tevd_eigval_v_opt_var1: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
//#endif
	
				//printf( "found internal deflation in column %d\n", r_val );
				// Set the off-diagonal element to zero.
				buff_e[ r_val*inc_e ] = 0.0;
				break;
			}
			else
			{
				double e_last    = buff_e[ (m_ATL-2)*inc_e ];
				double d_last_m1 = buff_d[ (m_ATL-2)*inc_d ];
				double d_last    = buff_d[ (m_ATL-1)*inc_d ];
				r_val = i;

				if ( MAC_Tevd_eigval_converged_opd( eps, safmin, d_last_m1, e_last, d_last ) )
				{
					//printf( "zeroing %22.15e\n", buff_e[ (m_ATL-2)*inc_e ] );
					buff_e[ (m_ATL-2)*inc_e ] = 0.0;
					break;
				}
			}
		}
*/

		// If the eigenvalue search did not result in any deflation, return.
		if ( r_val == FLA_FAILURE && k_iter == n_G )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_iteracc_v_opd_var1: failed to converge (m_A11 = %d) after %2d iters k_total=%d/%d\n", i, k, k_iter, n_G );
#endif
			*n_iter_perf = k_iter;
			return n_deflations;
		}

		// Update local counters according to the results of the eigenvalue
		// search.
		k_iter       += k;
		n_deflations += 1;


#ifdef PRINTF
		if ( r_val == i )
			printf( "FLA_Tevd_iteracc_v_opd_var3: found eig %22.15e in col %3d (n=%d) after %2d iters  k_total=%d/%d\n", buff_d[ r_val*inc_d ], ijTL+r_val, m_ATL, k, k_iter, n_G );
		else
			printf( "FLA_Tevd_iteracc_v_opd_var3: split occurred in col %3d (n=%d) after %2d iters  k_total=%d/%d\n", ijTL+r_val, m_ATL, k, k_iter, n_G );
#endif

		// If the most recent eigenvalue search put us at our limit
		// for accumulated Givens rotation sets, return.
		if ( k_iter == n_G )
		{
			*n_iter_perf = k_iter;
			return n_deflations;
		}


		// If r_val != i, then a split occurred somewhere within submatrix
		// ATL. Therefore, we must recurse with two subproblems.
		if ( r_val != i )
		{
			int       m_TLr = r_val + 1;
			int       m_BRr = m_ATL - m_TLr;
			int       ijTLr = 0;
			int       ijBRr = m_TLr;
			int       n_Gr  = n_G - k_iter;
			double*   dTL   = buff_d + (0    )*inc_d;
			double*   eTL   = buff_e + (0    )*inc_e;
			double*   puTL  = buff_pu+ (0    )*inc_pu;
			dcomplex* GT    = buff_G + (0    )*rs_G + (k_iter)*cs_G;
			double*   dBR   = buff_d + (ijBRr)*inc_d;
			double*   eBR   = buff_e + (ijBRr)*inc_e;
			double*   puBR  = buff_pu+ (ijBRr)*inc_pu;
			dcomplex* GB    = buff_G + (ijBRr)*rs_G + (k_iter)*cs_G;

			int       n_deflationsTL;
			int       n_deflationsBR;
			int       n_iter_perfTL;
			int       n_iter_perfBR;

#ifdef PRINTF
printf( "FLA_Tevd_iteracc_v_opd_var3: Internal deflation in col %d\n", ijTL+r_val );
printf( "FLA_Tevd_iteracc_v_opd_var3: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
printf( "FLA_Tevd_iteracc_v_opd_var3: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif
#ifdef PRINTF
printf( "FLA_Tevd_iteracc_v_opd_var3: recursing: m_TLr m_BRr: %d %d\n", m_TLr, m_BRr );
printf( "FLA_Tevd_iteracc_v_opd_var3:            ijTLr ijBRr: %d %d\n", ijTLr, ijBRr );
printf( "FLA_Tevd_iteracc_v_opd_var3:            GB(0,0) i,j: %d %d\n", ijTL + m_TLr+1, k_iter );
#endif
			n_deflationsTL = FLA_Tevd_iteracc_v_opd_var3( m_TLr,
			                                              m_U,
			                                              n_Gr,
			                                              ijTL + ijTLr,
			                                              dTL, inc_d,
			                                              eTL, inc_e,
			                                              buff_l, inc_l,
			                                              buff_ls, inc_ls,
			                                              puTL, inc_pu,
			                                              GT,  rs_G, cs_G,
			                                              &n_iter_perfTL );
			n_deflationsBR = FLA_Tevd_iteracc_v_opd_var3( m_BRr,
			                                              m_U,
			                                              n_Gr,
			                                              ijTL + ijBRr,
			                                              dBR, inc_d,
			                                              eBR, inc_e,
			                                              buff_l, inc_l,
			                                              buff_ls, inc_ls,
			                                              puBR, inc_pu,
			                                              GB,  rs_G, cs_G,
			                                              &n_iter_perfBR );

			*n_iter_perf = k_iter + max( n_iter_perfTL, n_iter_perfBR );

#ifdef PRINTF
printf( "FLA_Tevd_iteracc_v_opd_var3: num deflations: %d = (prev:%d, TL:%d, BR:%d)\n", n_deflations + n_deflationsTL + n_deflationsBR, n_deflations, n_deflationsTL, n_deflationsBR );
printf( "FLA_Tevd_iteracc_v_opd_var3: num iterations: %d = (prev:%d, TL:%d, BR:%d)\n", *n_iter_perf, k_iter, n_iter_perfTL, n_iter_perfBR );
#endif
			return n_deflations + n_deflationsTL + n_deflationsBR;
		}

		/*------------------------------------------------------------*/
	}

	// Skip 1x1 matrices (and submatrices) entirely.
	if ( m_A > 1 )
	{
		dcomplex* g1      = buff_G + (k_iter)*cs_G;

		double*   alpha11 = buff_d + (0  )*inc_d;
		double*   alpha21 = buff_e + (0  )*inc_e;
		double*   alpha22 = buff_d + (1  )*inc_d;
		double    lambda1;
		double    lambda2;

		double    gamma;
		double    sigma;

		// Find the eigenvalue decomposition of the remaining (or only) 2x2
		// submatrix.
		FLA_Hevv_2x2_opd( alpha11,
		                  alpha21,
		                  alpha22,
		                  &lambda1,
		                  &lambda2,
		                  &gamma,
		                  &sigma );

		// Store the eigenvalues.
		*alpha11 = lambda1;
		*alpha22 = lambda2;

		// Zero out the remaining subdiagonal element.
		*alpha21 = 0.0;

		// Store the rotation.
		g1[0].real = gamma;
		g1[0].imag = sigma;


		// Update the local counters.
		k_iter       += 1;
		n_deflations += 1;

#ifdef PRINTF
printf( "FLA_Tevd_iteracc_v_opd_var3: Hevv  eig %22.15e in col %3d (n=%d) after %2d iters  k_total=%d/%d\n", buff_d[ 1*inc_d ], ijTL+1, 2, 1, k_iter, n_G );
printf( "FLA_Tevd_iteracc_v_opd_var3: Hevv  eig %22.15e in col %3d (n=%d) after %2d iters  k_total=%d/%d\n", buff_d[ 0*inc_d ], ijTL+0, 2, 0, k_iter, n_G );
#endif
	}


	*n_iter_perf = k_iter;
	return n_deflations;
}
