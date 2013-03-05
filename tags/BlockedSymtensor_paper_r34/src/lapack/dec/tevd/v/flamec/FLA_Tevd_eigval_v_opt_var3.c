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


FLA_Error FLA_Tevd_eigval_v_ops_var3( int       m_A,
                                      int       m_U,
                                      int       n_G,
                                      scomplex* buff_G, int rs_G, int cs_G,
                                      float*    buff_d, int inc_d, 
                                      float*    buff_e, int inc_e,
                                      float*    buff_l, int inc_l,
                                      int*      buff_ls, int inc_ls,
                                      float*    buff_pu, int inc_pu,
                                      int*      n_iter )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_eigval_v_opd_var3( int       m_A,
                                      int       m_U,
                                      int       n_G,
                                      dcomplex* buff_G, int rs_G, int cs_G,
                                      double*   buff_d, int inc_d, 
                                      double*   buff_e, int inc_e,
                                      double*   buff_l, int inc_l,
                                      int*      buff_ls, int inc_ls,
                                      double*   buff_pu, int inc_pu,
                                      int*      n_iter )
{
	FLA_Error r_val;
	double    eps;
	double    safmin;
	double*   e_last;
	double*   d_last;
	double*   d_last_m1;
	double    shift;
	int       ij_shift;
	int       k;
	int       n_iter_allowed = n_G;

	// Query epsilon and safmin, which are used in the test for convergence.
	eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
	safmin = FLA_Mach_params_opd( FLA_MACH_SFMIN );

	// Initialize a pointer to the last sub-diagonal element and two
	// more to the last and second last
	e_last    = &buff_e[ (m_A-2)*inc_e ];
	d_last_m1 = &buff_d[ (m_A-2)*inc_d ];
	d_last    = &buff_d[ (m_A-1)*inc_d ];


	for ( k = 0; k < n_iter_allowed; ++k )
	{
		dcomplex* g1 = buff_G + (k  )*cs_G;

		/*------------------------------------------------------------*/

		// If we've converged, record k and return index of eigenvalue found.
		// The reason we check before the Francis step (rather than after)
		// is so we correctly handle situations where the last diagonal
		// element has already converged from previous eigenvalue searches
		// and thus no iteration is necessary. If we checked after the
		// Francis step, we would have unnecessarily executed an additional
		// Francis step's worth of rotations with a sub-optimal shift (since
		// it would be using a 2x2 that was not "centered" properly).
		if ( MAC_Tevd_eigval_converged_opd( eps, safmin, *d_last_m1, *e_last, *d_last ) )
		{
			*e_last = 0.0;
			*n_iter = k;
			return m_A - 1;
		}

		FLA_Tevd_find_perfshift_opd( m_A,
		                             m_U,
		                             buff_d, inc_d,
		                             buff_e, inc_e,
		                             buff_l, inc_l,
		                             buff_ls, inc_ls,
		                             buff_pu, inc_pu,
		                             &ij_shift );

		if ( ij_shift == -1 )
		{
			FLA_Wilkshift_tridiag_opd( *d_last_m1,
			                           *e_last,
			                           *d_last,
    			                       &shift );
		}
		else
		{
			shift = buff_l[ ij_shift*inc_l ];
		}

		// Perform a Francis step.
		r_val = FLA_Tevd_francis_v_opd_var1( m_A,
		                                     &shift,
		                                     g1,     rs_G,
		                                     buff_d, inc_d,
		                                     buff_e, inc_e );

		if ( ij_shift >= 0 &&
		     MAC_Tevd_eigval_converged_opd( eps, safmin, *d_last_m1, *e_last, *d_last ) )
		{
			buff_ls[ ij_shift * inc_ls ] = 1;
			*e_last = 0.0;
			*n_iter = k + 1;
			return m_A - 1;
		}

		// Check for internal deflation.
		if ( r_val != FLA_SUCCESS )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_eigval_v_opt_var3: Internal deflation in col %d, eig %d\n", r_val, m_A - 1 );
			printf( "FLA_Tevd_eigval_v_opt_var3: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
			printf( "FLA_Tevd_eigval_v_opt_var3: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif

			// Set the off-diagonal element to zero.
			buff_e[ r_val*inc_e ] = 0.0;

			*n_iter = k + 1;
			return r_val;
		}

		/*------------------------------------------------------------*/
	}

	*n_iter = n_iter_allowed;
	return FLA_FAILURE;
}

