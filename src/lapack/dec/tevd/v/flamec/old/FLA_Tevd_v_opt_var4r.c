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

FLA_Error FLA_Tevd_v_opt_var4r( FLA_Obj d, FLA_Obj e, FLA_Obj C, FLA_Obj S, FLA_Obj U )
{
	FLA_Datatype datatype;
	int          m_A, m_U, n_CS;
	int          inc_d;
	int          inc_e;
	int          rs_C, cs_C;
	int          rs_S, cs_S;
	int          rs_U, cs_U;

	datatype = FLA_Obj_datatype( U );

	m_A      = FLA_Obj_vector_dim( d );
	m_U      = FLA_Obj_length( U );
	n_CS     = FLA_Obj_width( C );

	inc_d    = FLA_Obj_vector_inc( d );
	inc_e    = FLA_Obj_vector_inc( e );
	
	rs_C     = FLA_Obj_row_stride( C );
	cs_C     = FLA_Obj_col_stride( C );

	rs_S     = FLA_Obj_row_stride( S );
	cs_S     = FLA_Obj_col_stride( S );

	rs_U     = FLA_Obj_row_stride( U );
	cs_U     = FLA_Obj_col_stride( U );


	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*    buff_d = FLA_FLOAT_PTR( d );
			float*    buff_e = FLA_FLOAT_PTR( e );
			float*    buff_C = FLA_FLOAT_PTR( C );
			float*    buff_S = FLA_FLOAT_PTR( S );
			float*    buff_U = FLA_FLOAT_PTR( U );

			FLA_Tevd_v_ops_var4r( m_A,
			                      m_U,
			                      n_CS,
			                      0,
			                      0,
			                      buff_d, inc_d,
			                      buff_e, inc_e,
			                      buff_C, rs_C, cs_C,
			                      buff_S, rs_S, cs_S,
			                      buff_U, rs_U, cs_U );

			break;
		}

		case FLA_DOUBLE:
		{
			double*   buff_d = FLA_DOUBLE_PTR( d );
			double*   buff_e = FLA_DOUBLE_PTR( e );
			double*   buff_C = FLA_DOUBLE_PTR( C );
			double*   buff_S = FLA_DOUBLE_PTR( S );
			double*   buff_U = FLA_DOUBLE_PTR( U );

			FLA_Tevd_v_opd_var4r( m_A,
			                      m_U,
			                      n_CS,
			                      0,
			                      0,
			                      buff_d, inc_d,
			                      buff_e, inc_e,
			                      buff_C, rs_C, cs_C,
			                      buff_S, rs_S, cs_S,
			                      buff_U, rs_U, cs_U );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Tevd_v_ops_var4r( int       m_A,
                                int       m_U,
                                int       n_CS,
                                int       ij,
                                int       k_iter,
                                float*    buff_d, int inc_d, 
                                float*    buff_e, int inc_e,
                                float*    buff_C, int rs_C, int cs_C,
                                float*    buff_S, int rs_S, int cs_S,
                                float*    buff_U, int rs_U, int cs_U )
{
	return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Tevd_v_opd_var4r( int       m_A,
                                int       m_U,
                                int       n_CS,
                                int       ij,
                                int       k_iter,
                                double*   buff_d, int inc_d, 
                                double*   buff_e, int inc_e,
                                double*   buff_C, int rs_C, int cs_C,
                                double*   buff_S, int rs_S, int cs_S,
                                double*   buff_G, int rs_G, int cs_G )
{
	FLA_Error r_val;
	int       i, k;

	// Iterate from back to front until all that is left is a 2x2.
	for ( i = m_A - 1; i > 1; --i )
	{
		int m_ATL = i + 1;

		/*------------------------------------------------------------*/

		// Find an eigenvalue of ATL submatrix.
		r_val = FLA_Eigval_v_opd_var2( m_ATL,
		                               n_CS,
		                               buff_C, rs_C, cs_C,
		                               buff_S, rs_S, cs_S,
		                               buff_d, inc_d,
		                               buff_e, inc_e,
		                               &k );

		if ( r_val == FLA_FAILURE )
		{
			printf( "FLA_Tevd_v_opd_var4r: failed to converge (iteration %d)\n", i );
			FLA_Abort();
		}

#ifdef PRINTF
		printf( "FLA_Tevd_v_opd_var4r: found eig  %14.11f in col %3d (n=%d) after %3d iterations\n", buff_d[ r_val*inc_d ], r_val, m_ATL, k );
#endif

		// Apply the resulting Givens rotations to update the Givens matrix.
		FLA_Apply_G_rf_opd_var1b( k,
		                          m_U,
		                          m_ATL,
		                          ij,
		                          k_iter,
		                          buff_C, rs_C, cs_C,
		                          buff_S, rs_S, cs_S,
		                          buff_G, rs_G, cs_G );

		k_iter += k;
		//k_weight += i * k;

		// If r_val != i, then the eigenvalue converged elsewhere within A.
		// Therefore, we must recurse with subproblems.
		if ( r_val != i )
		{
#ifdef PRINTF
			printf( "FLA_Tevd_v_opd_var4r: Internal deflation in col %d, eig %d\n", r_val, i );
			printf( "FLA_Tevd_v_opd_var4r: alpha11         = %23.19e\n", buff_d[r_val*inc_d] );
			printf( "FLA_Tevd_v_opd_var4r: alpha21 alpha22 = %23.19e %23.19e\n", buff_e[r_val*inc_e], buff_d[(r_val+1)*inc_d] );
#endif

			//printf( "FLA_Tevd_v_opd_var4r: iter so far = %3d\n", k_total );

			// An eigenvalue converged somewhere within the diagonal (not at
			// either the end), so we have to recurse with two subproblems.
			{
				int       m_ATL = r_val;
				int       m_ABR = m_A - r_val - 1;
				int       ijT   = 0;
				int       ijB   = m_ATL + 1;
				double*   dTL   = buff_d + (0      )*inc_d;
				double*   eTL   = buff_e + (0      )*inc_e;
				double*   CT    = buff_C + (0      )*rs_C;
				double*   ST    = buff_S + (0      )*rs_S;
				double*   GL    = buff_G + (0      )*cs_G;
				double*   dBR   = buff_d + (m_ATL+1)*inc_d;
				double*   eBR   = buff_e + (m_ATL+1)*inc_e;
				double*   CB    = buff_C + (m_ATL+1)*rs_C;
				double*   SB    = buff_S + (m_ATL+1)*rs_S;
				double*   GR    = buff_G + (m_ATL+1)*cs_G;

				FLA_Tevd_v_opd_var4r( m_ATL,
				                      m_U,
				                      n_CS,
				                      ijT,
				                      k_iter,
				                      dTL, inc_d,
				                      eTL, inc_e,
				                      CT,  rs_C, cs_C,
				                      ST,  rs_S, cs_S,
				                      GL,  rs_G, cs_G );
				FLA_Tevd_v_opd_var4r( m_ABR,
				                      m_U,
				                      n_CS,
				                      ijB,
				                      k_iter,
				                      dBR, inc_d,
				                      eBR, inc_e,
				                      CB,  rs_C, cs_C,
				                      SB,  rs_S, cs_S,
				                      GR,  rs_G, cs_G );
				return FLA_SUCCESS;
			}
		}

		/*------------------------------------------------------------*/
	}

	//printf( "FLA_Tevd_v_opd_var4r: total iter  = %3d\n", k_total );

	// Skip 1x1 matrices (and submatrices) entirely.
	if ( m_A > 1 )
	{
		double gamma;
		double sigma;

		// Find the eigenvalue decomposition of the remaining (or only) 2x2
		// submatrix.
		FLA_Hevv_2x2_opd( buff_d + (0  )*inc_d,
		                  buff_e + (0  )*inc_e,
		                  buff_d + (1  )*inc_d,
		                  buff_d + (0  )*inc_d,
		                  buff_d + (1  )*inc_d,
		                  &gamma,
		                  &sigma );

		// Update the eigenvectors.
		FLA_Apply_G_mx2_opd( m_U,
		                     &gamma,
		                     &sigma,
		                     buff_G + (0  )*cs_G, rs_G,
		                     buff_G + (1  )*cs_G, rs_G );

#ifdef PRINTF
		printf( "FLA_Tevd_v_opd_var4r: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 1*inc_d ], 1, 2, 1 );
		printf( "FLA_Tevd_v_opd_var4r: found eig  %14.11f in col %3d (n=%d) after %3d iterations (via Hevv)\n", buff_d[ 0*inc_d ], 0, 2, 1 );
#endif
	}

	return FLA_SUCCESS;
}
